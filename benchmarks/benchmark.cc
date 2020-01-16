#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <cstdint>
#include <chrono>
#include <algorithm>
#include <functional>
#include <string.h>
#include <sstream>

//#include <execution>

#ifdef __linux__
#include "linux-perf-events.h"
#endif

::std::uint64_t NowNanos() {
  return ::std::chrono::duration_cast<::std::chrono::nanoseconds>(
             ::std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

using namespace std;

vector<uint64_t> GenerateRandom64Fast(size_t count, uint64_t start) {
  vector<uint64_t> result(count);
  uint64_t index = start;
  auto genrand = [&index]() {
    // mix64
    uint64_t x = index++;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9L;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebL;
    x = x ^ (x >> 31);
    return x;
  };
  generate(result.begin(), result.end(), ref(genrand));
  return result;
}

uint64_t compareCount = 0;

bool compare_string(string i1, string i2) {
#ifdef MEASURE
    compareCount++;
#endif
    return i1 < i2;
}

bool compare(uint64_t i1, uint64_t i2) {
#ifdef MEASURE
    compareCount++;
#endif
    // return (i1 & 0x0f0f0f0f00f0f0f0fULL) < (i2 & 0x0f0f0f0f00f0f0f0fULL);
    return i1 < i2;
}

size_t branchy_search(const int *source, size_t n, int target) {
  size_t lo = 0;
  size_t hi = n;
  while (lo < hi) {
    size_t m = (lo + hi) / 2;
    if (target < source[m]) {
      hi = m;
    } else if (target > source[m]) {
      lo = m + 1;
    } else {
      return m;
    }
  }
  return hi;
}

template<typename T>
size_t branchfree_search(const T *source, size_t n, T target, bool (*compare)(T i1, T i2)) {
    const T *base = source;
    while (n > 1) {
        size_t half = n >> 1;
        base = compare(base[half], target) ? &base[half] : base;
        n -= half;
    }
    return compare(base[0], target) + base - source;
}

// TODO check branchfree search has fixed number of comparisons
const int groups = 256;
const int oversampling = 64;

template<typename T>
void stable_counting_comparison_sort(vector<T>& data, bool (*compare)(T i1, T i2)) {

    // number of elements to sort
    int size = data.size();

    if (size <= groups * oversampling) {
        // too small
        stable_sort(data.begin(), data.end(), compare);
        return;
    }

#ifdef MEASURE
#ifdef __linux__
  vector<int> evts;
  evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
  evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
  evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
  evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
  LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
  vector<unsigned long long> results;
  results.resize(evts.size());
  unified.start();
#endif
  auto start1 = NowNanos();
#endif

    // the samples are sentinels to group the data
    vector<T> sample2(oversampling * groups);
    copy_n(data.begin(), oversampling * groups, sample2.begin());
    // no need to use stable sort here, as it's only used to group items
    sort(sample2.begin(), sample2.end(), compare);

    vector<T> sample(groups - 1);
    for(int i = 0; i < groups - 1; i++) {
        sample[i] = sample2[i * oversampling + oversampling];
    }

    // (as in the regular counting sort)
    int* counts = new int[groups]();

    // for each entry, we remember in which group it belongs to
    int* pos = new int[size];

#ifdef MEASURE
  auto time1 = NowNanos() - start1;
  cout << "Time 1 (ns/key): " << (time1 / size) << endl;
  cout << "Comparisons: " << compareCount << endl;
#ifdef __linux__
  unified.end(results);
  printf("Cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
    results[0]*1.0/size,
    results[1]*1.0/size ,
    results[1]*1.0/results[0],
    results[2]*1.0/size,
    results[3]*1.0/size);
    unified.start();
#endif
  auto start2 = NowNanos();
#endif

    // calculate the groups, using binary search, and the counts
    for (int i = 0; i < size; i++) {
        // int x = branchfree_search(&sample[0], groups, data[i], compare);
        // int x = lower_bound(sample.begin(), sample.end(), data[i], compare) - sample.begin();
        int x = branchfree_search(&sample[0], groups - 1, data[i], compare);
        /*
        if (x != y) {
            cout << " " << data[i] << " x=" << x << " y=" << y << endl;
        }
        */
        pos[i] = x;
        counts[x]++;
    }

#ifdef MEASURE
  auto time2 = NowNanos() - start2;
  cout << "Time 2 (ns/key): " << (time2 / size) << endl;
  cout << "Comparisons: " << compareCount << endl;
#ifdef __linux__
  unified.end(results);
  printf("Cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
    results[0]*1.0/size,
    results[1]*1.0/size ,
    results[1]*1.0/results[0],
    results[2]*1.0/size,
    results[3]*1.0/size);
    unified.start();
#endif
  auto start3 = NowNanos();
#endif

/*
    for (int i = 1; i < groups; i++) {
    cout << "group " << i << " " << counts[i] << " = " << sample[i] << endl;
    }
*/

    // sum up all the counts so we know start positions
    for (int i = 1; i < groups; i++) {
        counts[i] = counts[i - 1] + counts[i];
    }
    // we want to remember the sizes
    int* counts2 = new int[size];
    std::copy(counts, counts + groups, counts2);

    // place entries in the right bucket
    vector<T> d2(size);
    for (int i = 0; i < size; i++) {
        d2[--counts[pos[i]]] = data[i];
    }

#ifdef MEASURE
  auto time3 = NowNanos() - start3;
  cout << "Time 3 (ns/key): " << (time3 / size) << endl;
#ifdef __linux__
  unified.end(results);
  printf("Cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
    results[0]*1.0/size,
    results[1]*1.0/size ,
    results[1]*1.0/results[0],
    results[2]*1.0/size,
    results[3]*1.0/size);
    unified.start();
#endif
  auto start4 = NowNanos();
#endif

    // copy the data back to the source array
    copy_n(d2.begin(), size, data.begin());

    // free up
    delete[] pos;
    delete[] counts;
    // sort the buckets
    int start = 0;
    int total = 0;
    for (int i = 0; i < groups; i++) {
        int next = counts2[i];
        stable_sort(data.begin() + start, data.begin() + next, compare);
        total += next - start;
        start = next;
    }
    delete[] counts2;

#ifdef MEASURE
  auto time4 = NowNanos() - start4;
  cout << "Time 4 (ns/key): " << (time4 / size) << endl;
  cout << "Comparisons: " << compareCount << endl;
#ifdef __linux__
  unified.end(results);
  printf("Cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
    results[0]*1.0/size,
    results[1]*1.0/size ,
    results[1]*1.0/results[0],
    results[2]*1.0/size,
    results[3]*1.0/size);
    unified.start();
#endif
#endif

}

template<typename T>
void counting_comparison_sort(vector<T>& data, bool (*compare)(T i1, T i2)) {
    int size = data.size();
    if (size <= groups * oversampling) {
        sort(data.begin(), data.end(), compare);
        return;
    }
    vector<T> sample2(oversampling * groups);
    copy_n(data.begin(), oversampling * groups, sample2.begin());
    sort(sample2.begin(), sample2.end(), compare);
    vector<T> sample(groups - 1);
    for(int i = 0; i < groups - 1; i++) {
        sample[i] = sample2[i * oversampling + oversampling];
    }
    int* counts = new int[groups]();
    int* pos = new int[size];
    for (int i = 0; i < size; i++) {
        int x = branchfree_search(&sample[0], groups - 1, data[i], compare);
        pos[i] = x;
        counts[x]++;
    }
    for (int i = 1; i < groups; i++) {
        counts[i] = counts[i - 1] + counts[i];
    }
    int* counts2 = new int[size];
    std::copy(counts, counts + groups, counts2);
    vector<T> d2(size);
    for (int i = 0; i < size; i++) {
        d2[--counts[pos[i]]] = data[i];
    }
    copy_n(d2.begin(), size, data.begin());
    delete[] pos;
    delete[] counts;
    int start = 0;
    int total = 0;
    for (int i = 0; i < groups; i++) {
        int next = counts2[i];
        sort(data.begin() + start, data.begin() + next, compare);
        total += next - start;
        start = next;
    }
    delete[] counts2;
}

int main(int argc, char * argv[]) {

    int algorithm = 0;
    if (argc > 1) {
        if (strcmp(argv[1], "stable_sort") == 0) {
            algorithm = 0;
        } else if (strcmp(argv[1], "stable_count") == 0) {
            algorithm = 1;
        } else if (strcmp(argv[1], "sort") == 0) {
            algorithm = 2;
        } else if (strcmp(argv[1], "count") == 0) {
            algorithm = 3;
        } else {
            cout << "Unknown algorithm, supported: stable_sort, stable_count, sort, count: " << argv[1] << endl;
            return 2;
        }
    }
    int size = 10 * 1000000;
    if (argc > 2) {
        stringstream input_string(argv[2]);
        input_string >> size;
        size *= 1000000;
        if (input_string.fail()) {
            cerr << "Invalid number: " << argv[2];
            return 2;
        }
    }

// #define STRING_DATA
#ifdef STRING_DATA
    vector<uint64_t> v1 = GenerateRandom64Fast(size, 0);
    vector<string> v(size);
    for (int i = 0; i < size; i++) {
        v[i] = "x" + to_string(v1[i]);
    }
    #define COMPARE compare_string
#else
    vector<uint64_t> v = GenerateRandom64Fast(size, 0);
    #define COMPARE compare
    uint64_t sum = 0;
    for (auto x : v) {
        sum += x;
    }
#endif

    // sort(v.begin(), v.end(), greater<int>());


    cout << "Size (million keys): " << (size / 1000000) << endl;
    if (algorithm == 0) {
        cout << "Algorithm: stable_sort" << endl;
    } else if (algorithm == 1) {
        cout << "Algorithm: stable_count" << endl;
        cout << "Groups: " << groups << endl;
        cout << "Oversampling: " << oversampling << endl;
    } else if (algorithm == 2) {
        cout << "Algorithm: sort" << endl;
    } else if (algorithm == 3) {
        cout << "Algorithm: count" << endl;
        cout << "Groups: " << groups << endl;
        cout << "Oversampling: " << oversampling << endl;
    }
    cout << "Sorting..." << endl;

#ifdef MEASURE
#ifdef __linux__
  vector<int> evts;
  evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
  evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
  evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
  evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
  LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
  vector<unsigned long long> results;
  results.resize(evts.size());
  unified.start();
#endif
#endif

    auto start_time = NowNanos();
    if (algorithm == 0) {
        stable_sort(v.begin(), v.end(), COMPARE);
    } else if (algorithm == 1) {
        stable_counting_comparison_sort(v, COMPARE);
    } else if (algorithm == 2) {
        sort(v.begin(), v.end(), COMPARE);
    } else if (algorithm == 3) {
        counting_comparison_sort(v, COMPARE);
    }
    auto time = NowNanos() - start_time;
    cout << "Time (ns/key): " << (time / size) << endl;

#ifdef MEASURE
#ifdef __linux__
  unified.end(results);
  printf("Cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
    results[0]*1.0/size,
    results[1]*1.0/size ,
    results[1]*1.0/results[0],
    results[2]*1.0/size,
    results[3]*1.0/size);
#endif
#endif

    int duplicates = 0;
    for (int i = 1; i < size; i++) {
        if (v[i - 1] == v[i]) {
            duplicates++;
        }
        if (!COMPARE(v[i - 1], v[i])) {
            cout << v[i - 1] << " > " << v[i] << endl;
        }
    }
    if (duplicates > 0) {
        cout << "Duplicate keys: " << duplicates << endl;
    }
#ifdef STRING_DATA
#else
    uint64_t sum2 = 0;
    for (auto x : v) {
        sum2 += x;
    }
    if (sum != sum2) {
        cout << "Incorrect sum: " << sum << " != " << sum2 << endl;
    }
#endif

#ifdef MEASURE
    cout << "Comparisons: " << compareCount << endl;
#endif

/*
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <numberOfEntries>" << endl;
    return 1;
  }
*/
}

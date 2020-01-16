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
    compareCount++;
    return i1 < i2;
}

bool compare(uint64_t i1, uint64_t i2) {
    compareCount++;
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
// auto start1 = NowNanos();

    // the samples are sentinels to group the data
    vector<T> sample2(oversampling * groups);
    copy_n(data.begin(), oversampling * groups, sample2.begin());
    // no need to use stable sort here, as it's only used to group items
    sort(sample2.begin(), sample2.end(), compare);

    vector<T> sample(groups - 1);
    for(int i = 0; i < groups - 1; i++) {
        sample[i] = sample2[i * oversampling + oversampling];
    }


//    copy_n(data.begin(), groups, sample.begin());
//    sort(sample.begin(), sample.end(), compare);

    // (as in the regular counting sort)
    int* counts = new int[groups]();

    // for each entry, we remember in which group it belongs to
    int* pos = new int[size];
// auto time1 = NowNanos() - start1;
// cout << "time1 " << (time1 / size) << endl;
// cout << "comp " << compareCount << endl;

// auto start2 = NowNanos();
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
        // x = x >= 0 ? x : (-x - 1);
        pos[i] = x;
        counts[x]++;
    }
// auto time2 = NowNanos() - start2;
// cout << "time2 " << (time2 / size) << endl;
// cout << "comp " << compareCount << endl;

/*
    for (int i = 1; i < groups; i++) {
    cout << "group " << i << " " << counts[i] << " = " << sample[i] << endl;
    }
*/

// auto start3 = NowNanos();
    // sum up all the counts so we know start positions
    for (int i = 1; i < groups; i++) {
        counts[i] = counts[i - 1] + counts[i];
    }
    // we want to remember the sizes
    int* counts2 = new int[size];
    std::copy(counts, counts + groups, counts2);
// auto time3 = NowNanos() - start3;
// cout << "time3 " << (time3 / size) << endl;

// auto start4 = NowNanos();
    // place entries in the right bucket
    vector<T> d2(size);
    for (int i = 0; i < size; i++) {
        d2[--counts[pos[i]]] = data[i];
    }
// auto time4 = NowNanos() - start4;
// cout << "time4 " << (time4 / size) << endl;

// auto start5 = NowNanos();
    // copy the data back to the source array
    copy_n(d2.begin(), size, data.begin());

    // free up
    delete[] pos;
    delete[] counts;
// auto time5 = NowNanos() - start5;
// cout << "time5 " << (time5 / size) << endl;

// auto start6 = NowNanos();
    // sort the buckets
//    stable_sort(data.begin(), data.end(), compare);

    int start = 0;
    int total = 0;
    for (int i = 0; i < groups; i++) {
        int next = counts2[i];
        stable_sort(data.begin() + start, data.begin() + next, compare);
        total += next - start;
// cout << "group " << i << " " << (next - start) << endl;

        start = next;
    }

// auto time6 = NowNanos() - start6;
// cout << "time6 " << (time6 / size) << endl;
// cout << "total " << total << endl;

    delete[] counts2;
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
            cout << "Unknown algorithm, supported: stable_sort, stable_count, sort: " << argv[1] << endl;
            return 2;
        }
    }
    int size = 10 * 1000000;
    if (argc > 2) {
        stringstream input_string(argv[2]);
        input_string >> size;
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


    if (algorithm == 0) {
        cout << "algorithm: stable_sort" << endl;
    } else if (algorithm == 1) {
        cout << "algorithm: stable_count" << endl;
    } else if (algorithm == 2) {
        cout << "algorithm: sort" << endl;
    } else if (algorithm == 3) {
        cout << "algorithm: count" << endl;
    }
    cout << "size " << size << endl;
    cout << "groups " << groups << endl;
    cout << "start" << endl;
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
    cout << "time " << (time / size) << endl;

#ifdef STRING_DATA
#else
    uint64_t last = 0;
    int duplicates = 0;
    for (auto x : v) {
        if (x == last) {
            duplicates++;
        }
        if (x < last) {
            cout << last << " < " << x << endl;
        }
        last = x;
    }
    uint64_t sum2 = 0;
    for (auto x : v) {
        sum2 += x;
    }
    if (sum != sum2) {
        cout << "incorrect sum " << sum << " != " << sum2 << endl;
    }
    cout << "duplicates " << duplicates << endl;
#endif

    cout << "comp " << compareCount << endl;

/*
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <numberOfEntries>" << endl;
    return 1;
  }
*/
}

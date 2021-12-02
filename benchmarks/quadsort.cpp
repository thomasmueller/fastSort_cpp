/*
	Copyright (C) 2014-2020 Igor van den Hoven ivdhoven@gmail.com
*/

/*
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
	quadsort 1.1.1.1
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <algorithm>


typedef int CMPFUNC (const void *a, const void *b);

template <CMPFUNC cmp>
void quad_swap32(void *array, void *swap, size_t nmemb)
{
	size_t offset;

	int *pta = (int*)array;
	int *pts = (int*)swap;

	for (offset = 0 ; offset + 4 <= nmemb ; offset += 4)
	{
		if (cmp(&pta[0], &pta[1]) > 0)
		{
			pts[0] = pta[1];
			pts[1] = pta[0];
		}
		else
		{
			pts[0] = pta[0];
			pts[1] = pta[1];
		}

		if (cmp(&pta[2], &pta[3]) > 0)
		{
			pts[2] = pta[3];
			pts[3] = pta[2];
		}
		else
		{
			pts[2] = pta[2];
			pts[3] = pta[3];
		}

		if (cmp(&pts[1], &pts[2]) <= 0)
		{
			*pta++ = pts[0];
			*pta++ = pts[1];
			*pta++ = pts[2];
			*pta++ = pts[3];
		}
		else if (cmp(&pts[0], &pts[3]) > 0)
		{
			*pta++ = pts[2];
			*pta++ = pts[3];
			*pta++ = pts[0];
			*pta++ = pts[1];
		}
		else if (cmp(&pts[0], &pts[2]) > 0)
		{
			*pta++ = pts[2];
			*pta++ = pts[0];

			if (cmp(&pts[1], &pts[3]) > 0)
			{
				*pta++ = pts[3];
				*pta++ = pts[1];
			}
			else
			{
				*pta++ = pts[1];
				*pta++ = pts[3];
			}
		}
		else
		{
			*pta++ = pts[0];
			*pta++ = pts[2];

			if (cmp(&pts[1], &pts[3]) > 0)
			{
				*pta++ = pts[3];
				*pta++ = pts[1];
			}
			else
			{
				*pta++ = pts[1];
				*pta++ = pts[3];
			}
		}
	}

	switch (nmemb - offset)
	{
		case 0:
		case 1:
			return;
		case 2:
			if (cmp(&pta[0], &pta[1]) > 0)
			{
				pts[0] = pta[0];
				pta[0] = pta[1];
				pta[1] = pts[0];
			}
			return;
		case 3:
			if (cmp(&pta[0], &pta[1]) > 0)
			{
				pts[0] = pta[0];
				pta[0] = pta[1];
				pta[1] = pts[0];
			}
			if (cmp(&pta[1], &pta[2]) > 0)
			{
				pts[0] = pta[1];
				pta[1] = pta[2];
				pta[2] = pts[0];
			}
			if (cmp(&pta[0], &pta[1]) > 0)
			{
				pts[0] = pta[0];
				pta[0] = pta[1];
				pta[1] = pts[0];
			}
			return;
		default:
			assert(nmemb < 1 && nmemb > 3);
	}
}

template <CMPFUNC cmp>
void quad_swap64(void *array, void *swap, size_t nmemb)
{
	size_t offset;
	long long *pts, *pta;

	pta = (long long*)array;
	pts = (long long*)swap;

	for (offset = 0 ; offset + 4 <= nmemb ; offset += 4)
	{
		if (cmp(&pta[0], &pta[1]) > 0)
		{
			pts[0] = pta[1];
			pts[1] = pta[0];
		}
		else
		{
			pts[0] = pta[0];
			pts[1] = pta[1];
		}

		if (cmp(&pta[2], &pta[3]) > 0)
		{
			pts[2] = pta[3];
			pts[3] = pta[2];
		}
		else
		{
			pts[2] = pta[2];
			pts[3] = pta[3];
		}

		if (cmp(&pts[1], &pts[2]) <= 0)
		{
			*pta++ = pts[0];
			*pta++ = pts[1];
			*pta++ = pts[2];
			*pta++ = pts[3];
		}
		else if (cmp(&pts[0], &pts[3]) > 0)
		{
			*pta++ = pts[2];
			*pta++ = pts[3];
			*pta++ = pts[0];
			*pta++ = pts[1];
		}
		else if (cmp(&pts[0], &pts[2]) > 0)
		{
			*pta++ = pts[2];
			*pta++ = pts[0];

			if (cmp(&pts[1], &pts[3]) > 0)
			{
				*pta++ = pts[3];
				*pta++ = pts[1];
			}
			else
			{
				*pta++ = pts[1];
				*pta++ = pts[3];
			}
		}
		else
		{
			*pta++ = pts[0];
			*pta++ = pts[2];

			if (cmp(&pts[1], &pts[3]) > 0)
			{
				*pta++ = pts[3];
				*pta++ = pts[1];
			}
			else
			{
				*pta++ = pts[1];
				*pta++ = pts[3];
			}
		}
	}

	switch (nmemb - offset)
	{
		case 0:
		case 1:
			return;
		case 2:
			if (cmp(&pta[0], &pta[1]) > 0)
			{
				pts[0] = pta[0];
				pta[0] = pta[1];
				pta[1] = pts[0];
			}
			return;
		case 3:
			if (cmp(&pta[0], &pta[1]) > 0)
			{
				pts[0] = pta[0];
				pta[0] = pta[1];
				pta[1] = pts[0];
			}
			if (cmp(&pta[1], &pta[2]) > 0)
			{
				pts[0] = pta[1];
				pta[1] = pta[2];
				pta[2] = pts[0];
			}
			if (cmp(&pta[0], &pta[1]) > 0)
			{
				pts[0] = pta[0];
				pta[0] = pta[1];
				pta[1] = pts[0];
			}
			return;
		default:
			assert(nmemb < 1 && nmemb > 3);
	}
}

template <CMPFUNC cmp>
void quad_sort32(void *array, void *swap, size_t nmemb)
{
	size_t offset, block = 4;
	int *pta, *pts, *c, *c_max, *d, *d_max, *end;

	end = (int*)array;
	end += nmemb;

	while (block < nmemb)
	{
		offset = 0;

		while (offset + block < nmemb)
		{
			pta = (decltype(pta))array;
			pta += offset;

			d_max = pta + block;

			if (cmp(d_max - 1, d_max) <= 0)
			{
				if (offset + block * 3 < nmemb)
				{
					d_max = pta + block * 3;

					if (cmp(d_max - 1, d_max) <= 0)
					{
						d_max = pta + block * 2;

						if (cmp(d_max - 1, d_max) <= 0)
						{
							offset += block * 4;
							continue;
						}
						pts = (decltype(pts))swap;

						c = pta;
						c_max = pta + block * 2;
						d = c_max;
						d_max = offset + block * 4 <= nmemb ? d + block * 2 : end;

						while (c < c_max)
							*pts++ = *c++;

						while (d < d_max)
							*pts++ = *d++;

						goto step3;
					}
					pts = (int*)swap;

					c = pta;
					c_max = pta + block * 2;

					while (c < c_max)
						*pts++ = *c++;

					goto step2;
				}
				else if (offset + block * 2 < nmemb)
				{
					d_max = pta + block * 2;

					if (cmp(d_max - 1, d_max) <= 0)
					{
						offset += block * 4;
						continue;
					}
					pts = (int*)swap;

					c = pta;
					c_max = pta + block * 2;

					while (c < c_max)
						*pts++ = *c++;

					goto step2;
				}
				else
				{
					offset += block * 4;
					continue;
				}
			}

			// step1:

			pts = (int*)swap;

			c = pta;
			c_max = pta + block;

			d = c_max;
			d_max = offset + block * 2 <= nmemb ? d + block : end;

			if (cmp(c_max - 1, d_max - 1) <= 0)
			{
				while (c < c_max)
				{
					while (cmp(c, d) > 0)
					{
						*pts++ = *d++;
					}
					*pts++ = *c++;
				}
				while (d < d_max)
					*pts++ = *d++;
			}
			else if (cmp(c, d_max - 1) > 0)
			{
				while (d < d_max)
					*pts++ = *d++;

				while (c < c_max)
					*pts++ = *c++;
			}
			else
			{
				while (d < d_max)
				{
					while (cmp(c, d) <= 0)
					{
						*pts++ = *c++;
					}
					*pts++ = *d++;
				}

				while (c < c_max)
				{
					*pts++ = *c++;
				}
			}

			step2:

			if (offset + block * 2 < nmemb)
			{
				c = pta + block * 2;

				if (offset + block * 3 < nmemb)
				{
					c_max = c + block;
					d = c_max;
					d_max = offset + block * 4 <= nmemb ? d + block : end;

					if (cmp(c_max - 1, d_max - 1) <= 0)
					{
						while (c < c_max)
						{
							while (cmp(c, d) > 0)
							{
								*pts++ = *d++;
							}
							*pts++ = *c++;
						}
						while (d < d_max)
							*pts++ = *d++;
					}
					else if (cmp(c, d_max - 1) > 0)
					{
						while (d < d_max)
							*pts++ = *d++;
						while (c < c_max)
							*pts++ = *c++;
					}
					else
					{
						while (d < d_max)
						{
							while (cmp(c, d) <= 0)
							{
								*pts++ = *c++;
							}
							*pts++ = *d++;
						}
						while (c < c_max)
							*pts++ = *c++;
					}
				}
				else
				{
					while (c < end)
						*pts++ = *c++;
				}
			}

			step3:

			pts = (int*)swap;

			c = pts;

			if (offset + block * 2 < nmemb)
			{
				c_max = c + block * 2;

				d = c_max;
				d_max = offset + block * 4 <= nmemb ? d + block * 2 : pts + nmemb - offset;

				if (cmp(c_max - 1, d_max - 1) <= 0)
				{
					while (c < c_max)
					{
						while (cmp(c, d) > 0)
						{
							*pta++ = *d++;
						}
						*pta++ = *c++;
					}
					while (d < d_max)
						*pta++ = *d++;
				}
				else if (cmp(c, d_max - 1) > 0)
				{
					while (d < d_max)
						*pta++ = *d++;
					while (c < c_max)
						*pta++ = *c++;
				}
				else
				{
					while (d < d_max)
					{
						while (cmp(d, c) > 0)
						{
							*pta++ = *c++;
						}
						*pta++ = *d++;
					}
					while (c < c_max)
						*pta++ = *c++;
				}
			}
			else
			{
				d_max = pts + nmemb - offset;

				while (c < d_max)
					*pta++ = *c++;
			}
			offset += block * 4;
		}
		block *= 4;
	}
}

template <CMPFUNC cmp>
void quad_sort64(void *array, void *swap, size_t nmemb)
{
	size_t offset, block = 4;
	long long *pta, *pts, *c, *c_max, *d, *d_max, *end;

	end = (decltype(end))array;
	end += nmemb;

	while (block < nmemb)
	{
		offset = 0;

		while (offset + block < nmemb)
		{
			pta = (decltype(pta))array;
			pta += offset;

			d_max = pta + block;

			if (cmp(d_max - 1, d_max) <= 0)
			{
				if (offset + block * 3 < nmemb)
				{
					d_max = pta + block * 3;

					if (cmp(d_max - 1, d_max) <= 0)
					{
						d_max = pta + block * 2;

						if (cmp(d_max - 1, d_max) <= 0)
						{
							offset += block * 4;
							continue;
						}
						pts = (decltype(pts))swap;

						c = pta;
						c_max = pta + block * 2;
						d = c_max;
						d_max = offset + block * 4 <= nmemb ? d + block * 2 : end;

						while (c < c_max)
							*pts++ = *c++;

						while (d < d_max)
							*pts++ = *d++;

						goto step3;
					}
					pts = (decltype(pts))swap;

					c = pta;
					c_max = pta + block * 2;

					while (c < c_max)
						*pts++ = *c++;

					goto step2;
				}
				else if (offset + block * 2 < nmemb)
				{
					d_max = pta + block * 2;

					if (cmp(d_max - 1, d_max) <= 0)
					{
						offset += block * 4;
						continue;
					}
					pts = (decltype(pts))swap;

					c = pta;
					c_max = pta + block * 2;

					while (c < c_max)
						*pts++ = *c++;

					goto step2;
				}
				else
				{
					offset += block * 4;
					continue;
				}
			}

			// step1:

			pts = (decltype(pts))swap;

			c = pta;
			c_max = pta + block;

			d = c_max;
			d_max = offset + block * 2 <= nmemb ? d + block : end;

			if (cmp(c_max - 1, d_max - 1) <= 0)
			{
				while (c < c_max)
				{
					while (cmp(c, d) > 0)
					{
						*pts++ = *d++;
					}
					*pts++ = *c++;
				}
				while (d < d_max)
					*pts++ = *d++;
			}
			else if (cmp(c, d_max - 1) > 0)
			{
				while (d < d_max)
					*pts++ = *d++;

				while (c < c_max)
					*pts++ = *c++;
			}
			else
			{
				while (d < d_max)
				{
					while (cmp(c, d) <= 0)
					{
						*pts++ = *c++;
					}
					*pts++ = *d++;
				}

				while (c < c_max)
				{
					*pts++ = *c++;
				}
			}

			step2:

			if (offset + block * 2 < nmemb)
			{
				c = pta + block * 2;

				if (offset + block * 3 < nmemb)
				{
					c_max = c + block;
					d = c_max;
					d_max = offset + block * 4 <= nmemb ? d + block : end;

					if (cmp(c_max - 1, d_max - 1) <= 0)
					{
						while (c < c_max)
						{
							while (cmp(c, d) > 0)
							{
								*pts++ = *d++;
							}
							*pts++ = *c++;
						}
						while (d < d_max)
							*pts++ = *d++;
					}
					else if (cmp(c, d_max - 1) > 0)
					{
						while (d < d_max)
							*pts++ = *d++;
						while (c < c_max)
							*pts++ = *c++;
					}
					else
					{
						while (d < d_max)
						{
							while (cmp(c, d) <= 0)
							{
								*pts++ = *c++;
							}
							*pts++ = *d++;
						}
						while (c < c_max)
							*pts++ = *c++;
					}
				}
				else
				{
					while (c < end)
						*pts++ = *c++;
				}
			}

			step3:

			pts = (long long*)swap;

			c = pts;

			if (offset + block * 2 < nmemb)
			{
				c_max = c + block * 2;

				d = c_max;
				d_max = offset + block * 4 <= nmemb ? d + block * 2 : pts + nmemb - offset;

				if (cmp(c_max - 1, d_max - 1) <= 0)
				{
					while (c < c_max)
					{
						while (cmp(c, d) > 0)
						{
							*pta++ = *d++;
						}
						*pta++ = *c++;
					}
					while (d < d_max)
						*pta++ = *d++;
				}
				else if (cmp(c, d_max - 1) > 0)
				{
					while (d < d_max)
						*pta++ = *d++;
					while (c < c_max)
						*pta++ = *c++;
				}
				else
				{
					while (d < d_max)
					{
						while (cmp(d, c) > 0)
						{
							*pta++ = *c++;
						}
						*pta++ = *d++;
					}
					while (c < c_max)
						*pta++ = *c++;
				}
			}
			else
			{
				d_max = pts + nmemb - offset;

				while (c < d_max)
					*pta++ = *c++;
			}
			offset += block * 4;
		}
		block *= 4;
	}
}


template <CMPFUNC cmp>
void quadsort(void *array, size_t nmemb, size_t size)
{
	void *swap;

	swap = malloc(nmemb * size);

	if (size == sizeof(int))
	{
		quad_swap32<cmp>(array, swap, nmemb);

		quad_sort32<cmp>(array, swap, nmemb);
	}
	else if (size == sizeof(long long))
	{
		quad_swap64<cmp>(array, swap, nmemb);

		quad_sort64<cmp>(array, swap, nmemb);
	}
	else
	{
		assert(size == 4 || size == 8);
	}

	free(swap);
}

static inline int cmp_int(const void * a, const void * b)
{
	return *(int *) a - *(int *) b;
}

static inline int cmp_str(const void * a, const void * b)
{
	return strcmp(*(const char **) a, *(const char **) b);
}

static inline int cmp_float(const void * a, const void * b)
{
	return *(float *) a - *(float *) b;
}

// benchmarking utilities

long long utime()
{
	struct timeval now_time;

	gettimeofday(&now_time, NULL);

	return now_time.tv_sec * 1000000LL + now_time.tv_usec;
}

void test_quad(int *z_array, int *r_array, int max, const char *desc)
{
	long long start, end, cnt;

	memcpy(z_array, r_array, max * sizeof(int));

	if (max <= 10) printf("\e[1;31m%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", z_array[0], z_array[1], z_array[2], z_array[3], z_array[4], z_array[5], z_array[6], z_array[7], z_array[8], z_array[9]);

	start = utime();

	quadsort<cmp_int>(z_array, max, sizeof(int));

	end = utime();

	if (max <= 10) printf("\e[1;32m%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", z_array[0], z_array[1], z_array[2], z_array[3], z_array[4], z_array[5], z_array[6], z_array[7], z_array[8], z_array[9]);

	printf("\e[0m        Quad Sort: sorted %d elements in %f seconds. (%s)\n", max, (end - start) / 1000000.0, desc);

	for (cnt = 1 ; cnt < max ; cnt++) if (z_array[cnt - 1] > z_array[cnt]) {printf("        Quad Sort: not properly sorted at index %lld. (%d vs %d\n", cnt, z_array[cnt - 1], z_array[cnt]);break;}
}

void test_quick(int *z_array, int *r_array, int max, const char *desc)
{
	long long start, end, cnt;

	memcpy(z_array, r_array, max * sizeof(int));

	if (max <= 10) printf("\e[1;31m%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", z_array[0], z_array[1], z_array[2], z_array[3], z_array[4], z_array[5], z_array[6], z_array[7], z_array[8], z_array[9]);

	start = utime();

	std::sort(z_array, z_array+max);

	end = utime();

	if (max <= 10) printf("\e[1;32m%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", z_array[0], z_array[1], z_array[2], z_array[3], z_array[4], z_array[5], z_array[6], z_array[7], z_array[8], z_array[9]);

	printf("\e[0m       Quick Sort: sorted %d elements in %f seconds. (%s)\n", max, (end - start) / 1000000.0, desc);

	for (cnt = 1 ; cnt < max ; cnt++) if (z_array[cnt - 1] > z_array[cnt]) {printf("       Quick Sort: not properly sorted at index %lld. (%d vs %d\n", cnt, z_array[cnt - 1], z_array[cnt]);break;}
}

const int groups = 256;
const int oversampling = 64;

template <typename T, CMPFUNC cmp>
size_t branchfree_search(const T *source, size_t n, T target) {
    const T *base = source;
    while (n > 1) {
        size_t half = n >> 1;
        base = cmp(&base[half], &target) <= 0 ? &base[half] : base;
        n -= half;
    }
    return (cmp(&base[0], &target) > 0 ? 1 : 0) + base - source;
}

template <CMPFUNC cmp>
void sample_sort(int *array, size_t size) {
    if (size <= groups * oversampling) {
        std::sort(array, array + size);
        return;
    }
    int* sample2 = new int[oversampling * groups];
    memcpy(array, sample2, oversampling * groups * sizeof(int));
    std::sort(sample2, sample2 + oversampling * groups);
    delete[] sample2;

    int* sample = new int[groups - 1];
    for(int i = 0; i < groups - 1; i++) {
        sample[i] = sample2[i * oversampling + oversampling];
    }

    int* counts = new int[groups]();
    int* pos = new int[size];
    for (size_t i = 0; i < size; i++) {
        int x = branchfree_search<int, cmp>(&sample[0], groups - 1, array[i]);
        pos[i] = x;
        counts[x]++;
    }
    for (int i = 1; i < groups; i++) {
        counts[i] = counts[i - 1] + counts[i];
    }
    int* counts2 = new int[size];
    std::copy(counts, counts + groups, counts2);

    int* d2 = new int[size];
    for (size_t i = 0; i < size; i++) {
        d2[--counts[pos[i]]] = array[i];
    }
    memcpy(d2, array, size * sizeof(int));

    delete[] d2;
    delete[] pos;
    delete[] counts;
    delete[] sample;

    int start = 0;
    int total = 0;
    for (int i = 0; i < groups; i++) {
        int next = counts2[i];
        std::sort(array + start, array + next);
        total += next - start;
        start = next;
    }
    delete[] counts2;
}

void test_sample(int *z_array, int *r_array, int max, const char *desc) {
	long long start, end, cnt;
	memcpy(z_array, r_array, max * sizeof(int));
	if (max <= 10) printf("\e[1;31m%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", z_array[0], z_array[1], z_array[2], z_array[3], z_array[4], z_array[5], z_array[6], z_array[7], z_array[8], z_array[9]);
	start = utime();
	sample_sort<cmp_int>(z_array, max);
	// std::sort(z_array, z_array+max);
	end = utime();
	if (max <= 10) printf("\e[1;32m%10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", z_array[0], z_array[1], z_array[2], z_array[3], z_array[4], z_array[5], z_array[6], z_array[7], z_array[8], z_array[9]);
	printf("\e[0m       Sample Sort: sorted %d elements in %f seconds. (%s)\n", max, (end - start) / 1000000.0, desc);
	for (cnt = 1 ; cnt < max ; cnt++) if (z_array[cnt - 1] > z_array[cnt]) {printf("       Sample Sort: not properly sorted at index %lld. (%d vs %d\n", cnt, z_array[cnt - 1], z_array[cnt]);break;}
}


int main(int argc, char **argv)
{
	static int max = 10000000;
	int cnt, rnd;
	rnd = 1;
	srand(rnd);
	int *z_array = (int*)malloc(max * sizeof(int));
	int *r_array = (int*)malloc(max * sizeof(int));

	srand(rnd);

	for (cnt = 0 ; cnt < max ; cnt++)
	{
		r_array[cnt] = rand();
	}

	test_quad(z_array, r_array, max, "random order");
	test_quick(z_array, r_array, max, "random order");
	test_sample(z_array, r_array, max, "random order");
	printf("\n");

	for (cnt = 0 ; cnt < max ; cnt++)
	{
		r_array[cnt] = cnt;
	}

	test_quad(z_array, r_array, max, "forward order");
	test_quick(z_array, r_array, max, "forward order");
	test_sample(z_array, r_array, max, "forward order");
	printf("\n");

	for (cnt = 0 ; cnt < max ; cnt++)
	{
		r_array[cnt] = max - cnt;
	}

	test_quad(z_array, r_array, max, "reverse order");
	test_quick(z_array, r_array, max, "reverse order");
	test_sample(z_array, r_array, max, "reverse order");
	printf("\n");

	srand(rnd);

	for (cnt = 0 ; cnt < max * 3 / 4 ; cnt++)
	{
		r_array[cnt] = cnt;
	}
	for (cnt = max * 3 / 4 ; cnt < max ; cnt++)
	{
		r_array[cnt] = rand();
	}

	test_quad(z_array, r_array, max, "random tail");
	test_quick(z_array, r_array, max, "random tail");
	test_sample(z_array, r_array, max, "random tail");
	printf("\n");

	free(z_array);
	free(r_array);

	return 0;
}

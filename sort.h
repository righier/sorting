#pragma once

#include <algorithm>
#include <cmath>

#include "utils.h"
#include "ssssort.h"

using namespace std;

class Sort {
public:
	virtual u64 sort(std::vector<u64> &v) = 0;
	virtual ~Sort() {}
};

static inline void step(u64 &prev, u64 x) {
	prev = (prev >> (prev & 1)) + x;
}

template<typename T>
static inline void eval(const T begin, const T end, u64 &output) {
	for (auto it = begin; it < end; it++) {
		step(output, *it);
	}
}

static inline void eval(u64 *__restrict v, int begin, int end, u64 &output) {
	for (int i = begin; i < end; i++) {
		step(output, v[i]);
	}
}

static inline u64 eval(vector<u64> &v) {
	u64 output = 0;
	eval(v.data(), 0, v.size(), output);
	return output;
}

static inline void insertionSort(u64 *__restrict__ v, int begin, int end) {
	for (int i = begin+1; i < end; i++) {
		u64 val = v[i];
		int j = i-1;
		for (; j >= begin && v[j] > val; j--) {
			v[j+1] = v[j];
		}
		v[j+1] = val;
	}
}

static inline void randomPivot(u64 *__restrict__ v, int begin, int end) {
	static Random rng(69);
	int pos = rng.getInt(begin, end);
	std::swap(v[pos], v[begin]);
}

static inline int partition(u64 *__restrict__ v, int begin, int end, u64 pivotVal) {
	while (begin < end) {
		if (v[begin] > pivotVal) {
			std::swap(v[begin], v[--end]);
		} else {
			++begin;
		}
	}
	return end;
}

class GccSort: public Sort {
public:
	u64 sort(std::vector<u64> &v) {
		std::sort(v.begin(), v.end());
		return eval(v);
	}
};

class SSort: public Sort {
public:
	u64 sort(std::vector<u64> &v) {
		ssssort::ssssort(v.begin(), v.end());
		return eval(v);
	}
};

class InsertionSort: public Sort {
public:
	u64 sort(std::vector<u64> &v) {
		insertionSort(v.data(), 0, v.size());
		return eval(v);
	}

};

class QuickSort: public Sort {

	void sortRec(u64 *__restrict__ v, int begin, int end, u64 &output) const {
		int size = end - begin;

		if (size == 1) {

			step(output, v[begin]);

		} else if (size > 1) {

			randomPivot(v, begin, end);
			u64 pivotVal = v[begin];
			int rightBegin = partition(v, begin+1, end, pivotVal);

			sortRec(v, begin+1, rightBegin, output);

			step(output, pivotVal);

			sortRec(v, rightBegin, end, output);

		}

	}

public:
	u64 sort(std::vector<u64> &v) {
		u64 output = 0;
		sortRec(v.data(), 0, v.size(), output);
		return output;
	}

};

class HybridSort: public Sort {
public:
	int threshold;

	void sortRec(u64 *__restrict__ v, int begin, int end, u64 &output) const {
		if (end - begin < threshold) {
			insertionSort(v, begin, end);
			eval(v, begin, end, output);
		} else {
			randomPivot(v, begin, end);
			u64 pivotVal = v[begin];
			int rightBegin = partition(v, begin+1, end, pivotVal);

			sortRec(v, begin+1, rightBegin, output);

			step(output, pivotVal);

			sortRec(v, rightBegin, end, output);
		}
	}

	HybridSort(int threshold = 50): threshold(threshold) { }

	u64 sort(std::vector<u64> &v) {
		u64 output = 0;
		sortRec(v.data(), 0, v.size(), output);
		return output;
	}
};

class BucketSort: public Sort {
	HybridSort hybrid;
	int numBuckets;
	int threshold;

	void sortRec(u64 *__restrict__ v, std::vector<std::vector<u64>> &buckets, int startBucket, int N, u64 &output) {
		if (N < threshold) {

			hybrid.sortRec(v, 0, N, output);

		} else {
			int endBucket = startBucket + numBuckets;
			if (endBucket > int(buckets.size())) {
				cout << "RESIZE" << endl;
				buckets.resize(endBucket);
			}

			u64 minVal = v[0];
			u64 maxVal = minVal;

			for (int i = 0; i < N; i++) {
				u64 x = v[i];
				minVal = min(minVal, x);
				maxVal = max(maxVal, x);
			}

			double irange = (double(maxVal - minVal) + 1.0) / double(numBuckets);


			// vector<vector<u64>> buckets(numBuckets);

			for (int i = 0; i < N; i++) {
				u64 x = v[i];
				int bucketId = int(double(x - minVal) / irange);
				buckets[startBucket + bucketId].push_back(x);
			}

			for (int i = startBucket; i < endBucket; i++) {
				auto &bucket = buckets[i];
				sortRec(bucket.data(), buckets, endBucket, bucket.size(), output);
				bucket.clear();
			}

			// for (auto &bucket: buckets) {
			// 	sortRec(bucket, output);
			// }

		}
	}

public:
	BucketSort(int numBuckets = 256, int threshold = 1000): numBuckets(numBuckets), threshold(threshold) { }

	u64 sort(std::vector<u64> &v) {
		u64 output = 0;

		int expectedRec = std::log(v.size()) / std::log(numBuckets);
		vector<vector<u64>> buckets(numBuckets * expectedRec);
		sortRec(v.data(), buckets, 0, v.size(), output);
		return output;
	}
};

class AdaptiveBucketSort: public Sort {
	HybridSort hybrid;
	int threshold;

	void sortRec(u64 *__restrict__ v, int N, u64 &output) {
		if (N < threshold) {

			hybrid.sortRec(v, 0, N, output);

		} else {

			int numBuckets = std::sqrt(N);

			u64 minVal = v[0];
			u64 maxVal = minVal;

			for (int i = 0; i < N; i++) {
				u64 x = v[i];
				minVal = min(minVal, x);
				maxVal = max(maxVal, x);
			}

			double irange = (double(maxVal - minVal) + 1.0) / double(numBuckets);


			vector<vector<u64>> buckets(numBuckets);

			for (int i = 0; i < N; i++) {
				u64 x = v[i];
				int bucketId = int(double(x - minVal) / irange);
				buckets[bucketId].push_back(x);
			}

			for (auto &bucket: buckets) {
				sortRec(bucket.data(), bucket.size(), output);
			}

		}
	}

public:
	AdaptiveBucketSort(int threshold = 1000): threshold(threshold) { }

	u64 sort(std::vector<u64> &v) {
		u64 output = 0;
		sortRec(v.data(), v.size(), output);
		return output;
	}
};

class BucketSort2: public Sort {
	HybridSort hybrid;
	int numBuckets;
	int threshold;

	void sortRec(u64 *__restrict__ v, u64 *__restrict__ v2, int begin, int end, vector<int> &sizes, int sizeStart, u64 &output) {
		int N = end - begin;
		if (N < threshold) {

			hybrid.sortRec(v, begin, end, output);

		} else {
			int sizeEnd = sizeStart + numBuckets;
			if ((int)sizes.size() < sizeEnd) {
				sizes.resize(sizeEnd);
			}

			// find min and max element
			u64 minVal = v[begin];
			u64 maxVal = minVal;
			for (int i = begin+1; i < end; i++) {
				u64 x = v[i];
				minVal = min(minVal, x);
				maxVal = max(maxVal, x);
			}

			double irange = (double(maxVal - minVal) + 1.0) / double(numBuckets);

			// count elements per bucket
			for (int i = begin; i < end; i++) {
				u64 x = v[i];
				int bucketId = int(double(x - minVal) / irange);
				sizes[sizeStart + bucketId]++;
			}

			// prefix sum
			int prev = begin;
			for (int i = sizeStart; i < sizeEnd; i++) {
				int x = sizes[i];
				sizes[i] = prev;
				prev += x;
			}

			// distribute elements
			for (int i = begin; i < end; i++) {
				u64 x = v[i];
				int bucketId = int(double(x - minVal) / irange);
				v2[sizes[sizeStart + bucketId]++] = x;
			}

			int newBegin = begin;
			for (int i = sizeStart; i < sizeEnd; i++) {
				int newEnd = sizes[i];
				sortRec(v2, v, newBegin, newEnd, sizes, sizeEnd, output);
				newBegin = newEnd;
			}

			std::fill(sizes.begin()+sizeStart, sizes.begin()+sizeEnd, 0);
		}
	}

public:
	BucketSort2(int numBuckets = 256, int threshold = 1000): numBuckets(numBuckets), threshold(threshold) { }

	u64 sort(std::vector<u64> &v) {
		u64 output = 0;
		std::vector<u64> v2(v.size());

		int expectedRec = std::log(v.size()) / std::log(numBuckets);
		std::vector<int> sizes(expectedRec);
		sortRec(v.data(), v2.data(), 0, v.size(), sizes, 0, output);
		return output;
	}
};


class BucketSort3: public Sort {
	HybridSort hybrid;
	static constexpr int logBuckets = 8;
	static constexpr int numBuckets = 1<<logBuckets;
	static constexpr int maskBucket = numBuckets - 1;
	int threshold;

	static inline int getId(int val, int shift) {
		return (val >> shift) & maskBucket;
	}

	void sortRec(u64 *__restrict__ v, u64 *__restrict__ v2, int begin, int end, vector<int> &sizes, int sizeStart, int shift, u64 &output) const {
		int N = end - begin;
		if (N < threshold) {

			hybrid.sortRec(v, begin, end, output);

		} else {
			int sizeEnd = sizeStart + numBuckets;
			if ((int)sizes.size() < sizeEnd) {
				sizes.resize(sizeEnd);
			}

			// count elements per bucket
			for (int i = begin; i < end; i++) {
				u64 x = v[i];
				int bucketId = getId(x, shift);
				sizes[sizeStart + bucketId]++;
			}

			// prefix sum
			int prev = begin;
			for (int i = sizeStart; i < sizeEnd; i++) {
				int x = sizes[i];
				sizes[i] = prev;
				prev += x;
			}

			// distribute elements
			for (int i = begin; i < end; i++) {
				u64 x = v[i];
				int bucketId = getId(x, shift);
				v2[sizes[sizeStart + bucketId]++] = x;
			}

			int newBegin = begin;
			for (int i = sizeStart; i < sizeEnd; i++) {
				int newEnd = sizes[i];
				sortRec(v2, v, newBegin, newEnd, sizes, sizeEnd, shift - 8, output);
				newBegin = newEnd;
			}

			std::fill(sizes.begin()+sizeStart, sizes.begin()+sizeEnd, 0);
		}
	}

public:
	BucketSort3(int threshold = 1000): threshold(threshold) { }

	u64 sort(std::vector<u64> &v) {
		u64 output = 0;
		std::vector<u64> v2(v.size());

		int shift = 31 - logBuckets;
		int expectedRec = std::log(v.size()) / std::log(numBuckets);
		std::vector<int> sizes(expectedRec);
		sortRec(v.data(), v2.data(), 0, v.size(), sizes, 0, shift, output);
		return output;
	}
};


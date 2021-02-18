#pragma once

#include "utils.h"
#include "sort.h"

typedef vector<u64>::iterator It;

constexpr int logBuckets = 8;
constexpr int numBuckets = 1 << logBuckets;
constexpr int numPivots = numBuckets - 1;

void sample(It begin, It end, size_t sample_size) {
	static Random rnd(31);

	assert(begin <= end);

	size_t max = end - begin;

	for (size_t i = 0; i < sample_size; i++) {
		size_t j = rnd.getULong(max);
		std::swap(*(begin+j), *(--end));
	}
}

void buildTree(const u64 *begin, const u64 *end, u64 *splitters, size_t pos) {
	const u64 *mid = begin + (end - begin) / 2;
	splitters[pos] = *mid;
	if (2 * pos < numPivots) {
		
	}
}

class SampleSort: public Sort {

public:
	u64 sort(std::vector<u64> &v) {

		return v.size();
	}

};
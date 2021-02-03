#pragma once

#include <vector>
#include <string>
#include <iostream>

#include "random.h"

typedef uint64_t u64;
typedef int64_t i64;
typedef uint32_t u32;
typedef uint8_t u8;

#define UNUSED(x) (void)(x)

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& v) {
    for (auto const &x : v) {
    	out << x << " ";
    }
    return out;
}

std::vector<u64> randomVector(Random &rng, int N) {
	std::vector<u64> v;
	v.reserve(N);
	for (int i = 0; i < N; i++) {
		v.push_back(rng.getPInt());
	}
	return v;
}

bool isSorted(const std::vector<u64> &v) {
	for (size_t i = 1; i < v.size(); i++) {
		if (v[i-1] > v[i]) {
			return false;
		}
	}
	return true;
}
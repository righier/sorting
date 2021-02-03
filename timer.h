#pragma once

#include <chrono>

template <typename Clock = std::chrono::high_resolution_clock>
class Timer {
	private:
		std::chrono::time_point<Clock> startTime;

	public:
		void start();
		double delta();
};

template <typename Clock>
void Timer<Clock>::start() {
	startTime = Clock::now();
}

template <typename Clock>
double Timer<Clock>::delta() {
	auto now = Clock::now();
	std::chrono::duration<double> dt = now - startTime;
	return dt.count();
}
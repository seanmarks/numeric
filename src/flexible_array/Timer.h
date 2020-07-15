// Simple timer class that uses std::chrono
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <exception>
#include <stdexcept>
#include <string>

class Timer
{
 public:
	using Clock     = std::chrono::high_resolution_clock;
	using TimePoint = std::chrono::time_point<Clock>;
	using Duration  = std::chrono::duration<double>;   // Period = std::ratio<1>  -->  seconds

	// Conversion factors
	static constexpr double s_to_ms = 1.0e3;  // s --> ms
	static constexpr double s_to_us = 1.0e6;  // s --> us (microseconds)

	Timer() {};

	Timer(const std::string& name):
		name_(name)
	{}

	// Starts the timer
	void start() {
		t_start_ = Clock::now();
		if ( is_on_ ) {
			throw std::runtime_error("Error: Timer \"" + name_ + "\" is already on.");
		}
		is_on_ = true;
	}

	// Stops the timer, and returns how long it was on
	// - Also stores timer duration for later retrieval
	double stop() {
		t_stop_ = Clock::now();
		if ( ! is_on_ ) {
			throw std::runtime_error("Error: Timer \"" + name_ + "\" was not on.");
		}
		is_on_ = false;
		duration_ = t_stop_ - t_start_;
		return duration_.count();
	}

	// Returns how long the timer has been on
	double get_time_elapsed_us() {
		TimePoint t_now = Clock::now();
		if ( ! is_on_ ) {
			throw std::runtime_error("Error: Timer \"" + name_ + "\" was not on.");
		}
		return std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_start_).count();
	}

	// Returns how long the timer was on
	double get_duration_() const {
		return duration_.count();
	}
	double get_duration_ms() const {
		return std::chrono::duration_cast<std::chrono::milliseconds>(duration_).count();
	}
	double get_duration_us() const {
		return std::chrono::duration_cast<std::chrono::microseconds>(duration_).count();
	}

 private:
	std::string name_ = "timer";

	bool is_on_ = false;
	TimePoint t_stop_, t_start_;             // [sec]
	Duration duration_ = Duration::zero();   // [sec]
};

#endif // ifndef TIMER_H

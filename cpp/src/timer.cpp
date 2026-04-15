#include "timer.h"

#include <chrono>

Timer::Timer() : start_time_(), end_time_(), is_running_(false) {}

void Timer::Reset() {
  is_running_ = false;
  start_time_ = Clock::time_point();
  end_time_ = Clock::time_point();
}

void Timer::Start() {
  start_time_ = Clock::now();
  is_running_ = true;
}

void Timer::Stop() {
  if (is_running_) {
    end_time_ = Clock::now();
    is_running_ = false;
  }
}

double Timer::GetTime() const {
  const auto end = is_running_ ? Clock::now() : end_time_;
  const std::chrono::duration<double> elapsed = end - start_time_;
  return elapsed.count();
}

#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
 public:
  Timer();

  void Reset();
  void Start();
  void Stop();
  double GetTime() const;

 private:
  using Clock = std::chrono::steady_clock;

  Clock::time_point start_time_;
  Clock::time_point end_time_;
  bool is_running_;
};

#endif

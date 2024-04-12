#pragma once


#include <chrono>
#include <iostream>

inline auto TimeStart() { return std::chrono::steady_clock::now(); }

inline void TimeStop(std::chrono::steady_clock::time_point start, const char* msg) {
    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << msg << " - completed in " << diff << " ms\n";
}
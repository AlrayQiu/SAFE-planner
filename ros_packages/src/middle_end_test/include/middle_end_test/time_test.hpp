#pragma once

#include <chrono>
#include <format>
#include <rclcpp/logger.hpp>
#include <string>
class TimeTest{
public:
TimeTest(int max):max_(max){}
    inline void Begin(){
        t_ = std::chrono::steady_clock::now();
    }
    inline void End(){
        auto epsilon = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_);
        seconds += epsilon.count() * 1e-9;
        window_++;
    }
    inline void Reset(){
        window_ = 0;   
        seconds = 0;
        t2_ = std::chrono::steady_clock::now();
    }
    inline bool Log(std::string& str_out){
        auto epsilon = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t2_);
        if(epsilon.count() < max_) return false;
        str_out = std::format("Avg:{}s , window size:{}", seconds / window_, window_);
        t2_ = std::chrono::steady_clock::now();
        return true;
    }
private:
    int window_ = 0;
    float seconds = 0;
    const float max_ = 1;
    std::chrono::steady_clock::time_point t_ = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point t2_ = std::chrono::steady_clock::now();
};
#pragma once
#include "trajectroies.hpp"
#include <algorithm>
#include <ctime>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace safe_planner::planner::trajectory
{
class MultiPoly final : public ITrajectory{
public:
MultiPoly() = default;
inline void set_param(const Eigen::MatrixX3d& b, const Eigen::VectorXd& t)
{
    builded = true;
    b_ = b;
    times_.clear();
    times_.reserve(t.size() + 1);
    times_.emplace_back(0);
    for(int i = 0, l = t.size();i < l;++i){
        times_.emplace_back(t(i) + times_[i].count());
    }
    int i = t.size() - 1;
    double t1 = t(i);
    double t2 = t1 * t1;
    double t3 = t2 * t1;
    double t4 = t3 * t1;
    double t5 = t4 * t1;
    Eigen::RowVector<double,6> beta(1,t1,t2,t3,t4,t5);
    tail_ = (beta * b_.block<6,3>(i * 6 ,0)).transpose();
}

inline Eigen::Vector3d pos(const seconds s) const {
    if(!builded) return Eigen::Vector3d::Zero();
    auto i = std::upper_bound(times_.begin(),times_.end(),s) - 1;
    if(i == times_.end() - 1) return tail_;
    if(i < times_.begin()) throw std::runtime_error("poly time err");
    auto id = std::distance(times_.begin(),i);
    auto t = (s - *i).count();
    auto t1 = t * t;
    auto t2 = t1 * t;
    auto t3 = t2 * t;
    auto t4 = t3 * t;
    Eigen::RowVector<double,6> beta(1, t, t1, t2, t3, t4);
    return (beta * b_.block<6,3>(id * 6, 0)).transpose();
};
inline double vel(const seconds s) const {
    if(!builded) return 0;
    auto i = std::upper_bound(times_.begin(),times_.end(),s) - 1;
    if(i == times_.end() - 1) return 0;
    if(i < times_.begin()) throw std::runtime_error("poly time err");
    auto id = std::distance(times_.begin(),i);
    auto t = (s - *i).count();
    auto t1 = t * t;
    auto t2 = t1 * t;
    auto t3 = t2 * t;
    Eigen::RowVector<double,6> beta(0, 1, 2 * t, 3 * t1, 4 * t2, 5 * t3);
    return (beta * b_.block<6,3>(id * 6, 0)).norm();
};
inline double acc(const seconds s) const {
    if(!builded) return 0;
    auto i = std::upper_bound(times_.begin(),times_.end(),s) - 1;
    if(i == times_.end() - 1) return 0;
    if(i < times_.begin()) throw std::runtime_error("poly time err");
    auto id = std::distance(times_.begin(),i);
    auto t = (s - *i).count();
    auto t1 = t * t;
    auto t2 = t1 * t;
    Eigen::RowVector<double,6> beta(0, 0, 2, 6 * t, 12 * t1, 20 * t2);
    return (beta * b_.block<6,3>(id * 6, 0)).norm();
};
inline double jerk(const seconds s) const {
    if(!builded) return 0;
    auto i = std::upper_bound(times_.begin(),times_.end(),s) - 1;
    if(i == times_.end() - 1) return 0;
    if(i < times_.begin()) throw std::runtime_error("poly time err");
    auto id = std::distance(times_.begin(),i);
    auto t = (s - *i).count();
    auto t1 = t * t;
    Eigen::RowVector<double,6> beta(0, 0, 0, 6, 24 * t, 60 * t1);
    return (beta * b_.block<6,3>(id * 6, 0)).norm();
};
inline double snap(const seconds s) const {
    if(!builded) return 0;
    auto i = std::upper_bound(times_.begin(),times_.end(),s) - 1;
    if(i == times_.end() - 1) return 0;
    if(i < times_.begin()) throw std::runtime_error("poly time err");
    auto id = std::distance(times_.begin(),i);
    auto t = (s - *i).count();
    Eigen::RowVector<double,6> beta(0, 0, 0, 0, 24, 120 * t);
    return (beta * b_.block<6,3>(id * 6, 0)).norm();
};
inline void get_time_range(seconds& from, seconds& to) const {
    if(!builded) return;
    from = *times_.begin();
    to = *times_.rbegin();
};

inline void control_point(std::vector<Eigen::Vector3d>& control){
    if(!builded) return;
    control.clear();
    for(int i = 0, l = times_.size() - 1;i < l;++i)
    {
        control.emplace_back(b_.row(i * 6).transpose());
    }
}

private:
bool builded = false;
Eigen::MatrixX3d b_;
Eigen::Vector3d tail_;
std::vector<seconds> times_;
};
}

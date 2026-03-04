#pragma once
#include "safe_planner/trajectory/trajectroies.hpp"
#include <Eigen/Eigen>
#include <algorithm>
#include <ctime>
#include <iterator>
#include <vector>

namespace safe_planner::planner::trajectory
{
class UniformBSpline final : public ITrajectory {
private:
Eigen::MatrixX3d Q_;
std::vector<seconds> time_;
double t1_;
double t2_;
double t3_;
std::tuple<int,double> get_s(const seconds t) const{
    auto c = std::upper_bound(time_.begin(),time_.end(),t) - 1;
    auto s = (t - *c).count() / ((c + 1)->count() - c->count());
    auto i = std::distance(time_.begin(),c);
    return {i,s};
}
public:
inline void set_param(const Eigen::MatrixX3d& Q, const float& t){
    Q_ = Q;
    time_.clear();
    auto st = seconds(t);
    time_.emplace_back(0);
    for(int i = 0, l = Q.rows() - 3;i < l;++i)
        time_.emplace_back(*time_.rbegin() + st);
    t1_ = t;
    t2_ = t1_ * t;
    t3_ = t2_ * t;
}

inline Eigen::Vector3d pos(const seconds t) const{
    if(t < *time_.begin() || t >= *time_.rbegin()) return Eigen::Vector3d::Zero();
    auto [i,s] = get_s(t);
    Eigen::RowVector4d S {1,s,s*s,s*s*s};
    return (S * M * Q_.block<4,3>(i,0)).transpose();
};


inline double vel(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    auto terr = t1_;
    Eigen::RowVector4d S {0,1.0 / terr,2.0 * s / terr,3.0 * s * s / terr};
    return (S * M * Q_.block<4,3>(i,0)).norm();
};


inline double acc(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    auto terr = t2_;
    Eigen::RowVector4d S {0, 0 ,2.0 / terr,6.0 * s / terr};
    return (S * M * Q_.block<4,3>(i,0)).norm();
};


inline double jerk(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    auto terr = t3_;
    Eigen::RowVector4d S {0, 0, 0, 6.0 / terr};
    return (S * M * Q_.block<4,3>(i,0)).norm();
};
inline double snap(const seconds) const{return 0;}

inline void get_time_range(seconds& from, seconds& to) const{
    from = *time_.begin();
    to = *time_.rbegin();
};

static inline const Eigen::Matrix4d M = (Eigen::Matrix4d() << 1./6., 4./6., 1./6., 0,
                                                             -3./6., 0    , 3./6 , 0,
                                                              3./6.,-6./6., 3./6., 0,
                                                             -1./6., 3./6., -3./6.,1./6.).finished();

static inline const Eigen::RowVector4d P0{1./6., 4./6., 1./6.,0};
static inline const Eigen::RowVector4d Pt{0, 1./6., 4./6., 1./6.};
static inline const Eigen::RowVector4d V0{-3./6., 0, 3./6., 0};
static inline const Eigen::RowVector4d Vt{0,-3./6., 0, 3./6.};
static inline const Eigen::RowVector4d A0{1, -2, 1, 0};
static inline const Eigen::RowVector4d At{0, 1, -2, 1};
static inline const Eigen::RowVector4d J{-1, 3, -3, 1};
};
}
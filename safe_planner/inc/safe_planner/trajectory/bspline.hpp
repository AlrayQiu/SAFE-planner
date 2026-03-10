#pragma once
#include "safe_planner/trajectory/trajectroies.hpp"
#include "safe_planner/trajectory/utils.hpp"
#include <Eigen/Eigen>
#include <algorithm>
#include <ctime>
#include <iterator>
#include <vector>

namespace safe_planner::planner::trajectory
{
template<int k>
class UniformBSpline_K final : public ITrajectory {
private:
Eigen::MatrixX3d Q_;
std::vector<seconds> time_;
double t1_;
double t2_;
double t3_;
double t4_;
std::tuple<int,double> get_s(const seconds t) const{
    auto c = std::upper_bound(time_.begin(),time_.end(),t) - 1;
    auto s = (t - *c).count() / ((c + 1)->count() - c->count());
    auto i = std::distance(time_.begin(),c);
    return {i,s};
}
public:

inline void control_point(std::vector<Eigen::Vector3d>& control){
    if(!build) return;
    control.clear();
    for(int i = 0;i < Q_.rows();++i){
        control.emplace_back(Q_.row(i).transpose());
    }
}
inline void set_param(const Eigen::MatrixX3d& Q, const float& t){
    build = true;
    Q_ = Q;
    time_.clear();
    auto st = seconds(t);
    time_.emplace_back(0);
    for(int i = 0, l = Q.rows() - k + 1;i < l;++i)
        time_.emplace_back(*time_.rbegin() + st);
    t1_ = t;
    t2_ = t1_ * t;
    t3_ = t2_ * t;
    t4_ = t2_ * t2_;
}

inline Eigen::Vector3d pos(const seconds t) const{
    if(t < *time_.begin() || t >= *time_.rbegin()) return Eigen::Vector3d::Zero();
    auto [i,s] = get_s(t);
    Eigen::RowVector<double,k> S = utils::make_S<k>(s, 0);
    return (S * M * Q_.block<k,3>(i,0)).transpose();
};


inline double vel(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    Eigen::RowVector<double,k> S = utils::make_S<k>(s, 1) / t1_;
    return (S * M * Q_.block<k,3>(i,0)).norm();
};


inline double acc(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    Eigen::RowVector<double,k> S = utils::make_S<k>(s, 2) / t2_;
    return (S * M * Q_.block<k,3>(i,0)).norm();
};


inline double jerk(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    Eigen::RowVector<double,k> S = utils::make_S<k>(s, 3) / t3_;
    return (S * M * Q_.block<k,3>(i,0)).norm();
};
inline double snap(const seconds t) const{
    if(t < *time_.begin() || t >= * time_.rbegin()) return 0;
    auto [i,s] = get_s(t);
    Eigen::RowVector<double,k> S = utils::make_S<k>(s, 3) / t4_;
    return (S * M * Q_.block<k,3>(i,0)).norm();
}

inline void get_time_range(seconds& from, seconds& to) const{
    from = *time_.begin();
    to = *time_.rbegin();
};
bool build = false;

static inline const Eigen::Matrix<double,k,k> M = utils::make_M<k>();

static inline const Eigen::RowVector<double,k> P0 = (utils::make_S<k>(0,0) * M);
static inline const Eigen::RowVector<double,k> Pt = (utils::make_S<k>(1,0) * M);
static inline const Eigen::RowVector<double,k> V0 = (utils::make_S<k>(0,1) * M);
static inline const Eigen::RowVector<double,k> Vt = (utils::make_S<k>(1,1) * M);
static inline const Eigen::RowVector<double,k> A0 = (utils::make_S<k>(0,2) * M);
static inline const Eigen::RowVector<double,k> At = (utils::make_S<k>(1,2) * M);
static inline const Eigen::RowVector<double,k> J0 = (utils::make_S<k>(0,3) * M);
static inline const Eigen::RowVector<double,k> Jt = (utils::make_S<k>(1,3) * M);
}; 
typedef UniformBSpline_K<4> UniformBSpline; 
}
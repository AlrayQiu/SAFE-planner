#pragma once

#include "i_optimizable.hpp"
#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/trajectory/bspline.hpp"
#include "uniform_bspline_opt.hpp"
#include <cmath>
#include <vector>
namespace safe_planner::planner::middle_end::optimizer {
template <class TMap>
class ParameterizedUrbsOpt final : public IOptimizeable{
public:
inline ParameterizedUrbsOpt(
    trajectory::impl::UniformBSplineImpl& urbs,
    const esdf::ImplicitESDF<TMap>& esdf,
    const TMap& map)
: urbs_(urbs)
, urbs_opt_(urbs,esdf,map)
, t_index(urbs.inner_points_count_ * 3)
{
    dj_dp_.resize(urbs.inner_points_count_,3);
}
inline void get_j(double& j) const{
    j = j_;
};
inline void get_g(Eigen::VectorXd& g) const{
    for(int i = 0;i < urbs_.inner_points_count_;++i){
        g.segment<3>(i * 3) = (dj_dp_.row(i)).transpose();
    }
    g(t_index) = dj_dtau_;
};
inline void set_x(const Eigen::VectorXd& x){
    points_.clear();
    for(int i = 0;i < urbs_.inner_points_count_;++i){
        points_.emplace_back(x.segment<3>(i*3));
    }
    urbs_.set_param(points_, tau_to_t(x(t_index)));
    reset(x);
};
inline void init_x(Eigen::VectorXd& x) const{
    x.resize(t_index + 1);
    for(int i = 0;i < urbs_.inner_points_count_;++i){
        x.segment<3>(i * 3) = (
            trajectory::UniformBSpline::P0 *
            urbs_.Q_.block<4,3>(2 + i,0)).transpose();
    }
    x(t_index) = t_to_tau(urbs_.t_);
};

inline void reset(const Eigen::VectorXd& x){
    j_ = 0;
    dj_dtau_ = 0;
    dj_dp_.setZero();
    urbs_opt_.reset();
    urbs_opt_.get_j(j_);
    urbs_opt_.dj_dq_ = urbs_.lu_.transpose().solve(urbs_opt_.dj_dq_);
    dj_dp_ = urbs_opt_.dj_dq_.block(3,0,urbs_.inner_points_count_,3);
    dj_dtau_ = urbs_opt_.dj_dt_ - prop_c_to_dt(urbs_opt_.dj_dq_);
    
    dt_to_dtau(dj_dtau_, x(t_index));
    // only_for_test(x);
    auto x1 = x.segment(0,t_index);
    tikhonov_normalize(x1);
    // std::cerr << urbs_.t_ << ' ' << x(t_index) <<' ' << tau_to_t(x(t_index)) << std::endl;
}

private:
inline void only_for_test(const Eigen::VectorXd& x){
    Eigen::Vector3f nearp;
    for(int i = 0;i < urbs_.inner_points_count_;++i){
        auto vp = (x.segment<3>(i*3) - urbs_.points[i]).squaredNorm() - 0.25;
        if(vp < 0) continue;
        j_ += vp * vp * vp;
        dj_dp_.row(i) += 6 * vp * vp * (x.segment<3>(i*3) - urbs_.points[i]).transpose();
    }
}

inline void tikhonov_normalize(const Eigen::VectorXd& x){
    j_ += tikhonov_epsilon / 2 * x.squaredNorm();
    for(int i = 0;i < urbs_.inner_points_count_;++i){
        dj_dp_.row(i) += tikhonov_epsilon * x.segment<3>(i * 3).transpose();
    }
}

inline double prop_c_to_dt(const Eigen::MatrixX3d& G){
    const auto Q1 = urbs_.Q_.topRows(4);
    const auto Q2 = urbs_.Q_.bottomRows(4);

    const auto t1 = urbs_.t_;
    const auto t2 = t1 * t1; 
    const auto t3 = t2 * t1; 
    const Eigen::Matrix<double,4,3> MQ = (trajectory::UniformBSpline::M * Q2);
    const Eigen::RowVector3d dv0 = -trajectory::UniformBSpline::V0 / t2 * Q1;
    const Eigen::RowVector3d da0 = -2 * trajectory::UniformBSpline::A0 / t3 * Q1;
    const Eigen::RowVector3d dpt = Eigen::RowVector4d{0, -1, -2,-3} / t1 * MQ;
    const Eigen::RowVector3d dvt = Eigen::RowVector4d{0, -1, -4 * 1,-9 * 1} / t2 * MQ;
    const Eigen::RowVector3d dat = Eigen::RowVector4d{0, 0, -4, -18} / t3 * MQ;

    return G.row(1).dot(dv0) + G.row(2).dot(da0) 
         + G.row(urbs_.control_points_count_ - 3).dot(dpt)
         + G.row(urbs_.control_points_count_ - 2).dot(dvt)
         + G.row(urbs_.control_points_count_ - 1).dot(dat); 
}

inline void dt_to_dtau(auto& dj_dt ,const auto& tau) {
    double dtau_dt;
    if(tau > 0){
        dtau_dt = tau + 1.0;
    }else{
        auto den_sqrt = (0.5 * tau - 1.0) * tau + 1.0;
        dtau_dt = (1 - tau) / (den_sqrt * den_sqrt);
    }
    dj_dt = (dj_dt + time_weight * urbs_.segment_count_) * dtau_dt;
    j_ += time_weight * urbs_.segment_count_ * urbs_.t_;
}

inline double tau_to_t(const double tau) const{
    if(tau > 0){
        return (0.5 * tau + 1) * tau + 1;
    }else{
        return 1 / ((0.5 * tau - 1) * tau + 1);
    }
}
inline double t_to_tau(const double t) const{
    if(t > 1){
        return std::sqrt(2 * t - 1) - 1;
    }else{
        return 1 - std::sqrt(2 / t - 1);
    }
}


double j_;
double dj_dtau_;
Eigen::MatrixX3d dj_dp_;
const double tikhonov_epsilon = 1e-6; 
const double time_weight = 1;
trajectory::impl::UniformBSplineImpl& urbs_;
UniformBSplineOpt<TMap> urbs_opt_;
const int t_index;

std::vector<Eigen::Vector3d> points_{};
};
}
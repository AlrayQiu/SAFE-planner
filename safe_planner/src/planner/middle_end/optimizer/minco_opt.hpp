#pragma once

#include "i_optimizable.hpp"
#include "../trajectory/minco_impl.hpp"
#include "polynomial_opt.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

namespace safe_planner::planner::middle_end::optimizer {
template<class TMap>
class MincoOpt final : public IOptimizeable{
public:
inline MincoOpt(
    trajectory::impl::MincoImpl& minco,
    const esdf::ImplicitESDF<TMap>& esdf,
    const TMap& map)
: minco_(minco)
, poly_(minco_.poly_, esdf, map){
    inner_points_.resize(3, minco_.inner_points_count_);
    t_.resize(minco_.pieces_num_);
    dj_dp_.resize(3 * minco_.inner_points_count_);
    dj_dtau_.resize(minco_.pieces_num_);
}

inline void init_x(Eigen::VectorXd& x) const
{
    x.resize(minco_.inner_points_count_ * 3 + minco_.pieces_num_);
    for(int i = 0;i < minco_.inner_points_count_; ++i){
        x.segment<3>(i * 3) = minco_.poly_.b_.row(6 * (i + 1)).transpose();
    }
    for(int i = 0;i < minco_.pieces_num_;++i){
        x(i + minco_.inner_points_count_ * 3) = t_to_tau(minco_.poly_.t1_(i));
    }
}

inline void get_j(double& j) const{ j = j_;};
inline void get_g(Eigen::VectorXd& g) const {
    g.segment(0, minco_.inner_points_count_ * 3) = dj_dp_;
    g.segment(minco_.inner_points_count_ * 3, minco_.pieces_num_) = dj_dtau_;
};
inline void set_x(const Eigen::VectorXd& x){
    for(int i = 0;i < minco_.inner_points_count_;++i){
        inner_points_.col(i) = x.segment<3>(i * 3);
    }
    tau_ = x.segment(3 * minco_.inner_points_count_, minco_.pieces_num_);
    for(int i = 0;i < minco_.pieces_num_;++i)
    {
        t_(i) = tau_to_t(tau_(i));
    }
    minco_.set_param(inner_points_, t_);
    
    reset();
};

inline void reset(){
    poly_.reset();
    minco_.m_.solve_adj(poly_.dj_dc_);
    j_ = poly_.j_;
    add_porp_c_to_p();
    add_prop_c_to_dt();
    dt_to_dtau(poly_.dj_dt_, tau_);
    dj_dp_.setZero();
    j_ += t_.sum() * weight_time_;
}

private:
inline void add_porp_c_to_p(){
    for(int i = 0;i < minco_.inner_points_count_;++i)
        dj_dp_.segment<3>(i * 3) = poly_.dj_dc_.row(i * 6 + 5).transpose();
}
inline void add_prop_c_to_dt(){
    Eigen::MatrixXd B1(6, 3), B2(3, 3);
    Eigen::RowVector3d negVel, negAcc, negJer, negSnp, negCrk;
    const int N = minco_.pieces_num_;
    for (int i = 0; i < N - 1; i++)
    {
        // calculate \frac{\partial E_i}{\partial T_i}. E_i is 6*6 matrix
        // According the matrix layout, E_i = [\beta`3(T) \beta`4(T) \beta(T) \beta(T) \beta`1(T) \beta`2]^T
        // so \frac{\partial E_i}{\partial T_i} = [\beta`4(T) \beta`5(T) \beta`1(T) \beta`1(T) \beta`2(T) \beta`3]^T
        // multiply by c_i and get the [Snap, Crk, Vel, Vel, Acc, Jerk]
        negVel = -(minco_.poly_.b_.row(i * 6 + 1) +
                   2.0 * minco_.poly_.t1_(i) * minco_.poly_.b_.row(i * 6 + 2) +
                   3.0 * minco_.poly_.t2_(i) * minco_.poly_.b_.row(i * 6 + 3) +
                   4.0 * minco_.poly_.t3_(i) * minco_.poly_.b_.row(i * 6 + 4) +
                   5.0 * minco_.poly_.t4_(i) * minco_.poly_.b_.row(i * 6 + 5));
        negAcc = -(2.0 * minco_.poly_.b_.row(i * 6 + 2) +
                   6.0 * minco_.poly_.t1_(i) * minco_.poly_.b_.row(i * 6 + 3) +
                   12.0 * minco_.poly_.t2_(i) * minco_.poly_.b_.row(i * 6 + 4) +
                   20.0 * minco_.poly_.t3_(i) * minco_.poly_.b_.row(i * 6 + 5));
        negJer = -(6.0 * minco_.poly_.b_.row(i * 6 + 3) +
                   24.0 * minco_.poly_.t1_(i) * minco_.poly_.b_.row(i * 6 + 4) +
                   60.0 * minco_.poly_.t2_(i) * minco_.poly_.b_.row(i * 6 + 5));
        negSnp = -(24.0 * minco_.poly_.b_.row(i * 6 + 4) +
                   120.0 * minco_.poly_.t1_(i) * minco_.poly_.b_.row(i * 6 + 5));
        negCrk = -120.0 * minco_.poly_.b_.row(i * 6 + 5);

        B1 << negSnp, negCrk, negVel, negVel, negAcc, negJer;
        // Eq.(10). 
        poly_.dj_dt_(i) += B1.cwiseProduct(poly_.dj_dc_.template block<6, 3>(6 * i + 3, 0)).sum();
    }
    // The E_M should be specially handled because it's 3*6
    negVel = -(minco_.poly_.b_.row(6 * N - 5) +
               2.0 * minco_.poly_.t1_(N - 1) * minco_.poly_.b_.row(6 * N - 4) +
               3.0 * minco_.poly_.t2_(N - 1) * minco_.poly_.b_.row(6 * N - 3) +
               4.0 * minco_.poly_.t3_(N - 1) * minco_.poly_.b_.row(6 * N - 2) +
               5.0 * minco_.poly_.t4_(N - 1) * minco_.poly_.b_.row(6 * N - 1));
    negAcc = -(2.0 * minco_.poly_.b_.row(6 * N - 4) +
               6.0 * minco_.poly_.t1_(N - 1) * minco_.poly_.b_.row(6 * N - 3) +
               12.0 * minco_.poly_.t2_(N - 1) * minco_.poly_.b_.row(6 * N - 2) +
               20.0 * minco_.poly_.t3_(N - 1) * minco_.poly_.b_.row(6 * N - 1));
    negJer = -(6.0 * minco_.poly_.b_.row(6 * N - 3) +
               24.0 * minco_.poly_.t1_(N - 1) * minco_.poly_.b_.row(6 * N - 2) +
               60.0 * minco_.poly_.t2_(N - 1) * minco_.poly_.b_.row(6 * N - 1));
    B2 << negVel, negAcc, negJer;
    poly_.dj_dt_(N - 1) += B2.cwiseProduct(poly_.dj_dc_.template block<3, 3>(6 * N - 3, 0)).sum();
    return;
}

inline void dt_to_dtau(const auto& dj_dt ,const auto& tau){
    double dtau_dt;
    for(int i = 0, l = dj_dtau_.size();i < l;++i){
        if(tau(i) > 0){
            dtau_dt = tau(i) + 1.0;
        }else{
            auto den_sqrt = (0.5 * tau(i) - 1.0) * tau(i) + 1.0;
            dtau_dt = (1 - tau(i)) / (den_sqrt * den_sqrt);
        }
        dj_dtau_(i) = (dj_dt(i) + weight_time_) * dtau_dt;

    }
}

inline double tau_to_t(const double tau){
    if(tau > 0){
        return (0.5 * tau + 1) * tau + 1;
    }else{
        return 1 / ((0.5 * tau - 1) * tau + 1);
    }
}
inline double t_to_tau(const double t){
    if(t > 1){
        return std::sqrt(2 * t - 1) - 1;
    }else{
        return 1 - std::sqrt(2 / t - 1);
    }
}

trajectory::impl::MincoImpl& minco_;
optimizer::PolynomialOpt<TMap> poly_;
Eigen::Matrix3Xd inner_points_;
Eigen::VectorXd t_;
Eigen::VectorXd tau_;
Eigen::VectorXd dj_dp_;
Eigen::VectorXd dj_dtau_;

double j_;
const double weight_time_ = 1;
};
}
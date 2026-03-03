#pragma once

#include <Eigen/Eigen>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "../../utils/banded_system.hpp"
#include "polynomial_impl.hpp"
#include "safe_planner/trajectory/multi_poly.hpp"

namespace safe_planner::planner::trajectory::impl {
class MincoImpl final{
public:

MincoImpl() = delete;

inline void set_param(const Eigen::Matrix3Xd& inner_points, const Eigen::VectorXd& t){
    poly_.set_t(t);
    m_.reset();
    poly_.b_.setZero();
    m_(0,0) = 1.0;
    m_(1,1) = 1.0;
    m_(2,2) = 2.0;
    poly_.b_.row(0) = head_.col(0).transpose();
    poly_.b_.row(1) = head_.col(1).transpose();
    poly_.b_.row(2) = head_.col(2).transpose();


    if(inner_points_count_ != inner_points.cols()) throw std::runtime_error("minco set_param error with innerpoint");
    if(pieces_num_ != t.size()) throw std::runtime_error("minco set_param error with times");

    for(int i = 0;i < inner_points_count_;++i){
        const int curr = 6 * i;

        m_(curr + 3, curr + 3) = 6.0;
        m_(curr + 3, curr + 4) = 24.0 * poly_.t1_(i);
        m_(curr + 3, curr + 5) = 60.0 * poly_.t2_(i);
        m_(curr + 3, curr + 9) = -6.0;

        m_(curr + 4, curr + 4) = 24.0;
        m_(curr + 4, curr + 5) = 120.0 * poly_.t1_(i);
        m_(curr + 4, curr + 10) = -24.0;

        m_(curr + 5, curr) = 1.0;
        m_(curr + 5, curr + 1) = poly_.t1_(i);
        m_(curr + 5, curr + 2) = poly_.t2_(i);
        m_(curr + 5, curr + 3) = poly_.t3_(i);
        m_(curr + 5, curr + 4) = poly_.t4_(i);
        m_(curr + 5, curr + 5) = poly_.t5_(i);

        m_(curr + 6, curr) = 1.0;
        m_(curr + 6, curr + 1) = poly_.t1_(i);
        m_(curr + 6, curr + 2) = poly_.t2_(i);
        m_(curr + 6, curr + 3) = poly_.t3_(i);
        m_(curr + 6, curr + 4) = poly_.t4_(i);
        m_(curr + 6, curr + 5) = poly_.t5_(i);
        m_(curr + 6, curr + 6) = -1.0;
        
        m_(curr + 7, curr + 1) = 1.0;
        m_(curr + 7, curr + 2) = 2.0 * poly_.t1_(i);
        m_(curr + 7, curr + 3) = 3.0 * poly_.t2_(i);
        m_(curr + 7, curr + 4) = 4.0 * poly_.t3_(i);
        m_(curr + 7, curr + 5) = 5.0 * poly_.t4_(i);
        m_(curr + 7, curr + 7) = -1.0;

        m_(curr + 8, curr + 2) = 2.0;
        m_(curr + 8, curr + 3) = 6.0 * poly_.t1_(i);
        m_(curr + 8, curr + 4) = 12.0 *poly_. t2_(i);
        m_(curr + 8, curr + 5) = 20.0 *poly_. t3_(i);
        m_(curr + 8, curr + 8) = -2.0;

        poly_.b_.row(curr + 5) = inner_points.col(i).transpose();
    }

    const int curr = 6 * pieces_num_;
    m_(curr - 3, curr - 6) = 1.0;
    m_(curr - 3, curr - 5) = poly_.t1_(inner_points_count_);
    m_(curr - 3, curr - 4) = poly_.t2_(inner_points_count_);
    m_(curr - 3, curr - 3) = poly_.t3_(inner_points_count_);
    m_(curr - 3, curr - 2) = poly_.t4_(inner_points_count_);
    m_(curr - 3, curr - 1) = poly_.t5_(inner_points_count_);

    m_(curr - 2, curr - 5) = 1.0;
    m_(curr - 2, curr - 4) = 2.0 * poly_.t1_(inner_points_count_);
    m_(curr - 2, curr - 3) = 3.0 * poly_.t2_(inner_points_count_);
    m_(curr - 2, curr - 2) = 4.0 * poly_.t3_(inner_points_count_);
    m_(curr - 2, curr - 1) = 5.0 * poly_.t4_(inner_points_count_);

    m_(curr - 1, curr - 4) = 2.0;
    m_(curr - 1, curr - 3) = 6.0 *  poly_.t1_(inner_points_count_);
    m_(curr - 1, curr - 2) = 12.0 * poly_.t2_(inner_points_count_);
    m_(curr - 1, curr - 1) = 20.0 * poly_.t3_(inner_points_count_);

    poly_.b_.row(curr - 3) = tail_.col(0).transpose();
    poly_.b_.row(curr - 2) = tail_.col(1).transpose();
    poly_.b_.row(curr - 1) = tail_.col(2).transpose();
    m_.factorize_lu();
    m_.solve(poly_.b_);
}


inline MincoImpl(
    const std::vector<Eigen::Vector3d>& control_points, 
    const Eigen::Matrix<double,3,2>& head_va,
    const Eigen::Matrix<double,3,2>& tail_va):
    poly_(),
    pieces_num_(poly_.pieces_num_)
{   
    static std::vector<double> times{};
    static std::vector<Eigen::Vector3d> points{};
    static Eigen::VectorXd t;
    static Eigen::Matrix3Xd p;

    times.clear();
    points.clear();
    points.emplace_back(control_points[0]);
    
    for(int i = 1, l = control_points.size() - 1;i < l;++i){
        const auto perr = control_points[i] - *points.rbegin();
        const auto dis = (perr).norm();
        if(dis < init_min_distance_)
            continue;
        const auto j = static_cast<int>((dis / init_step_distance_) + 1);
        const auto step = 1.0 / j;
        
        for(int k = 1;k <= j;++k){
            points.emplace_back(*points.rbegin() + perr * step);
            times.emplace_back(dis * step);
        }
    }
    
    const auto perr = *control_points.rbegin() - *points.rbegin();
    const auto dis = (perr).norm();
    const auto j = static_cast<int>((dis / init_step_distance_) + 1);
    const auto step = 1.0 / j;
    for(int k = 1;k <= j;++k){
        points.emplace_back(*points.rbegin() + perr * step);
        times.emplace_back(dis * step);
    }
    *times.begin() *= 3;
    *times.rbegin() *= 3;
    if(times.size() == 1){
        points.emplace_back(*control_points.rbegin());
        times.emplace_back(0.1);
    }
    
    pieces_num_ = times.size();
    inner_points_count_ = pieces_num_ - 1;
    
    t.setZero(pieces_num_);
    p.setZero(3,inner_points_count_);
    for(int i = 1, l = inner_points_count_;i <= l; ++i){
        p.col(i - 1) = points[i];
        t(i - 1) = times[i - 1];
    }
    t(inner_points_count_) = times[inner_points_count_];
    head_.col(0) = *control_points.begin();
    head_.block<3,2>(0,1) = head_va;
    tail_.col(0) = *control_points.rbegin();
    tail_.block<3,2>(0,1) = tail_va;
    m_.create(6 * pieces_num_, 6, 6);
    poly_.b_.resize(6 * pieces_num_, 3);
    set_param(p, t);
}

inline void build(MultiPoly& poly){
    poly.set_param(poly_.b_, poly_.t1_);
}

const double init_min_distance_ = 0.01;
const double init_step_distance_ = 1;

utils::BandedSystem m_;
Eigen::Matrix3d head_;
Eigen::Matrix3d tail_;
PolynomialImpl poly_;

int inner_points_count_;
int& pieces_num_;
};
}
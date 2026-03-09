#pragma once
#include "safe_planner/trajectory/bspline.hpp"
#include <Eigen/Eigen>
#include <cmath>
#include <vector>
namespace safe_planner::planner::trajectory::impl{
class UniformBSplineImpl{
public:
inline UniformBSplineImpl(
    const std::vector<Eigen::Vector3d>& control_points, 
    const Eigen::Matrix<double,3,2>& head_va,
    const Eigen::Matrix<double,3,2>& tail_va)
    {
    head_.col(0) = *control_points.begin();
    head_.block<3,2>(0,1) = head_va;
    tail_.col(0) = *control_points.rbegin();
    tail_.block<3,2>(0,1) = tail_va;

    points.clear();
    points.emplace_back(control_points[0]);
    
    for(int i = 1, l = control_points.size();i < l;++i){
        const Eigen::Vector3d perr = control_points[i] - points.back();
        const auto dis = (perr).norm();
        if(dis < init_min_distance_)
            continue;
        const auto j = static_cast<int>((dis / init_step_distance_) + 1);
        const auto step = 1.0 / j;
        
        for(int k = 1;k <= j;++k){
            points.emplace_back(points.back() + perr * step);
        }
    }
    
    points.erase(points.begin());
    points.pop_back();
    
    inner_points_count_ = std::max<int>(points.size(),0);
    segment_count_ = inner_points_count_ + 3;
    control_points_count_ = inner_points_count_ + 6;

    
    M_.setZero(control_points_count_,control_points_count_);
    M_.block<1,4>(0,0) = UniformBSpline::P0;
    for(int i = 0;i < inner_points_count_;i++){
        M_.block<1,4>(3 + i, i + 2) = UniformBSpline::P0;
    }
    M_.block<1,4>(3 + inner_points_count_, inner_points_count_ + 2) = UniformBSpline::Pt;

    set_param(points, 1);
}
int segment_count_;
int inner_points_count_ = 0;
int control_points_count_;

double t_ = 1;

Eigen::Matrix3d head_;
Eigen::Matrix3d tail_;


Eigen::MatrixX3d Q_;
Eigen::MatrixXd M_;
Eigen::PartialPivLU<Eigen::MatrixXd> lu_;

inline void build(UniformBSpline& spline){
    spline.set_param(Q_, t_);
}

inline void set_param(const std::vector<Eigen::Vector3d>& inner_points,const double t){
    Q_.setZero(control_points_count_,3);
    Q_.block<3,3>(0,0) = head_.transpose();

    for(int i = 0;i < inner_points_count_;i++){
        Q_.row(3 + i) = inner_points[i].transpose();
    }
    Q_.block<3,3>(3 + inner_points_count_, 0) = tail_.transpose();
    t_ = t;
    const auto t1 = t;
    const auto t2 = t1 * t1;
    M_.block<1,4>(1,0) = UniformBSpline::V0 / t1;
    M_.block<1,4>(2,0) = UniformBSpline::A0 / t2;
    M_.block<1,4>(4 + inner_points_count_, inner_points_count_ + 2) = UniformBSpline::Vt / t1;
    M_.block<1,4>(5 + inner_points_count_, inner_points_count_ + 2) = UniformBSpline::At / t2;
    lu_.compute(M_);
    Q_ = lu_.solve(Q_);
    // std::cerr << M_.fullPivLu().rank() << " " << M_.rows() << std::endl;
} 

std::vector<Eigen::Vector3d> points; 

private:
const double init_min_distance_ = 0.01;
const double init_step_distance_ = 0.5;

};
}
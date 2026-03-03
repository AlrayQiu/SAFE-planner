#pragma once
#include "safe_planner/trajectory/bspline.hpp"
#include <Eigen/Eigen>
#include <vector>
namespace safe_planner::planner::trajectory::impl{
class UniformBSplineImpl{
public:
inline UniformBSplineImpl(
    const std::vector<Eigen::Vector3d>& control_points, 
    const Eigen::Matrix<double,3,2>& head_va,
    const Eigen::Matrix<double,3,2>& tail_va){
    head_.col(0) = *control_points.begin();
    head_.block<3,2>(0,1) = head_va;
    tail_.col(0) = *control_points.rbegin();
    tail_.block<3,2>(0,1) = tail_va;

    static std::vector<Eigen::Vector3d> points; 
    points.clear();
    points.emplace_back(control_points[0]);
    
    for(int i = 1, l = control_points.size();i < l;++i){
        const auto perr = control_points[i] - *points.rbegin();
        const auto dis = (perr).norm();
        if(dis < init_min_distance_)
            continue;
        const auto j = static_cast<int>((dis / init_step_distance_) + 1);
        const auto step = 1.0 / j;
        
        for(int k = 1;k <= j;++k){
            points.emplace_back(*points.rbegin() + perr * step);
        }
    }
    
    points.erase(points.begin());
    points.pop_back();
    inner_points_count_ = points.size();
    set_param(points);
}

const double init_min_distance_ = 0.01;
const double init_step_distance_ = 1;

Eigen::Matrix3d head_;
Eigen::Matrix3d tail_;

Eigen::MatrixX3d Q_;

void build(UniformBSpline& spline){
    spline.set_param(Q_, 1);
}

private:
void set_param(const std::vector<Eigen::Vector3d>& inner_points){
    Q_.setZero(inner_points_count_ + 6,3);
    M_.setZero(inner_points_count_ + 6,inner_points_count_ + 6);
    M_.block<1,4>(0,0) = UniformBSpline::P0;
    M_.block<1,4>(1,0) = UniformBSpline::V0;
    M_.block<1,4>(2,0) = UniformBSpline::A0;
    Q_.block<3,3>(0,0) = head_.transpose();

    for(int i = 0;i < inner_points_count_;i++){
        M_.block<1,4>(3 + i, i + 2) = UniformBSpline::P0;
        Q_.row(3 + i) = inner_points[i].transpose();
    }
    M_.block<1,4>(3 + inner_points_count_, inner_points_count_ + 2) = UniformBSpline::Pt;
    M_.block<1,4>(4 + inner_points_count_, inner_points_count_ + 2) = UniformBSpline::Vt;
    M_.block<1,4>(5 + inner_points_count_, inner_points_count_ + 2) = UniformBSpline::At;
    Q_.block<3,3>(3 + inner_points_count_, 0) = tail_.transpose();
    Q_ = M_.fullPivLu().solve(Q_);
    // std::cerr << M_.fullPivLu().rank() << " " << M_.rows() << std::endl;
} 

int inner_points_count_ = 0;

Eigen::MatrixXd M_;
};
}
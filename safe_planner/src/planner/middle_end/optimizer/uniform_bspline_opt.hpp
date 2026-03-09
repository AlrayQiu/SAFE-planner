#pragma once

#include "ego_utils.hpp"
#include "i_optimizable.hpp"
#include "../trajectory/uniform_bspline_impl.hpp"
#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/trajectory/bspline.hpp"

#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

namespace safe_planner::planner::middle_end::optimizer {
template<class TMap>
class UniformBSplineOpt final : public IOptimizeable{
public:
inline UniformBSplineOpt(
    trajectory::impl::UniformBSplineImpl& urbs,
    const esdf::ImplicitESDF<TMap>& esdf,
    const TMap& map
)
: urbs_(urbs)
, esdf_(esdf)
, map_(map)
, ego_(urbs.inner_points_count_ + 4,map)
{
    dj_dq_.resize(urbs.Q_.rows(),urbs.Q_.cols());
    cps.resize(urbs.Q_.cols(), urbs.inner_points_count_ + 4);
    cps.col(0) = urbs.head_.col(0);
    cps.col(urbs.inner_points_count_ + 3) = urbs.tail_.col(0);
    gradp.resize(cps.rows(), cps.cols());
}

inline void get_j(double& j) const{
    j = j_;
}
inline void get_g(Eigen::VectorXd& g) const{
    for(int i = 0;i < urbs_.control_points_count_;++i){
        g.segment<3>(i * 3) = dj_dq_.row(i).transpose();    
    }
}
inline void set_x(const Eigen::VectorXd& x){
    for(int i = 0;i < urbs_.control_points_count_;++i){
        urbs_.Q_.row(i) = x.segment<3>(i * 3).transpose();    
    }
    reset();
}
inline void init_x(Eigen::VectorXd& x) const{
    x.resize(dj_dq_.size());
    for(int i = 0;i < urbs_.control_points_count_;++i){
        x.segment<3>(i * 3) = urbs_.Q_.row(i).transpose();    
    }
}

inline void reset(){
    j_ = 0;
    dj_dq_.setZero();
    dj_dt_ = 0;
    ++iter_num;
    minimal_jerk();
    collision_cost();
    feasibility_cost();
}

double j_;
Eigen::MatrixX3d dj_dq_;
double dj_dt_;
private:
inline void minimal_jerk()
{
    const double t = urbs_.t_;
    const double t3 = t * t * t;
    const double t6 = t3 * t3;
    const Eigen::RowVector4d  J = trajectory::UniformBSpline::J / t3;
    const Eigen::Matrix4d jTj = JTJ / t6;

    for(int i = 0;i < urbs_.segment_count_;++i){
        auto Q = urbs_.Q_.block<4,3>(i,0);
        const auto j =  omega_smooth_ * (J * Q).squaredNorm();
        j_ += j * t;

        dj_dq_.block<4,3>(i,0) += (2 * omega_smooth_ * t * jTj * Q);
        dj_dt_ += -5 * j; 
    } 

    
}
inline void feasibility_cost(){
    const auto t1 = urbs_.t_;
    const auto t2 = t1 * t1;
    const auto t3 = t2 * t1;
    const auto step = 1.0 / K;
    Eigen::RowVector3d dv;
    Eigen::RowVector3d da;
    double cv,ca;
    for(int i = 0, l = urbs_.segment_count_ - 1;i <= l;++i){
        auto s1 = 0.0;
        const auto Q = urbs_.Q_.block<4,3>(i,0);

        for(int j = 0;j <= K; ++j)
        {
            const auto s2 = s1 * s1;

            const auto omega = 
                ((j == 0 || j == K) && (i != 0 && i != l)) 
                ? 0.5 : 1.0;
            const Eigen::Matrix<double,4,3> MQ = (trajectory::UniformBSpline::M * Q);
            const Eigen::RowVector4d B2 = 
                Eigen::RowVector4d{0, 1, 2 * s1,3 * s2} / t1 * trajectory::UniformBSpline::M;
            const Eigen::RowVector4d B3 = 
                Eigen::RowVector4d{0, 0, 2, 6 * s1} / t2 * trajectory::UniformBSpline::M;
            const Eigen::RowVector3d v = (B2 * Q);
            const Eigen::RowVector3d a = (B3 * Q);

            if(feasibility_v(v, dv, cv)){
                j_ += step * omega * cv;
                dj_dq_.block<4,3>(i,0) += step * omega * B2.transpose() * dv;
                const Eigen::RowVector3d dv_dt = (Eigen::RowVector4d{0, -1, -4 * s1,-9 * s2} / t2 * MQ);
                dj_dt_ += step * omega * dv_dt.dot(dv.transpose());
            }
            
            if(feasibility_a(a, da, ca)){
                j_ += step * omega * ca;
                dj_dq_.block<4,3>(i,0) += step * omega * B3.transpose() * da;
                const Eigen::RowVector3d da_dt = (Eigen::RowVector4d{0, 0, -4, -18 * s1} / t3 * MQ );
                dj_dt_ += step * omega * da_dt.dot(da.transpose());
            }

            s1 += step;
        }
    }
}

inline void collision_cost(){
    double costp = 0;
    // for(int i = 0;i < urbs_.inner_points_count_ + 2;++i){
    //     cps.col(i + 1) = (trajectory::UniformBSpline::P0 * urbs_.Q_.block<4,3>(1 + i,0)).transpose();
    // }
    // ego_.collision_cost(cps, gradp, costp);
    // j_ += omega_safety_ * costp;
    // for(int i = 1;i < urbs_.inner_points_count_ + 3;++i){
    //     dj_dq_.block<4,3>(i,0) += omega_safety_ * trajectory::UniformBSpline::P0.transpose() * gradp.col(i).transpose();
    // }
    
    for(int i = 0;i < urbs_.inner_points_count_ + 2;++i){
        cps.col(i + 1) = (trajectory::UniformBSpline::P0 * urbs_.Q_.block<4,3>(1 + i,0)).transpose();
    }
    ego_.collision_cost(cps, gradp, costp);
    j_ += omega_safety_ * costp;
    for(int i = 1;i < urbs_.inner_points_count_ + 3;++i){
        dj_dq_.block<4,3>(i,0) += omega_safety_ * trajectory::UniformBSpline::P0.transpose() * gradp.col(i).transpose();
    }
    // std::cerr << omega_safety_ * costp << std::endl;
    
}

inline bool feasibility_v(const Eigen::RowVector3d& v, Eigen::RowVector3d& gradv,double& costv){
    double v_squnorm = v.squaredNorm();
    double vpen = v_squnorm - max_vel_square;
    if (vpen > 0)
    {
      gradv = omega_feasib_ * 6 * vpen * vpen * v;
      costv = omega_feasib_ * vpen * vpen * vpen;
      return true;
    }
    return false;
}
inline bool feasibility_a(const Eigen::RowVector3d& a, Eigen::RowVector3d& grada,double& costa){
    double a_squnorm = a.squaredNorm();
    double apen = a_squnorm - max_acc_square;
    if (apen > 0)
    {
      grada = omega_feasib_ * 6 * apen * apen * a;
      costa = omega_feasib_ * apen * apen * apen;
      return true;
    }
    return false;
}

Eigen::Matrix3Xd cps;
Eigen::Matrix3Xd gradp;

const esdf::ImplicitESDF<TMap>&         esdf_;
const TMap&                             map_;
trajectory::impl::UniformBSplineImpl&   urbs_;
Ego<TMap> ego_;

int iter_num = 0;
const double max_vel_square = 1;
const double max_acc_square = 1;

const double omega_feasib_ = 10;
const double omega_smooth_ = 1e-2;
const double omega_safety_ = 1e4;

const int K = 10;
static const inline Eigen::Matrix4d JTJ = (trajectory::UniformBSpline::J.transpose() * trajectory::UniformBSpline::J).eval(); 
};
}
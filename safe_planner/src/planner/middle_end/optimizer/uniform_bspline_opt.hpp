#pragma once

#include "i_optimizable.hpp"
#include "../trajectory/uniform_bspline_impl.hpp"
#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/map/imap.hpp"
#include "safe_planner/trajectory/bspline.hpp"

#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>

namespace safe_planner::planner::middle_end::optimizer {
template<class TMap>
class UniformBSplineOpt final : public IOptimizeable{
public:
UniformBSplineOpt(
    trajectory::impl::UniformBSplineImpl& urbs,
    const esdf::ImplicitESDF<TMap>& esdf,
    const TMap& map
)
: urbs_(urbs)
, esdf_(esdf)
, map_(map)
{
    dj_dq_.resize(urbs.Q_.rows(),urbs.Q_.cols());
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

void reset(){
    j_ = 0;
    dj_dq_.setZero();
    dj_dt_ = 0;
    minimal_jerk();
    samples_cost();
}

double j_;
Eigen::MatrixX3d dj_dq_;
double dj_dt_;
private:
void minimal_jerk()
{
    const double t = urbs_.t_;
    const double t3 = t * t * t;
    const double t6 = t3 * t3;
    const auto J = trajectory::UniformBSpline::J / t3;
    const auto jTj = JTJ / t6;

    for(int i = 0;i < urbs_.segment_count_;++i){
        auto Q = urbs_.Q_.block<4,3>(i,0);
        const auto j =  omega_smooth_ * (J * Q).squaredNorm();
        j_ += j;

        dj_dq_.block<4,3>(i,0) += (2 * omega_smooth_ * jTj * Q);
        dj_dt_ += omega_smooth_ * -6 * j / t; 
    }
}
void samples_cost(){
    const auto t1 = urbs_.t_;
    const auto t2 = t1 * t1;
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
            const auto s3 = s2 * s1;

            const auto omega = 
                ((j == 0 || j == K) && (i != 0 && i != l)) 
                ? 0.5 : 1.0;
            const Eigen::RowVector4d B1 = 
                Eigen::RowVector4d{1, s1, s2, s3} * trajectory::UniformBSpline::M;
            const Eigen::RowVector4d B2 = 
                Eigen::RowVector4d{0, 1, 2 * s1,3 * s2} / t1 * trajectory::UniformBSpline::M;
            const Eigen::RowVector4d B3 = 
                Eigen::RowVector4d{0, 0, 2, 6 * s1} / t2 * trajectory::UniformBSpline::M;
            const auto p = (B1 * Q).cast<float>();
            const auto v = (B2 * Q);
            const auto a = (B3 * Q);

            if(map_.check_point_d(p) != IMap::State::Safe){
                const auto [l,g] = esdf_.get_esdf(p, 0);
                j_ += omega_safety_ * omega * step * l * l;
                dj_dq_.block<4,3>(i,0) += 2 * omega_safety_ * omega * step * (l * g.template cast<double>() * B1).transpose();
            }

            if(feasibility_v(v, dv, cv)){
                j_ += step * omega * cv;
                dj_dq_.block<4,3>(i,0) += step * omega * B2.transpose() * dv;
                const auto dv_dt = (step * omega * -1 * B2 * Q / t1);
                dj_dt_ += dv_dt.dot(dv.transpose());
            }
            
            if(feasibility_a(a, da, ca)){
                j_ += step * omega * ca;
                dj_dq_.block<4,3>(i,0) += step * omega * B3.transpose() * da;
                const auto da_dt = (step * omega * -2 * B3 * Q / t1);
                dj_dt_ += da_dt.dot(da.transpose());
            }

            s1 += step;
        }
    }
}

bool feasibility_v(const Eigen::RowVector3d& v, Eigen::RowVector3d& gradv,double& costv){
    double v_squnorm = v.squaredNorm();
    double vpen = v_squnorm - max_vel_square;
    if (vpen > 0)
    {
      gradv = wei_feas_mod_ * 6 * vpen * vpen * v;
      costv = wei_feas_mod_ * vpen * vpen * vpen;
      return true;
    }
    return false;
}
bool feasibility_a(const Eigen::RowVector3d& a, Eigen::RowVector3d& grada,double& costa){
    double a_squnorm = a.squaredNorm();
    double apen = a_squnorm - max_acc_square;
    if (apen > 0)
    {
      grada = wei_feas_mod_ * 6 * apen * apen * a;
      costa = wei_feas_mod_ * apen * apen * apen;
      return true;
    }
    return false;
}

trajectory::impl::UniformBSplineImpl&   urbs_;
const esdf::ImplicitESDF<TMap>&         esdf_;
const TMap&                             map_;

const double omega_smooth_ = 1;
const double omega_safety_ = 0;

const double max_vel_square = 9;
const double max_acc_square = 9;

const double wei_feas_mod_ = 10;

const int K = 10;
static const inline auto JTJ = trajectory::UniformBSpline::J.transpose() * trajectory::UniformBSpline::J; 
};
}
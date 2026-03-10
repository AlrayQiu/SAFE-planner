#pragma once

#include "i_optimizable.hpp"
#include "safe_planner/esdf/implicit_esdf.hpp"
#include "../trajectory/polynomial_impl.hpp"
#include "safe_planner/map/imap.hpp"
#include <stdexcept>

namespace safe_planner::planner::middle_end::optimizer {
template<class TMap>
class PolynomialOpt final : public IOptimizeable{
public:
PolynomialOpt(
    trajectory::impl::PolynomialImpl& poly,
    const esdf::ImplicitESDF<TMap>& esdf,
    const TMap& map
)
: poly_(poly)
, esdf_(esdf)
, map_(map)
{
    dj_dc_.resize(poly_.b_.rows(),poly_.b_.cols());
    dj_dt_.resize(poly_.t1_.size());
}

inline void get_j(double&) const{
    throw std::runtime_error("No impl");
}
inline void get_g(Eigen::VectorXd&) const{throw std::runtime_error("No impl");};
inline void set_x(const Eigen::VectorXd&){
throw std::runtime_error("No impl");
}
inline void init_x(Eigen::VectorXd&) const{
throw std::runtime_error("No impl");
};

inline void reset(){
    j_ = 0;
    dj_dc_.setZero();
    dj_dt_.setZero();

    minimal_jerk();
    j_ *= 1e-4;
    dj_dc_.array() *= 1e-4;
    dj_dt_.array() *= 1e-4;
    samples_cost();
}
inline void minimal_jerk(){
    for(int i = 0;i < poly_.pieces_num_;++i){
        j_ += 36.0 * poly_.b_.row(6 * i + 3).squaredNorm() * poly_.t1_(i)
            + 144.0 * poly_.b_.row(6 * i + 3).dot(poly_.b_.row(6 * i + 4)) * poly_.t2_(i) 
            + 240.0 * poly_.b_.row(6 * i + 3).dot(poly_.b_.row(6 * i + 5)) * poly_.t3_(i)
            + 192.0 * poly_.b_.row(6 * i + 4).squaredNorm() * poly_.t3_(i)
            + 720.0 * poly_.b_.row(6 * i + 4).dot(poly_.b_.row(6 * i + 5)) * poly_.t4_(i)
            + 720.0 * poly_.b_.row(6 * i + 5).squaredNorm() * poly_.t5_(i);
        
        dj_dc_.row(i * 6 + 3) += 72.0 * poly_.b_.row(6 * i + 3) * poly_.t1_(i) 
                               + 144.0 * poly_.b_.row(6 * i + 4) * poly_.t2_(i)
                               + 240.0 * poly_.b_.row(6 * i + 5) * poly_.t3_(i);

        
        dj_dc_.row(i * 6 + 4) += 144.0 * poly_.b_.row(6 * i + 3) * poly_.t2_(i) 
                               + 384.0 * poly_.b_.row(6 * i + 4) * poly_.t3_(i)
                               + 720.0 * poly_.b_.row(6 * i + 5) * poly_.t4_(i);

        dj_dc_.row(i * 6 + 5) += 240.0 * poly_.b_.row(6 * i + 3) * poly_.t3_(i) 
                               + 720.0 * poly_.b_.row(6 * i + 4) * poly_.t4_(i)
                               + 1440.0 * poly_.b_.row(6 * i + 5) * poly_.t5_(i);

        dj_dt_(i) += 36.0 * poly_.b_.row(6 * i + 3).squaredNorm()
                   + 288.0 * poly_.b_.row(6 * i + 3).dot(poly_.b_.row(6 * i + 4)) * poly_.t1_(i) 
                   + 720.0 * poly_.b_.row(6 * i + 3).dot(poly_.b_.row(6 * i + 5)) * poly_.t2_(i)
                   + 576.0 * poly_.b_.row(6 * i + 4).squaredNorm() * poly_.t2_(i)
                   + 2880.0 * poly_.b_.row(6 * i + 4).dot(poly_.b_.row(6 * i + 5)) * poly_.t3_(i)
                   + 3600.0 * poly_.b_.row(6 * 4 + 5).squaredNorm() * poly_.t4_(i); 
        
    }
}

inline void samples_cost(){
    Eigen::Matrix3Xd points;
    points.resize(3, poly_.pieces_num_ * sample_count_ + 1);
    int id_p = -1;
    Eigen::Vector<double,6> beta0{};
    Eigen::Vector<double,6> beta1{};
    Eigen::Vector<double,6> beta2{};
    Eigen::Vector<double,6> beta3{};
    Eigen::Vector3d gradp;
    Eigen::Vector3d gradv;
    Eigen::Vector3d grada;
    double costp;
    double costv;
    double costa;
    for(int i = 0;i < poly_.pieces_num_;++i){
        double s1 = 0;
        double step = poly_.t1_(i) / sample_count_;
        auto ci = poly_.b_.block<6,3>(i * 6, 0).transpose();
        for(int j = 0;j <= sample_count_;++j){
            if(j != sample_count_ || (j == sample_count_ && i== poly_.pieces_num_ - 1))
                ++id_p;
            auto s2 = s1 * s1;
            auto s3 = s2 * s1;
            auto s4 = s3 * s1;
            auto s5 = s4 * s1;
            beta0 << 1,s1,s2,s3,s4,s5;
            beta1 << 0,1,2 * s1,3 * s2,4 * s3,5 * s4;
            beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3; 
            beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
            auto alpha = 1.0 * j / sample_count_;
            auto pos = ci * beta0;
            auto vel = ci * beta1;
            auto acc = ci * beta2;
            auto jer = ci * beta3;
            const auto omg = (j == 0 || j == sample_count_) ? 0.5 : 1.0;
            points.col(id_p) = pos;

            if(obstacle_cost(pos, gradp, costp))
            {
                auto gradViolaPc = beta0 * gradp.transpose();  
                auto gradViolaPt = (alpha * gradp.transpose() * vel).trace(); 
                dj_dc_.block<6,3>(i* 6, 0) += omg * step * gradViolaPc;
                dj_dt_(i) += omg * (costp / sample_count_ + step * gradViolaPt);
                j_ += omg * step * costp;
            }
            if (feasibility_cost_v(vel, gradv, costv))
            {
                auto gradViolaVc = beta1 * gradv.transpose();      
                auto gradViolaVt = (alpha * gradv.transpose() * acc).trace();   
                dj_dc_.block<6, 3>(i * 6, 0) += omg * step * gradViolaVc;   
                dj_dt_(i) += omg * (costv / sample_count_ + step * gradViolaVt); 
                j_ += omg * step * costv;
            }

            if (feasibility_cost_a(acc, grada, costa))
            {
                auto gradViolaAc = beta2 * grada.transpose();
                auto gradViolaAt = (alpha * grada.transpose() * jer).trace();    
                dj_dc_.block<6, 3>(i * 6, 0) += omg * step * gradViolaAc; 
                dj_dt_(i) += omg * (costa / sample_count_ + step * gradViolaAt); 
                j_ += omg * step * costa;
            }

        }
    }
    Eigen::Matrix3Xd gdp{};
    double cost;
    distance_cost(points,gdp,cost);
    j_ += cost;
    id_p = -1;
    for (int i = 0; i < poly_.pieces_num_; ++i) // N control point. N pieces
    {
        auto step = poly_.t1_(i) / sample_count_;  
        auto s1 = 0.0; 
        
        for (int j = 0; j <= sample_count_; ++j)
        {
            if (j != sample_count_ || (j == sample_count_ && i == poly_.pieces_num_ - 1)) 
            {
                ++id_p; // a control point(MINCO mid point) has K finely check points
            }
            auto s2 = s1 * s1;   //t^2
            auto s3 = s2 * s1;   //t^3
            auto s4 = s2 * s2;   //t^4
            auto s5 = s4 * s1;   //t^5
            beta0 << 1.0, s1, s2, s3, s4, s5;
            beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
            auto alpha = 1.0 / sample_count_ * j;
            auto vel = poly_.b_.block<6, 3>(i * 6, 0).transpose() * beta1;

            auto omg = (j == 0 || j == sample_count_) ? 0.5 : 1.0;
            auto gradViolaPc = beta0 * gdp.col(id_p).transpose();       // \frac{\partial k}{\partial c}
            auto gradViolaPt = (alpha * gdp.col(id_p).transpose() * vel).trace(); // \frac{\partial k}{\partial T}
            dj_dc_.block<6, 3>(i * 6, 0) += omg * gradViolaPc;
            dj_dt_(i) += omg * (gradViolaPt);
            s1 += step;
        }
    }
}  
inline void distance_cost(const Eigen::Matrix3Xd &ps,
                    Eigen::Matrix3Xd &gdp,
                    double &var)
{
    int N = ps.cols() - 1;  
    Eigen::MatrixXd dps = ps.rightCols(N) - ps.leftCols(N); // delta
    Eigen::VectorXd dsqrs = dps.colwise().squaredNorm().transpose();
    double dquarsum = dsqrs.squaredNorm();
    double dquarmean = dquarsum / N;
    var = wei_sqrvar_ * (dquarmean);
    gdp.resize(3, N + 1);
    gdp.setZero();
    for (int i = 0; i <= N; i++)
    {
        if (i != 0)
        { 
            gdp.col(i) += wei_sqrvar_ * (4.0 * (dsqrs(i - 1)) / N * dps.col(i - 1));
        }
        if (i != N)
        {
            gdp.col(i) += wei_sqrvar_ * (-4.0 * (dsqrs(i)) / N * dps.col(i));
        }
    }
    return;
}


inline bool obstacle_cost(const Eigen::Vector3d& p, Eigen::Vector3d& dp, double& cost){
    if(map_.check_point_d(p.template cast<float>()) == IMap::State::Safe) return false;

    const auto[d,g] = esdf_.get_esdf(p.cast<float>(), 0);
    dp = 2 * g.template cast<double>() * d;
    cost = d * d;
    return true; 
}

inline bool feasibility_cost_v(const Eigen::Vector3d &v,
                        Eigen::Vector3d &gradv,
                        double &costv){
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
inline bool feasibility_cost_a(const Eigen::Vector3d &a,
                        Eigen::Vector3d &grada,
                        double &costa)
{
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
const double max_vel_square = 9;
const double max_acc_square = 9;
const double wei_feas_mod_ = 0.1;
const double wei_sqrvar_ = 0;
const int sample_count_ = 100;
Eigen::MatrixX3d dj_dc_;
Eigen::VectorXd  dj_dt_;
double j_;
trajectory::impl::PolynomialImpl& poly_;
const esdf::ImplicitESDF<TMap>& esdf_;
const TMap& map_;
};

}
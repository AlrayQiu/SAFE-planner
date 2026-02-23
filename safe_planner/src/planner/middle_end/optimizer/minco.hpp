#pragma once
#include "i_optimizable.hpp"
#include "safe_planner/trajectory/minco.hpp"
#include "safe_planner/trajectory/utils/math.hpp"

namespace  safe_planner::planner::optimizable {
class MincoOpt final : public IOptimizeable{
public:
MincoOpt() = delete;
MincoOpt(trajectory::Minco& minco)
: minco_(minco)
, time_min_(*minco.times.begin())
, time_err_(*minco.times.rbegin() - time_min_)
, dJ_dc_(minco.params.rows(), minco.params.cols())
, dJ_dt_(minco.times.size() - 1){};
inline void clear_j_g(){
    j_ = 0;
    dJ_dc_.setZero();
    dJ_dt_.setZero();
    const auto& times = minco_.times;
    const auto& param = minco_.params;
    time_min_ = *times.begin(); 
    time_err_ = *times.rbegin() - time_min_; 
    const auto M = minco_.times.size() - 1;
    minimal_snap(M,param,times);
    dynamic(0.1,param,times);
};
void add_j(const float j){j_ += j;}
void add_g(const float ratio, const int i, const Eigen::Vector3f& dk_dp){
    const auto t = ratio * time_err_ + time_min_;
    const auto index = minco_.times[i];
    Eigen::Vector<float, 6> b;
    auto tc = t.count(),ic = index.count(),
    alpha = (tc - ic) / ((minco_.times[i + 1].count() - ic));
    trajectory::utils::power_base<6>((t.count() - index.count()), 0, b);
    dJ_dc_.block<6,3>(i * 6, 0) += b * (dk_dp.transpose());
    Eigen::RowVector<float, 6> b1;
    for(int j = 0;j < 6;++j){
        b1(j) = b(j) * j;
    }
    dJ_dt_(i) += alpha * ((b1 * minco_.params.block<6,3>(i * 6,0)).dot(dk_dp));
}
void enum_p(const float ratio, std::function<void(const int, const float,const Eigen::Vector3f&)> foo){
    auto time = time_min_.count();
    auto i_t = minco_.times.begin(); 
    int i = 0;
    float t2 = (i_t + 1)->count() - time;
    while(i_t + 1 != minco_.times.end()){
        Eigen::RowVector<float, 6> B;
        const auto t1 = time - i_t->count();
        trajectory::utils::power_base<6>(t1, 0, B);
        foo(i, t1 / t2, (B * minco_.params.block<6,3>(i*6,0)).transpose());
        time += ratio;
        if(time < (i_t + 1)->count()) continue;
        ++i_t;
        ++i;
        t2 = i_t->count() - (i_t - 1)->count();
    }
}
void get_j(float& j) const {j = j_ + time_wight * time_err_.count() ;}
void get_g(const Eigen::VectorXf& x,Eigen::VectorXf& g) const{
    const auto M = minco_.times.size() - 1;
    const auto L = (M - 1) * 3;
    const auto m = minco_.get_M().transpose();
    const Eigen::MatrixXf G = m.lu().solve(dJ_dc_);
    for(size_t i = 1; i < M; i++){
        g.block<3,1>((i - 1) * 3, 0) = G.block<1,3>(i * 6 - 3, 0).transpose();
        g(L + i - 1) = dJ_dt_(i - 1) - tr_gg_dei_dt_ci(i, G);
        g(L + i - 1) = (g(L + i - 1) + time_wight) * dt_dtau(x(L + i - 1));
    }
    g(L + M - 1) = dJ_dt_(M - 1) - tr_gg_dem_dt_ci(M, G);
    g(L + M - 1) = (g(L + M - 1) + time_wight) * dt_dtau(x(L + M - 1));
}
void set_x(const Eigen::VectorXf& x) const{
    const auto M = minco_.times.size() - 1;
    const auto L = (M - 1) * 3;
    Eigen::VectorXf xcpy = x;
    for(auto i = L;i < L + M;i++)
        tau_to_t(x(i), xcpy(i));
    minco_.set_param(xcpy);
};
void init_x(Eigen::VectorXf& x) const{
    const int M = minco_.times.size() - 2;
    const int L = (minco_.times.size() - 2) * 3;
    x.setZero(L + M + 1);
    for(int i = 0;i < M;i++){
        x.block<3,1>(i*3, 0) = minco_.params.row((i + 1) * 6);
    }
    float tau;
    for(int i = 0;i < minco_.times.size() - 1;i++){
        t_to_tau((minco_.times[i + 1] - minco_.times[i]).count() ,tau);
        x(L + i) = tau;
    }
};

private:
inline void tau_to_t(const float& tau,float& t) const{
    if(tau > 0)
        t = tau * tau * 0.5f + tau + 1;
    else 
        t = 1 / (tau * tau * 0.5f - tau + 1);
}
inline void t_to_tau(const float& t,float& tau) const{
    if(t > 1)
        tau = std::sqrt(2 * t - 1) - 1;
    else
        tau = 1 - std::sqrt(2 / t - 1);
}
inline float dt_dtau(const float tau) const{
    if(tau > 0)
        return tau + 1;
    else 
    {
        auto den_sqrt = (0.5f * tau - 1.0f) * tau + 1.0;
        return (1 - tau) / (den_sqrt * den_sqrt);
    }
}
inline void calcu_minimalsnap_j(const auto& c3,const auto& c4, const auto& c5, const Eigen::Vector<float,5> t, float& j)
{
    j +=    36 * c3.squaredNorm() * t(0) + 
            144 * c3.dot(c4) * t(1) +
            240 * c3.dot(c5) * t(2) +
            192 * c4.squaredNorm() * t(2) +
            720 * c4.dot(c5) * t(3) +
            720 * c5.squaredNorm() * t(4);
}
inline void minimalsnap_dj_dc(const auto& c3,const auto& c4, const auto& c5, const Eigen::Vector<float,5> t, const int id, auto& dj_dc){
    dj_dc.row(id) +=        72 * c3 * t(0) + 
                            144 * c4 * t(1) + 
                            240 * c5 * t(2);
    dj_dc.row(id + 1) +=    144 * c3 * t(1) +
                            384 * c4 * t(2) +
                            720 * c5 * t(3);
    dj_dc.row(id + 2) +=    240 * c3 * t(2) +
                            720 * c4 * t(3) +
                            1440 * c5 * t(4);
}
inline void minimal_snap_dj_dt(const auto& c3,const auto& c4, const auto& c5, const Eigen::Vector<float,5> t, float& dj_dt) const{
    
    dj_dt +=    36 * c3.squaredNorm() + 
                288 * c3.dot(c4) * t(0) +
                720 * c3.dot(c5) * t(1) +
                576 * c4.squaredNorm() * t(1) +
                2880 * c4.dot(c5) * t(2) +
                3600 * c5.squaredNorm() * t(3);
}

inline float tr_gg_dei_dt_ci(const int i, const auto& G) const{
    Eigen::Vector<float,6> b;
    auto curr = i * 6 - 3;
    const auto t = (minco_.times[i] - minco_.times[i - 1]).count();
    trajectory::utils::power_base<6>(t, 1, b);
    Eigen::Matrix<float,6,6> Ei{};
    for(int j = 0;j < 6;++j){
        if(t == 0)
        {
            Ei(0,j) = Ei(1,j) = b(j);
            Ei(2,j) = Ei(1,j) * (j - 1);
            Ei(3,j) = Ei(2,j) * (j - 2);
            Ei(4,j) = Ei(3,j) * (j - 3);
            Ei(5,j) = Ei(4,j) * (j - 4);
        } else{
            Ei(0,j) = Ei(1,j) = b(j);
            Ei(2,j) = Ei(1,j) * (j - 1) / t;
            Ei(3,j) = Ei(2,j) * (j - 2) / t;
            Ei(4,j) = Ei(3,j) * (j - 3) / t;
            Ei(5,j) = Ei(4,j) * (j - 4) / t;
        }
    }
    return (G.template block<6,3>(curr, 0).transpose() * Ei * minco_.params.block<6,3>((i - 1) * 6,0)).trace();
}
inline float tr_gg_dem_dt_ci(const int m, const auto& G) const{
    Eigen::Vector<float,6> b;
    auto curr = m * 6 - 3;
    const auto t = (minco_.times[m] - minco_.times[m - 1]).count();
    trajectory::utils::power_base<6>(t, 1, b);
    Eigen::Matrix<float,3,6> Em{};
    Em.setZero();
    for(int j = 0;j < 6;++j){
        if(t == 0)
        {
            Em(0,j) = b(j);
            Em(1,j) = Em(0,j) * (j - 1);
            Em(2,j) = Em(1,j) * (j - 2);
        }else{
            Em(0,j) = b(j);
            Em(1,j) = Em(0,j) * (j - 1) / t;
            Em(2,j) = Em(1,j) * (j - 2) / t;
        }
    }
    return (G.template block<3,3>(curr, 0).transpose() * Em * minco_.params.block<6,3>((m - 1) * 6,0)).trace();
}

inline int index(const int row, const int col) const{
    return row * 3 + col;
}

inline void minimal_snap(int M, const Eigen::MatrixX3f& param, const auto& times){
    for(int i = 0;i < M;i++){
        const auto id3 = i * 6 + 3;
        const auto id4 = i * 6 + 4;
        const auto id5 = i * 6 + 5;
        const auto t = times[i + 1].count() - times[i].count();
        const Eigen::Vector<float,5> t_vec{t, t * t, t * t * t, t * t * t * t, t * t * t * t * t};
        calcu_minimalsnap_j(param.row(id3), param.row(id4), param.row(id5), t_vec, j_);
        minimalsnap_dj_dc(param.row(id3), param.row(id4), param.row(id5), t_vec, id3, dJ_dc_);
        minimal_snap_dj_dt(param.row(id3), param.row(id4), param.row(id5), t_vec, dJ_dt_(i));
    }
}
inline void dynamic(const float rate, const auto& param, const auto& times){
    auto time = time_min_.count();
    auto i_t = times.begin(); 
    while(i_t + 1 != times.end()){
        const float t = (time - i_t->count());
        const float ratio = t / ((i_t + 1)->count() - i_t->count());
        const int i = std::distance(times.begin(), i_t);
        static Eigen::RowVector<float, 6> b1{};
        static Eigen::RowVector<float, 6> b2{};
        static Eigen::RowVector<float, 6> b3{};
        trajectory::utils::power_base<6>(t, 1, b1);
        for(int j = 0;j < 6;++j){
            if(t == 0){
                b2(j) = b1(j) * (j - 1);
                b3(j) = b2(j) * (j - 2);    
            }else{
                b2(j) = b1(j) * (j - 1) / t;
                b3(j) = b2(j) * (j - 2) / t;
            }
        }
        const auto& ci = param.template block<6,3>(i * 6,0);
        const Eigen::RowVector3f v = (b1 * ci);
        const Eigen::RowVector3f a = (b2 * ci);
        dynamic_v(i, ratio, v, b1, b2, ci);
        dynamic_a(i, ratio, a, b2, b3, ci);
        time += rate;
        if(time > (i_t + 1)->count()) ++i_t;
    }
}
inline void dynamic_v(
    const int i, const float alpha,
    const Eigen::RowVector3f& v,
    const Eigen::RowVector<float,6>& b1,
    const Eigen::RowVector<float,6>& b2,
    const auto& c_i){
    const auto l2_v = v.squaredNorm();
    if(l2_v <= vmax_square) return;
    const auto err_v = 2 * (l2_v - vmax_square);
    dJ_dc_.block<6,3>(i * 6, 0) += err_v * b1.transpose() * v;
    dJ_dt_(i) += alpha * err_v * b2 * c_i * v.transpose(); 
}
inline void dynamic_a(
    const int i, const float alpha,
    const Eigen::RowVector3f& a,
    const Eigen::RowVector<float,6>& b2,
    const Eigen::RowVector<float,6>& b3,
    const auto& c_i){
    const auto l2_a = a.squaredNorm();
    if(l2_a <= amax_square) return;
    const auto err_v = 2 * (l2_a - vmax_square);
    dJ_dc_.block<6,3>(i * 6, 0) += err_v * b2.transpose() * a;
    dJ_dt_(i) += alpha * err_v * b3 * c_i * a.transpose(); 
}

float j_ = 0;
trajectory::Minco& minco_;
trajectory::Minco::seconds time_min_;
trajectory::Minco::seconds time_err_;
Eigen::MatrixX3f dJ_dc_{};
Eigen::VectorXf dJ_dt_{};


const int vmax_square = 9;
const int amax_square = 9;
const float time_wight = 20;
};
}
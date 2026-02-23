#pragma once

#include "trajectroies.hpp"
#include "utils/math.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <chrono>
#include <iterator>
#include <vector>

namespace safe_planner::planner::trajectory {
class Minco final : public ITrajectory{
public:

Minco() = default;
inline Eigen::Vector3f pos(const seconds s) const{
    if(s >= *times.rbegin()) {
        Eigen::RowVector<float, 6> B;
        auto t = times.rbegin()->count() - (times.rbegin() + 1) -> count();
        trajectory::utils::power_base<6>(t, 0, B);
        return (B * params.block<6,3>((M_ - 1)*6,0)).transpose();
    }
    Eigen::RowVector<float, 6> B;
    auto t = std::upper_bound(times.begin(), times.end(), s) - 1;
    auto time = s.count() - t->count();
    if(t >= times.end()) t = times.end() - 1;
    utils::power_base<6>(time, 0, B);
    auto i = std::distance(times.begin(),t);
    // std::cerr << s <<  times[i] << std::endl;
    return (B * params.block<6,3>(i*6,0)).transpose();
};
inline std::vector<Eigen::Vector3f> control_point(std::vector<Eigen::Vector3f>& s) const{
    s.clear();
    for(int i = 1; i < M_; i++){
        s.push_back(params.block<6,3>(i*6,0).row(0));
    }
    return s;
};
void get_time_range(seconds& from, seconds& to) const
{
    from = seconds{0};
    to = *times.rbegin();
};

const Eigen::MatrixXf& get_M() const{
    return m_;
}

void set_param(const Eigen::VectorXf& control){
    const auto L = (M_ - 1) * s();
    for(int i = 1; i < M_; i++){
        b_.block<1,3>(i * s2() - s(),0) = control.block<3,1>((i - 1) * s(),0).transpose();
        times[i] = std::chrono::duration<float>(control(L + i - 1)) + times[i - 1];
    }
    times[M_] = std::chrono::duration<float>(control(L + M_ - 1)) + times[M_ - 1];
    build_m();
    params = m_.lu().solve(b_);
}


Minco(
    const std::vector<Eigen::Vector3f>& control_points, 
    const std::vector<seconds>& control_seconds,
    const Eigen::Matrix<float,2,3>& begin_va,
    const Eigen::Matrix<float,2,3>& end_va)
    : times(control_seconds)
    , M_(control_points.size() - 1)
{
    const auto L = M_ * s2();
    b_.setZero(L,3);
    m_.setZero(L,L);
    build_m();
    b_.block<1,3>(0,0) = control_points[0].transpose();
    b_.block<2,3>(1,0) = begin_va;
    for(int i = 1; i < M_; i++){
        b_.block<1,3>(i * s2() - s(),0) = control_points[i].transpose();
    }
    b_.block<1,3>(L - s(),0) = control_points[M_].transpose();
    b_.block<2,3>(L - s() + 1,0) = end_va;
    params.resize(M_ * s2(),3);
    params = m_.lu().solve(b_);
}

Eigen::MatrixX3f params;
std::vector<seconds> times;
private:
using VecP = Eigen::Vector<float, 6>;
using MatP = Eigen::Matrix<float, 3, 6>;

#pragma region M矩阵构建
void build_m(){
    // F
    m_(0, 0) = 1;
    m_(1, 1) = 1;
    m_(2, 2) = 2;
    
    VecP b;
    
    // Ei Fi
    for(int i = 1; i < M_; ++i){
        auto curr_row = i * s2() - s();
        auto crr_col = s2() * (i - 1);
        const auto t = (times[i] - times[i - 1]).count();
        utils::power_base<6>(t, 0, b);
        for(int j = 0;j < s2();++j){
            auto b_0_j = b(j);
            auto b_1_j = b_0_j * j;
            auto b_2_j = b_1_j * (j - 1);
            auto b_3_j = b_2_j * (j - 2);
            auto b_4_j = b_3_j * (j - 3);
            if(t!=0){
                b_1_j /= t;
                b_2_j /= t;
                b_3_j /= t;
                b_4_j /= t;
            }
            m_(curr_row,      crr_col + j) = b_0_j;
            m_(curr_row + 1,  crr_col + j) = b_0_j;
            
            m_(curr_row + 2, crr_col + j) = b_1_j;
            m_(curr_row + 3, crr_col + j) = b_2_j;
            m_(curr_row + 4, crr_col + j) = b_3_j;
            m_(curr_row + 5, crr_col + j) = b_4_j;
        }
        insert_f(m_,curr_row,s2() * i);
    }
    // EM
    const auto tm = (times[M_] - times[M_ - 1]).count();
    utils::power_base<6>(tm, 0, b);
    for(int j = 0;j < s2();++j){
        const auto curr = M_ * s2() - s();
        const auto b_0_j = b(j);
        const auto b_1_j = b_0_j * j / tm;
        const auto b_2_j = b_1_j * (j - 1) / tm;
        m_(curr + 0, (M_ - 1) * s2() + j) = b_0_j;
        m_(curr + 1, (M_ - 1) * s2() + j) = b_1_j;
        m_(curr + 2, (M_ - 1) * s2() + j) = b_2_j;
    }
}
inline void insert_f(Eigen::MatrixXf& m,const int row,const int col){
    m(row + 1, col + 0) = -1;
    m(row + 2, col + 1) = -1;
    m(row + 3, col + 2) = -2;
    m(row + 4, col + 3) = -6;
    m(row + 5, col + 4) = -24;
}
#pragma endregion

constexpr inline auto size_(){return times.size();}
constexpr inline int  s(){return 3;}
constexpr inline int  s2(){return s() * 2;}

int M_;
Eigen::MatrixXf m_;
Eigen::MatrixX3f b_;
};
}
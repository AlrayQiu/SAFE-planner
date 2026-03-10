#pragma once
#include <Eigen/Eigen>

namespace safe_planner::planner::trajectory::utils
{
    constexpr static inline long long factorial(int n){
        long long r = 1;
        for(int i = 2; i <= n; ++i) r *= i;
        return r;
    }

    constexpr static inline double binom(int n, int k){
        if(k < 0 || k > n) return 0;
        double r = 1;
        for(int i = 1; i <= k; ++i)
            r = r * (n - (k - i)) / i;
        return r;
    }

    // 计算 M(k) 的单个元素 m(i,j)
    constexpr static inline double M_element(int k, int i, int j){
        double C1 = binom(k - 1, k - 1 - i);
        double sum = 0;

        for(int s = j; s <= k - 1; ++s){
            double C2 = binom(k, s - j);
            double term = C2;

            // (k - s - 1)^(k - 1 - i)
            double base = k - s - 1;
            double power = 1;
            for(int p = 0; p < (k - 1 - i); ++p)
                power *= base;

            term *= power;
            if((s - j) % 2 == 1) term = -term;

            sum += term;
        }

        return double(C1 * sum) / double(factorial(k - 1));
    }

    // constexpr 生成 M(k)
    template<int k>
    constexpr static inline Eigen::Matrix<double, k, k> make_M(){
        Eigen::Matrix<double, k, k> M;
        for(int i = 0; i < k; ++i)
            for(int j = 0; j < k; ++j)
                M(i,j) = M_element(k, i, j);
        return M;
    }
    

    constexpr static inline double pow_const(double s, int p){
        double r = 1.0;
        for(int i = 0; i < p; ++i) r *= s;
        return r;
    }
    constexpr static inline double nth_derivative_term(int p, int n, double s){
        if(n > p) return 0.0; // 高阶导数为 0
        double coef = 1.0;
        for(int i = 0; i < n; ++i)
            coef *= (p - i);
        return coef * pow_const(s, p - n);
    }
    template<int k>
    constexpr static inline Eigen::RowVector<double, k> make_S(double s, int n){
        Eigen::RowVector<double, k> S;
        for(int p = 0; p < k; ++p)
            S[p] = nth_derivative_term(p, n, s);
        return S;
    }
}
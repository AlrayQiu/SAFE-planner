#pragma once
#include <Eigen/Eigen>
namespace safe_planner::planner::trajectory::impl{
class PolynomialImpl{

public:
void set_t(const Eigen::VectorXd& t){
    t1_ = t;
    t2_ = t1_.cwiseProduct(t);
    t3_ = t2_.cwiseProduct(t);
    t4_ = t3_.cwiseProduct(t);
    t5_ = t4_.cwiseProduct(t);
}
void set_b(const Eigen::MatrixX3d& b){
    b_ = b;
}
void set_pieces_num(const int n){
    pieces_num_ = n;
}

int pieces_num_ = 0;

Eigen::MatrixX3d b_;
Eigen::VectorXd t1_;
Eigen::VectorXd t2_;
Eigen::VectorXd t3_;
Eigen::VectorXd t4_;
Eigen::VectorXd t5_;
};
}
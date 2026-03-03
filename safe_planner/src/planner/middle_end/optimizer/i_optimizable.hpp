
#pragma once
#include <Eigen/Eigen>
namespace safe_planner::planner
{
class IOptimizeable{
public:
virtual ~IOptimizeable() = default;
virtual void get_j(double&) const = 0;
virtual void get_g(Eigen::VectorXd&) const = 0;
virtual void set_x(const Eigen::VectorXd&) = 0;
virtual void init_x(Eigen::VectorXd&) const = 0;
};
}
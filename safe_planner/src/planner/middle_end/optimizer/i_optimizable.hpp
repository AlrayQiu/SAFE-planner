
#pragma once
#include <Eigen/Eigen>
namespace safe_planner::planner
{
class IOptimizeable{
public:
virtual ~IOptimizeable() = default;
virtual void clear_j_g() = 0;
virtual void add_j(const float j) = 0;
virtual void add_g(const float ratio, const int, const Eigen::Vector3f& g) = 0;
virtual void enum_p(const float ratio, std::function<void(const int, const float,const Eigen::Vector3f&)>) = 0;
virtual void get_j(float&) const = 0;
virtual void get_g(const Eigen::VectorXf& x,Eigen::VectorXf&) const = 0;
virtual void set_x(const Eigen::VectorXf&) const = 0;
virtual void init_x(Eigen::VectorXf&) const = 0;
};
}
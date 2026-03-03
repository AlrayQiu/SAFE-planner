#pragma once

#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/map/imap.hpp"
#include "safe_planner/trajectory/bspline.hpp"
#include "safe_planner/trajectory/multi_poly.hpp"
#include <concepts>
#include <memory>
#include <vector>

namespace safe_planner::planner
{
namespace middle_end{
enum class Strategy{

};
}
template<class TMap>
requires std::derived_from<TMap, IMap>
class MiddleEnd{
public:
    MiddleEnd(const TMap&, const esdf::ImplicitESDF<TMap>&);
    MiddleEnd() = delete;
    ~MiddleEnd();

    void optimize(
        const std::vector<Eigen::Vector3d>& ref_path,
        const Eigen::Matrix<double,3,2>& begin_va,
        const Eigen::Matrix<double,3,2>& end_va,
        trajectory::UniformBSpline& traj);
    void test_optimizer(
        trajectory::MultiPoly& traj,
        std::vector<Eigen::Vector3d>& gradients);
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
}
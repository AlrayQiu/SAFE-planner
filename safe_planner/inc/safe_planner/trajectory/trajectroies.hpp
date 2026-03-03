#pragma once

#include <Eigen/Eigen>
#include <chrono>
namespace safe_planner::planner
{
class ITrajectory{
public:
using seconds = std::chrono::duration<double>;
virtual ~ITrajectory() = default;

virtual Eigen::Vector3d pos(const seconds) const = 0;
virtual void get_time_range(seconds& from, seconds& to) const = 0;
};
}

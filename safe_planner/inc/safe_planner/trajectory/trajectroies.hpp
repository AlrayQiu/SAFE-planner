#pragma once

#include <Eigen/Eigen>
#include <chrono>
namespace safe_planner::planner
{
class ITrajectory{
public:
using seconds = std::chrono::duration<float>;
virtual ~ITrajectory() = default;

virtual Eigen::Vector3f pos(const seconds) const = 0;
virtual void get_time_range(seconds& from, seconds& to) const = 0;
};
}

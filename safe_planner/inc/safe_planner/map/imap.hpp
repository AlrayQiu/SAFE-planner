#pragma once

#include <Eigen/Eigen>

namespace safe_planner{

class IMap{
public:
enum class State{
    Safe,
    Unsafe  
};

virtual ~IMap() = default;

virtual State check_line_i(const Eigen::Vector3i& form, const Eigen::Vector3i& to) const= 0;
virtual State check_line_d(const Eigen::Vector3f& form, const Eigen::Vector3f& to) const= 0;
virtual State check_point_i (const Eigen::Vector3i& index) const = 0;
virtual State check_point_d (const Eigen::Vector3f& pos)   const = 0;

virtual void get_map_bound_i(Eigen::Vector3i& min,Eigen::Vector3i& max) const = 0;
virtual void get_map_bound_d(Eigen::Vector3f& min,Eigen::Vector3f& max) const = 0;
};
}
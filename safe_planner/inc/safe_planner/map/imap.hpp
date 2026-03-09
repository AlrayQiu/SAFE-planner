#pragma once

#include <Eigen/Eigen>

namespace safe_planner{


namespace cell {
enum class CenterPosition : bool {
    center_in_cornor = false,
    center_in_center = true
};
}

static inline const cell::CenterPosition center_position = cell::CenterPosition::center_in_center;

class IMap{
public:
enum class State{
    Safe,
    Unsafe  
};

virtual ~IMap() = default;

virtual State check_line_d(const Eigen::Vector3f& form, const Eigen::Vector3f& to) const= 0;
virtual State check_point_d (const Eigen::Vector3f& pos)   const = 0;

virtual void get_map_bound_d(Eigen::Vector3f& min,Eigen::Vector3f& max) const = 0;
};

class IGridMap : public IMap{
public:
virtual ~IGridMap() = default;

virtual void pos_to_grid(const Eigen::Vector3f &pos, Eigen::Vector3f &index) const = 0;
virtual void index_to_pos(const Eigen::Vector3i &index, Eigen::Vector3f& pos) const = 0;
virtual void pos_to_index(const Eigen::Vector3f &index, Eigen::Vector3i& pos) const = 0;
virtual void get_map_bound_i(Eigen::Vector3i& min,Eigen::Vector3i& max) const = 0;
virtual State check_line_i(const Eigen::Vector3i& form, const Eigen::Vector3i& to) const= 0;
virtual State check_point_i (const Eigen::Vector3i& index) const = 0;
virtual float resolution() const = 0;
};
}
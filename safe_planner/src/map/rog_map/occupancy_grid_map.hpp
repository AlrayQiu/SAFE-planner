#pragma once

#include "safe_planner/map/rog_map.hpp"

#include <Eigen/Eigen>
#include <cmath>
#include <cstddef>
#include <pcl/PCLPointField.h>
#include <vector>

namespace safe_planner::map::rog_map {
    
namespace occupancy_grid_map {

enum class GridType{
    UNDEFINED = 0,
    UNKNOWN = 1,
    KNOWN_FREE = 2,
    OCCUPIED = 3
};



struct ProblisticMap final{
    using Cfg = Config::ProblisticConfig;
    inline ProblisticMap(
        const Eigen::Vector3i &half_map_size_i,
        const Cfg & config)
    : occupancy_buffer_{}
    , config_{config}
    , map_size_{half_map_size_i * 2 + Eigen::Vector3i::Constant(1)}
    {
        occupancy_buffer_.resize(static_cast<size_t>(map_size_.prod()));
    }

    inline void reset_local_map(){
        std::fill(occupancy_buffer_.begin(), occupancy_buffer_.end(), 0);
    }
    inline void reset_cell(const int& id){
        occupancy_buffer_[static_cast<size_t>(id)] = 0;
    }

    inline bool is_occupied(const float& log_odds_value) const{
        return log_odds_value >= config_.l_occu;
    }
    inline bool is_free(const float& log_odds_value) const{
        return log_odds_value <= config_.l_free;
    }

    std::vector<float> occupancy_buffer_;
    const Cfg config_;
    const Eigen::Vector3i map_size_;
};
}



}

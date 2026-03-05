#pragma once

#include "occupancy_grid_map.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <Eigen/Eigen>

namespace safe_planner::map::rog_map {

template <cell::CenterPosition center = cell::CenterPosition::center_in_cornor>
class InflationMap {
public:
    inline InflationMap(const Eigen::Vector3i &half_map_size_i, const int inflation_radius)
    : inflation_core_(calculate_inflation_core(inflation_radius))
    {
        auto map_size = half_map_size_i * 2 + Eigen::Vector3i::Constant(1);
        occupancy_buffer_.resize(static_cast<size_t>(map_size.prod()), 0);
    };

    inline void reset_cell(const int id){
        occupancy_buffer_[static_cast<size_t>(id)] = 0;
    };
    inline void reset_local_map(){
        std::fill(occupancy_buffer_.begin(), occupancy_buffer_.end(), 0);
    };

    inline const std::vector<Eigen::Vector3i>& get_local_update_buffer() const noexcept{
        return inflation_core_;
    }

    inline void update_cell(const int id,const occupancy_grid_map::GridType from, const occupancy_grid_map::GridType to){
        if(from == occupancy_grid_map::GridType::OCCUPIED && to == occupancy_grid_map::GridType::KNOWN_FREE)
            --occupancy_buffer_[id];
        else if (to == occupancy_grid_map::GridType::OCCUPIED)
            ++occupancy_buffer_[id];
    }

    inline bool is_occupied(const int index) const{
        return occupancy_buffer_[index] > 0;
    }
    inline bool is_free(const int index){
        return occupancy_buffer_[index] <= 0; 
    }

private:
    inline constexpr std::vector<Eigen::Vector3i> calculate_inflation_core(const int inflation_radius){
        int radius = inflation_radius;
        if constexpr(center == cell::CenterPosition::center_in_center) ++radius;
        
        const int side = radius + inflation_radius; 
        const std::size_t count = static_cast<std::size_t>(side) 
                                * static_cast<std::size_t>(side) 
                                * static_cast<std::size_t>(side);
        std::vector<Eigen::Vector3i> core;
        core.reserve(count);

        for(int i = -inflation_radius; i < radius; i++)
            for(int j = -inflation_radius; j < radius; j++)
                for(int k = -inflation_radius; k < radius; k++)
                    core.emplace_back(i,j,k);

        return core;
    }

    std::vector<int> occupancy_buffer_;
    const std::vector<Eigen::Vector3i> inflation_core_;
};

}

#pragma once

#include <Eigen/Eigen>

#include <map/rog_map.hpp>

namespace safe_planner::planner {

namespace rrt_star {
    enum class SearchAs : bool{
        free = false,
        occupied = true
    };
    struct SearchConfig{
        double step_size = 0.5;
        double goal_sample_rate = 0.1;
        double max_iter = 500;
        SearchAs unknown = SearchAs::free;
        SearchAs out_of_bound = SearchAs::free;
    };
}

template <class T>
requires std::is_same_v<std::decay_t<T>, rrt_star::SearchConfig>
static void rrt_star_search(
    const Eigen::Vector3d   &start,
    const Eigen::Vector3d   &goal,
    const map::ROGMap       &ogm,
    T&& config,
    std::vector<Eigen::Vector3d>  &path_out) {

    path_out.clear();                     
    return;
}
}
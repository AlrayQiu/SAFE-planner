#pragma once

#include <Eigen/Eigen>
#include <type_traits>
namespace safe_planner::map::utils::math {
template <typename T>
static inline constexpr int sign(T&& val) noexcept{
    if(val == std::decay_t<T>(0)) return 0;
    return (val < std::decay_t<T>(0)) ? -1 : 1;
}


static Eigen::Vector3f line_box_intersect_point(const Eigen::Vector3f &pt, const Eigen::Vector3f &pos,
                                    const Eigen::Vector3f &box_min, const Eigen::Vector3f &box_max) {
    Eigen::Vector3f diff = pt - pos;
    Eigen::Vector3f max_tc = box_max - pos;
    Eigen::Vector3f min_tc = box_min - pos;

    float min_t = 1000000;

    for (int i = 0; i < 3; ++i) {
        if (fabs(diff(i)) > 0) continue;

        float t1 = max_tc(i) / diff(i);
        if (t1 > 0 && t1 < min_t)
            min_t = t1;

        float t2 = min_tc(i) / diff(i);
        if (t2 > 0 && t2 < min_t)
            min_t = t2;
    }

    return pos + (min_t - 1e-3) * diff;
}

}

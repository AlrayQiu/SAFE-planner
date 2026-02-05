#pragma once

#include <Eigen/Eigen>
#include <memory>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace safe_planner::map {

namespace rog_map {

struct Config{
    const double virtual_ceil_height{std::numeric_limits<double>::infinity()};
    const double virtual_ground_height{ -std::numeric_limits<double>::infinity()};
    const Eigen::Vector3i half_local_update_box_i{100,100,100};
    const int   downsample_num{1};
    const float intensity_thresh{0};
    const int   update_pcd_batch_size{1};

    const int    inflation_radius{1};
    const float  uva_nearby_zone{0.5};
    const struct ProblisticConfig{
        const float l_miss = -0.27;
        const float l_hit = 0.28;
        const float l_free = -1;
        const float l_occu = 1;
        const float l_min = -2;
        const float l_max = 2;
    } prob_config;
};
}

class ROGMap final
{
public:
    ROGMap() = delete;
    ROGMap(
        const Eigen::Vector3i &half_map_size_i,
        const double &resolution,
        const bool &map_sliding_en,
        const double &sliding_thresh,
        const Eigen::Vector3f &fix_map_origin,
        rog_map::Config &&map_config);
    ~ROGMap();


    void update(const pcl::PointCloud<pcl::PointXYZINormal>& pcd, const Eigen::Vector3f robot_position);
    void get_occupied_points(pcl::PointCloud<pcl::PointXYZI> &out_points) const;
    void get_local_scale(Eigen::Vector3f& position, Eigen::Vector3f& map_min, Eigen::Vector3f& map_max) const;

private:  
    class ROGMapImpl;
    std::unique_ptr<ROGMapImpl> impl_;
};

}

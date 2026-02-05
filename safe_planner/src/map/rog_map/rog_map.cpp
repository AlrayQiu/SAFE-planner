#include "safe_planner/map/rog_map.hpp"
#include "inflation_map.hpp"
#include "sliding_map.hpp"
#include "occupancy_grid_map.hpp"

#include "utils/math.hpp"
#include "utils/raycaster.hpp"

#include "pcl/point_types.h"
#include "pcl/point_cloud.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <cmath>
#include <memory>
#include <queue>

namespace safe_planner::map {
using namespace rog_map;
using namespace occupancy_grid_map;

static const CenterPosition center_position = utils::math::cell::CenterPosition::center_in_center;

class ROGMapMemory final : public rog_map::ICellMemoryReseter
{
public:
    ROGMapMemory(const ROGMapMemory&) = delete;
    ROGMapMemory(ROGMapMemory&&) = delete;
    ROGMapMemory() = delete;
    inline ~ROGMapMemory() = default;
    inline ROGMapMemory(ProblisticMap& ogm,InflationMap<center_position>& inf_map) : ogm_(ogm),inf_map_(inf_map){};
    inline void reset_local_map() {
        ogm_.reset_local_map();
        inf_map_.reset_local_map();
    };
    inline void reset_cell(const int& id) {
        ogm_.reset_cell(id);
        inf_map_.reset_cell(id);
    };

private:
    ProblisticMap&  ogm_;
    InflationMap<center_position>& inf_map_;
};
    
class ROGMap::ROGMapImpl{ 
public:
    ROGMapImpl(
        const Vector3i &half_map_size_i,
        const double &resolution,
        const bool &map_sliding_en,
        const double &sliding_thresh,
        const Vector3f &fix_map_origin,
        Config &&map_config
    ) 
    : config_(std::move(map_config))
    , ogm_(half_map_size_i, config_.prob_config)
    , inf_map_(half_map_size_i,config_.inflation_radius)
    , memory_(ogm_, inf_map_)
    , sliding_map_(
        half_map_size_i,
        resolution,
        map_sliding_en,
        sliding_thresh,
        fix_map_origin,
        memory_)
    , ray_caster_(resolution)
    {
        raycast_data_.hit_count_buffer.resize(static_cast<size_t>(sliding_map_.sliding_config_.map_vox_num), 0);
        raycast_data_.operation_cnt.resize(static_cast<size_t>(sliding_map_.sliding_config_.map_vox_num), 0);
        raycast_data_.raycast_box_range_min = sliding_map_.local_map_bound_min_d_;
        raycast_data_.raycast_box_range_max = sliding_map_.local_map_bound_max_d_;
    }
    
    /*
    * Exam:
    * ```
    * pcl::PointCloud<pcl::PointXYZ> occupied_points;
    * problistic_map.get_occupied_points(occupied_points);
    * ```
    */
    void get_occupied_points(pcl::PointCloud<pcl::PointXYZI>& out_points) const{
        if(!out_points.empty())
            throw std::runtime_error("Output point cloud for problity map must be empty");   
        auto expt_size = static_cast<size_t>(sqrt(sliding_map_.sliding_config_.half_map_size_i.prod()));
        out_points.reserve(expt_size);
        
        Vector3f p{};
        for(int i = sliding_map_.local_map_bound_min_i_.x(); i < sliding_map_.local_map_bound_max_i_.x(); i++)
        for(int j = sliding_map_.local_map_bound_min_i_.y(); j < sliding_map_.local_map_bound_max_i_.y(); j++)
        for(int k = sliding_map_.local_map_bound_min_i_.z(); k < sliding_map_.local_map_bound_max_i_.z(); k++)
        {
            Vector3i global_index{i,j,k};
            auto index = sliding_map_.get_hash_id_from_global_index(global_index);
            auto val = ogm_.occupancy_buffer_[index];
            if(!ogm_.is_occupied(val)) continue;
            
            sliding_map_.global_index_to_pos(global_index, p);                    
            out_points.push_back({p.x(),p.y(),p.z(),val});
        }
    }


    inline void update(const pcl::PointCloud<pcl::PointXYZINormal>& points, const Vector3f &robot_pos){
        
        if(robot_pos.z() > config_.virtual_ceil_height)     return;
        if(robot_pos.z() < config_.virtual_ground_height)   return;

        
        if (raycast_data_.processed_batch_num == 0 && 
            ((robot_pos - sliding_map_.local_map_origin_d_).norm() > sliding_map_.sliding_config_.sliding_thresh))
            slide_all_map(robot_pos);
        
        
        
        update_local_box(robot_pos);
        raycast(points, robot_pos);

        raycast_data_.processed_batch_num++;
        if(raycast_data_.processed_batch_num >= config_.update_pcd_batch_size){
            raycast_data_.processed_batch_num = 0;
            update_map_from_cache();
        }
    }

    void get_local_scale(Eigen::Vector3f& position, Eigen::Vector3f& map_min, Eigen::Vector3f& map_max) const{
        position = robot_position_d_;
        map_min = sliding_map_.local_map_bound_min_d_;
        map_max = sliding_map_.local_map_bound_max_d_;
    }
    

    ~ROGMapImpl() = default;
private:
    inline void update_local_box(const Vector3f &robot_pos){
        robot_position_d_ = robot_pos;
        std::lock_guard<std::mutex> lck(raycast_data_.range_lock_);

        Vector3i robot_pos_i;
        sliding_map_.pos_to_global_index(robot_pos, robot_pos_i);
        Vector3i local_updatebox_min_i, local_updatebox_max_i;

        local_updatebox_max_i = robot_pos_i + config_.half_local_update_box_i;
        local_updatebox_min_i = robot_pos_i - config_.half_local_update_box_i;

        sliding_map_.global_index_to_pos(local_updatebox_min_i, raycast_data_.raycast_box_range_min);
        sliding_map_.global_index_to_pos(local_updatebox_max_i, raycast_data_.raycast_box_range_max);

        // the local update box must inside the local map
        raycast_data_.raycast_box_range_max = raycast_data_.raycast_box_range_max.cwiseMin(
            sliding_map_.local_map_bound_max_d_);
        raycast_data_.raycast_box_range_min = raycast_data_.raycast_box_range_min.cwiseMax(
            sliding_map_.local_map_bound_min_d_);
    }

    inline void slide_all_map(const Vector3f &robot_pos){
        sliding_map_.map_sliding(robot_pos);
    };

    inline void check_state(const float& l_odd_log, GridType& type) noexcept{        
        if (ogm_.is_occupied(l_odd_log)) {
            type = GridType::OCCUPIED;
        }
        else if (ogm_.is_free(l_odd_log)) {
            type = GridType::KNOWN_FREE;
        }
        else {
            type = GridType::UNKNOWN;
        }
    }

    inline void update_hit(const Vector3i& global_index, const int& hash_id, const int& hit_num){
        
        float& ret = ogm_.occupancy_buffer_[hash_id];
        GridType from_type = GridType::UNDEFINED;
        check_state(ret, from_type);

        const auto& cfg = config_.prob_config;
        ret += cfg.l_hit * hit_num;
        if (ret > cfg.l_max) {
            ret = cfg.l_max;
        }

        GridType to_type;
        check_state(ret, to_type);

        if (from_type == to_type) return;

        for(const auto& p : inf_map_.get_local_update_buffer())
            inf_map_.update_cell(
                sliding_map_.get_hash_id_from_global_index(p + global_index),
                from_type,
                to_type);
    } 

    inline void update_miss(const Vector3i& global_index, const int& hash_id, const int& hit_num){
        float& ret = ogm_.occupancy_buffer_[hash_id];
        GridType from_type;
        check_state(ret, from_type);
        
        const auto& cfg = config_.prob_config;
        ret += cfg.l_miss * hit_num;
        if (ret < cfg.l_min) {
            ret = cfg.l_min;
        }

        GridType to_type;
        check_state(ret, to_type);
        // Catch the jump edge
        if (from_type == to_type) return; 
        for(const auto& p : inf_map_.get_local_update_buffer())
            inf_map_.update_cell(
                sliding_map_.get_hash_id_from_global_index(p + global_index),
                from_type,
                to_type);
            
    }

    inline void raycast(const pcl::PointCloud<pcl::PointXYZINormal>& points, const Vector3f &robot_pos){

        Vector3f box_min, box_max;
        {
            std::lock_guard<std::mutex> lock(raycast_data_.range_lock_);
            box_max = raycast_data_.raycast_box_range_max;
            box_min = raycast_data_.raycast_box_range_min;
        }
        
        std::vector<Vector3f> raycasting_pcd{};
        raycasting_pcd.reserve(points.size());

        int sample_cnt{0};
        for(const auto& point : points){
            // intensity filter
            if(config_.intensity_thresh > 0 && point.intensity < config_.intensity_thresh)
                continue;

            
            // downsample
            if (sample_cnt++ % config_.downsample_num)
                continue;

            Vector3f p(point.x, point.y, point.z);
            Vector3i point_index_global;            
            bool     hitted{true};
            
            if(p.z() > config_.virtual_ceil_height){
                hitted = false;
                const double dz = p.z() - robot_pos.z();
                const double pc = config_.virtual_ceil_height - robot_pos.z();
                p = robot_pos + (p - robot_pos).normalized() * pc / dz;
            }
            else if (p.z() < config_.virtual_ground_height) {
                hitted = false;
                const double dz = p.z() - robot_pos.z();
                const double pg = config_.virtual_ground_height - robot_pos.z();
                p = robot_pos + (p - robot_pos).normalized() * pg / dz;
            }

            // near uva
            if ((p - robot_position_d_).norm() < config_.uva_nearby_zone) hitted = false;

            // Bounding box
            if (((p - box_min).minCoeff() < 0) || ((p - box_max).maxCoeff() > 0)) {
                p = line_box_intersect_point(p,
                                        robot_pos,
                                        box_min,
                                        box_max);
                hitted = false;
            }

            raycasting_pcd.push_back(p);

            if(!hitted) continue;
            sliding_map_.pos_to_global_index(p, point_index_global);
            insert_update_candidate(point_index_global, true);
        }

        for(const auto& p : raycasting_pcd){
            if(!ray_caster_.set_input(robot_position_d_, p))
                continue;
            Eigen::Vector3i ray_point;
            while(ray_caster_.next(ray_point) && sliding_map_.inside_local_map(ray_point)){ 
                insert_update_candidate(ray_point, false);
            }
        }
    }

    void insert_update_candidate(const Vector3i& id_g, bool is_hit) {
        const auto& hash_id = sliding_map_.get_hash_id_from_global_index(id_g);
        raycast_data_.operation_cnt[hash_id]++;
        if (raycast_data_.operation_cnt[hash_id] == 1) 
            raycast_data_.changed_cell_global_id_list.push(id_g);
        if (is_hit)
            raycast_data_.hit_count_buffer[hash_id]++;           
    } 

    void update_map_from_cache(){
        while (!raycast_data_.changed_cell_global_id_list.empty()) {
            const auto& global_id = raycast_data_.changed_cell_global_id_list.front();
            const auto& hash_id = sliding_map_.get_hash_id_from_global_index(global_id);
            const int hit_count = raycast_data_.hit_count_buffer[hash_id];
            const int total_count = raycast_data_.operation_cnt[hash_id];

            if(hit_count > 0)
                update_hit(global_id ,hash_id, hit_count);
            else
                update_miss(global_id,hash_id, total_count);
            raycast_data_.hit_count_buffer[hash_id] = 0;
            raycast_data_.operation_cnt[hash_id] = 0;
            raycast_data_.changed_cell_global_id_list.pop();
        }
    }

    struct RaycastData{
        int processed_batch_num{0};
        std::vector<int> hit_count_buffer;
        std::vector<int> operation_cnt;
        std::queue<Vector3i>  changed_cell_global_id_list;
        Vector3f raycast_box_range_min{};
        Vector3f raycast_box_range_max{};
        std::mutex range_lock_;
    } raycast_data_;

    Vector3f robot_position_d_;
    
    Config  config_;
    
    ProblisticMap                   ogm_;
    InflationMap<center_position>   inf_map_;
    
    ROGMapMemory                                memory_;
    SlidingMap<ROGMapMemory,center_position>    sliding_map_;

    utils::raycaster::Raycaster<center_position>    ray_caster_;
};

ROGMap::ROGMap(
        const Vector3i &half_map_size_i,
        const double &resolution,
        const bool &map_sliding_en,
        const double &sliding_thresh,
        const Vector3f &fix_map_origin,
        Config &&map_config)
        :impl_(std::make_unique<ROGMapImpl>(    
            half_map_size_i,
            resolution,
            map_sliding_en,
            sliding_thresh,
            fix_map_origin,
            std::move(map_config)
        )){

}

ROGMap::~ROGMap() = default;

void ROGMap::update(const pcl::PointCloud<pcl::PointXYZINormal>& pcd, const Eigen::Vector3f robot_position)
{
    impl_->update(pcd, robot_position);
}

void ROGMap::get_occupied_points(pcl::PointCloud<pcl::PointXYZI> &out_points) const{
    return impl_->get_occupied_points(out_points);
}

void ROGMap::get_local_scale(Eigen::Vector3f& position, Eigen::Vector3f& map_min, Eigen::Vector3f& map_max) const{
    impl_->get_local_scale(position, map_min, map_max);
}
}
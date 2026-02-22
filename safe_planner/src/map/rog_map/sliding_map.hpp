/*
 * sliding_map.hpp
 *
 *  Created on: June 15, 2022
 * =============================================
 * Sliding Map implementation
 * Refer to "ROG-Map: A Robust and Efficient 3D Mapping Framework for
 *    Autonomous Robotics" (IROS 2022)
 * Paper link: https://doi.org/10.1109/IROS58592.2024.10802303
 */
#pragma once

#include <Eigen/Eigen>
#include <iostream>

#include "utils/math.hpp"

namespace safe_planner::map::rog_map {

class ICellMemoryReseter{
public:
    virtual void reset_local_map() = 0;
    virtual void reset_cell(const int id) = 0;
};

using namespace Eigen;
using namespace utils::math;
using namespace utils::math::cell;
template<class T, CenterPosition center = CenterPosition::center_in_cornor> 
requires std::derived_from<T, ICellMemoryReseter>
class SlidingMap{
public:
    inline SlidingMap(
        const Vector3i &half_map_size_i,
        const double &resolution,
        const bool map_sliding_en,
        const double &sliding_thresh,
        const Vector3f &fix_map_origin,
        T& cell_memory_reseter 
    ) : reseter_(cell_memory_reseter)
    {
        sliding_config_.resolution = resolution;
        sliding_config_.resolution_inv = 1.0 / resolution;
        sliding_config_.map_sliding_en = map_sliding_en;
        sliding_config_.sliding_thresh = sliding_thresh;
        sliding_config_.fix_map_origin = fix_map_origin;
        sliding_config_.half_map_size_i = half_map_size_i;
        sliding_config_.map_size_i = 2 * sliding_config_.half_map_size_i + Vector3i::Constant(1);
        sliding_config_.map_vox_num = sliding_config_.map_size_i.prod();
        if (!map_sliding_en) {
            local_map_origin_d_ = fix_map_origin;
            pos_to_global_index(local_map_origin_d_, local_map_origin_i_);
        }
        map_sliding({0,0,0});
    }
    
    inline SlidingMap() = delete;
    inline SlidingMap(const SlidingMap&) = delete;
    inline SlidingMap& operator=(const SlidingMap&) = delete;
    inline ~SlidingMap() = default;

                
    inline void print_map_info()
    {
        std::cout<<"half_map_size_i: "<<sliding_config_.half_map_size_i.transpose()<<std::endl;
        std::cout<<"resolution: "<<sliding_config_.resolution<<std::endl;
        std::cout<<"map_sliding_en: "<<sliding_config_.map_sliding_en<<std::endl;
        std::cout<<"sliding_thresh: "<<sliding_config_.sliding_thresh<<std::endl;
        std::cout<<"fix_map_origin: "<<sliding_config_.fix_map_origin.transpose()<<std::endl;
        std::cout<<"\tresolution: " << sliding_config_.resolution<< std::endl;
        std::cout<<"\tmap_sliding_en: " << sliding_config_.map_sliding_en<< std::endl;
        std::cout<<"\tlocal_map_size_i: " << sliding_config_.map_size_i.transpose()<< std::endl;
        std::cout<<"\tlocal_map_size_d: " << sliding_config_.map_size_i.template cast<double>().transpose() * sliding_config_.resolution<< std::endl;
    }


    inline void map_sliding(const Vector3f &odom){
        /// Compute the shifting index
        Vector3i new_origin_i;
        pos_to_global_index(odom, new_origin_i);
        Vector3f new_origin_d = new_origin_i.cast<float>() * sliding_config_.resolution;
        /// Compute the delta shift
        Vector3i shift_num = new_origin_i - local_map_origin_i_;
        for (long unsigned int i = 0; i < 3; i++) {
            if (std::fabs(shift_num[i]) > sliding_config_.map_size_i[i]) {
                // Clear all map
                reseter_.reset_local_map();
                update_local_map_origin_and_bound(new_origin_d, new_origin_i);
                return;
            }
        }
        static auto normalize = [](int x, int a, int b) -> int {
            int range = b - a + 1;
            int y = (x - a) % range;
            return (y < 0 ? y + range : y) + a;
        };

        /// Clear the memory out of the map size
        for (int i = 0; i < 3; i++) {
            if (shift_num[i] == 0) {
                continue;
            }
            int min_id_g = -sliding_config_.half_map_size_i(i) + local_map_origin_i_(i);
            int min_id_l = min_id_g % sliding_config_.map_size_i(i);
            std::vector<int> clear_id;
            if (shift_num(i) > 0) {
                /// forward shift, the min id should be cut
                for (int k = 0; k < shift_num(i); k++) {
                    int temp_id = min_id_l + k;
                    temp_id = normalize(temp_id, -sliding_config_.half_map_size_i(i), sliding_config_.half_map_size_i(i));
                    clear_id.push_back(temp_id);
                }
            } else {
                /// backward shift, the max should be shifted
                for (int k = -1; k >= shift_num(i); k--) {
                    int temp_id = min_id_l + k;
                    temp_id = normalize(temp_id, -sliding_config_.half_map_size_i(i), sliding_config_.half_map_size_i(i));
                    clear_id.push_back(temp_id);
                }
            }

            if (clear_id.empty()) {
                continue;
            }
            clear_memory_out_of_map(clear_id, i);
        }

        update_local_map_origin_and_bound(new_origin_d, new_origin_i);
    }

    inline bool inside_local_map(const Vector3f &global_pos) const{
        Vector3i global_id;
        pos_to_global_index(global_pos, global_id);
        return inside_local_map(global_id);
    }

    inline bool inside_local_map(const Vector3i &global_index) const{
        if (((global_index - local_map_origin_i_).cwiseAbs() - sliding_config_.half_map_size_i).maxCoeff() > 0) {
            return false;
        }
        return true;
    }

    struct SlidingConfig {
        double      resolution{0.0};
        double      resolution_inv{0.0};
        double      sliding_thresh{0.0};
        bool        map_sliding_en{false};
        Vector3f    fix_map_origin{};
        Vector3i    visualization_range_i{};
        Vector3i    map_size_i{};
        Vector3i    half_map_size_i{};
        int         virtual_ceil_height_id_g{0};
        int         virtual_ground_height_id_g{0};
        int         map_vox_num{0};
    } sliding_config_;

    

    Vector3f local_map_origin_d_, local_map_bound_min_d_, local_map_bound_max_d_;
    Vector3i local_map_origin_i_, local_map_bound_min_i_, local_map_bound_max_i_;

    inline void clear_memory_out_of_map(const std::vector<int> &clear_ids, const int axis){
        std::vector<int> ids{axis, (axis + 1) % 3, (axis + 2) % 3};
        for (const auto &idd: clear_ids) {
            for (int x = -sliding_config_.half_map_size_i(ids[1]); x <= sliding_config_.half_map_size_i(ids[1]); x++) {
                for (int y = -sliding_config_.half_map_size_i(ids[2]); y <= sliding_config_.half_map_size_i(ids[2]); y++) {
                    Vector3i temp_clear_id;
                    temp_clear_id(ids[0]) = idd;
                    temp_clear_id(ids[1]) = x;
                    temp_clear_id(ids[2]) = y;
                    reseter_.reset_cell(get_local_index_hash(temp_clear_id));
                }
            }
        }
    }

    inline int get_local_index_hash(const Vector3i &index_in) const{  
        Vector3i id = index_in + sliding_config_.half_map_size_i;
        return  id(0) * sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2) +
                id(1) * sliding_config_.map_size_i(2) + id(2);
    }

    inline void pos_to_global_index(const Vector3f &pos, Vector3i &index) const{
        if constexpr (center == CenterPosition::center_in_center)
            index = (sliding_config_.resolution_inv * pos + pos.cwiseSign() * 0.5).template cast<int>();

        else if constexpr (center == CenterPosition::center_in_cornor)
            index = (pos.array() * sliding_config_.resolution_inv).floor().template cast<int>();
    }

    inline void pos_to_hash(const double &pos, int &index) const{
        if constexpr (center == CenterPosition::center_in_center)
            index = static_cast<int>((sliding_config_.resolution_inv * pos + sign(pos) * 0.5));
        else if constexpr (center == CenterPosition::center_in_cornor)
            index = std::floor(pos * sliding_config_.resolution_inv);
    }

    void global_index_to_pos(const Vector3i &global_index, Vector3f &pos) const{
        if constexpr (center == CenterPosition::center_in_center)
            pos = global_index.template cast<float>() * sliding_config_.resolution;
        else if constexpr (center == CenterPosition::center_in_cornor) 
            pos = (global_index.template cast<float>() + Vector3f(0.5, 0.5, 0.5)) * sliding_config_.resolution;
    }

    void global_index_to_local_index(const Vector3i &global_index, Vector3i &local_index) const{
        for (int i = 0; i < 3; ++i) {
            /// [eq. (7) in paper] Compute the i_k
            local_index(i) = global_index(i) % sliding_config_.map_size_i(i);
            /// [eq. (8) in paper] Normalize the local index
            if (local_index(i) > sliding_config_.half_map_size_i(i)) {
                local_index(i) -= sliding_config_.map_size_i(i);
            } else if (local_index(i) < -sliding_config_.half_map_size_i(i)) {
                local_index(i) += sliding_config_.map_size_i(i);
            }
        }
    }

    /* Only used in clearMemoryOutOfMap and */
    inline void local_index_to_global_index(const Vector3i &local_index, Vector3i &global_index) const{
        for (int i = 0; i < 3; ++i) {
            int min_id_g = -sliding_config_.half_map_size_i(i) + local_map_origin_i_(i);
            int min_id_l = min_id_g % sliding_config_.map_size_i(i);
            min_id_l -= min_id_l > sliding_config_.half_map_size_i(i) ? sliding_config_.map_size_i(i) : 0;
            min_id_l += min_id_l < -sliding_config_.half_map_size_i(i) ? sliding_config_.map_size_i(i) : 0;
            int cur_dis_to_min_id = local_index(i) - min_id_l;
            cur_dis_to_min_id =
                    (cur_dis_to_min_id) < 0 ? (sliding_config_.map_size_i(i) + cur_dis_to_min_id) : cur_dis_to_min_id;
            int cur_id = cur_dis_to_min_id + min_id_g;
            global_index(i) = cur_id;
        }
    }

    inline void local_index_to_pos(const Vector3i &local_index, Vector3f &pos) const{
        if constexpr (center == CenterPosition::center_in_center) {
            for (int i = 0; i < 3; ++i) {
                int min_id_g = -sliding_config_.half_map_size_i(i) + local_map_origin_i_(i);

                int min_id_l = min_id_g % sliding_config_.map_size_i(i);
                min_id_l -= min_id_l > sliding_config_.half_map_size_i(i) ? sliding_config_.map_size_i(i) : 0;
                min_id_l += min_id_l < -sliding_config_.half_map_size_i(i) ? sliding_config_.map_size_i(i) : 0;

                int cur_dis_to_min_id = local_index(i) - min_id_l;
                cur_dis_to_min_id =
                        (cur_dis_to_min_id) < 0 ? (sliding_config_.map_size_i(i) + cur_dis_to_min_id) : cur_dis_to_min_id;
                int cur_id = cur_dis_to_min_id + min_id_g;
                pos(i) = cur_id * sliding_config_.resolution;
            }
        }

        if constexpr (center == CenterPosition::center_in_cornor) {
            for (int i = 0; i < 3; ++i) {
                int min_id_g = -sliding_config_.half_map_size_i(i) + local_map_origin_i_(i);
    
                int min_id_l = min_id_g % sliding_config_.map_size_i(i);
                min_id_l -= min_id_l > sliding_config_.half_map_size_i(i) ? sliding_config_.map_size_i(i) : 0;
                min_id_l += min_id_l < -sliding_config_.half_map_size_i(i) ? sliding_config_.map_size_i(i) : 0;
    
                int cur_dis_to_min_id = local_index(i) - min_id_l;
                cur_dis_to_min_id =
                        (cur_dis_to_min_id) < 0 ? (sliding_config_.map_size_i(i) + cur_dis_to_min_id) : cur_dis_to_min_id;
                int cur_id = cur_dis_to_min_id + min_id_g;
                pos(i) = (cur_id + 0.5) * sliding_config_.resolution;
            }
        }
    }

    inline void hash_id_to_local_index(const int hash_id, Vector3i &local_index) const{                        
        local_index(0) = hash_id / (sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2));
        local_index(1) = (hash_id - local_index(0) * sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2)) / sliding_config_.map_size_i(2);
        local_index(2) = hash_id - local_index(0) * sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2) - local_index(1) * sliding_config_.map_size_i(2);
        local_index -= sliding_config_.half_map_size_i;
    }

    inline void hash_id_to_pos(const int hash_id, Vector3f &pos) const{
        Vector3i id;
        hash_id_to_local_index(hash_id, id);
        local_index_to_pos(id, pos);
    }

    inline void hash_id_to_global_index(const int hash_id, Vector3i &global_index) const{        
        Vector3i id;
        id(0) = hash_id / (sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2));
        id(1) = (hash_id - id(0) * sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2)) / sliding_config_.map_size_i(2);
        id(2) = hash_id - id(0) * sliding_config_.map_size_i(1) * sliding_config_.map_size_i(2) - id(1) * sliding_config_.map_size_i(2);
        id -= sliding_config_.half_map_size_i;
        local_index_to_global_index(id, global_index);
    }

    inline int get_hash_index_from_pos(const Vector3f &pos) const{
        Vector3i id_g, id_l;
        pos_to_global_index(pos, id_g);
        global_index_to_local_index(id_g, id_l);
        return get_local_index_hash(id_l);
    }

    inline int get_hash_id_from_global_index(const Vector3i &global_index) const{
        Vector3i local_index;
        global_index_to_local_index(global_index, local_index);
        return get_local_index_hash(local_index);
    }

    inline void update_local_map_origin_and_bound(const Vector3f &new_origin_d, const Vector3i &new_origin_i){
                                                // update local map origin and local map bound
        local_map_origin_i_ = new_origin_i;
        local_map_origin_d_ = new_origin_d;
    
        local_map_bound_max_i_ = local_map_origin_i_ + sliding_config_.half_map_size_i;
        local_map_bound_min_i_ = local_map_origin_i_ - sliding_config_.half_map_size_i;

    
        // the double map bound only consider the closed cell center
        global_index_to_pos(local_map_bound_min_i_, local_map_bound_min_d_);
        global_index_to_pos(local_map_bound_max_i_, local_map_bound_max_d_);
    }


    private:
        T& reseter_;
};

typedef SlidingMap<ICellMemoryReseter> SlidingMapCommon;
}
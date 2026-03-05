/*
 * raycaster.hpp
 *
 *  Created on: Jun 8, 2021
 *
 * 算法使用bresham算法实现3D空间中的射线投射
 * 实现参考 https://github.com/hku-mars/ROG-Map
*/

#pragma once

#include <Eigen/Eigen>
#include <cmath>
#include <stdexcept>

#include "math.hpp"
#include "safe_planner/map/imap.hpp"

namespace safe_planner::map::utils::raycaster {
using namespace cell;
template <CenterPosition center = CenterPosition::center_in_cornor>
class Raycaster {
public:
    inline Raycaster(const double& resolution_in)
        : resolution_(resolution_in) {
            if(resolution_in < 0)
                throw std::runtime_error("resolution should be positive");
        };
    Raycaster() = delete;
    ~Raycaster() = default;

    bool set_input(const Eigen::Vector3f &start, const Eigen::Vector3f &end) noexcept {
        start_x_d_ = start.x();
        start_y_d_ = start.y();
        start_z_d_ = start.z();

        end_x_d_ = end.x();
        end_y_d_ = end.y();
        end_z_d_ = end.z();

        // Calculate expand dir and index
        pos2index(start_x_d_, start_x_i_);
        pos2index(start_y_d_, start_y_i_);
        pos2index(start_z_d_, start_z_i_);

        pos2index(end_x_d_, end_x_i_);
        pos2index(end_y_d_, end_y_i_);
        pos2index(end_z_d_, end_z_i_);

        int delta_X = end_x_i_ - start_x_i_;
        int delta_Y = end_y_i_ - start_y_i_;
        int delta_Z = end_z_i_ - start_z_i_;

        cur_ray_pt_id_x_ = start_x_i_;
        cur_ray_pt_id_y_ = start_y_i_;
        cur_ray_pt_id_z_ = start_z_i_;

        expand_dir_x_ = static_cast<int>(math::sign(delta_X));
        expand_dir_y_ = static_cast<int>(math::sign(delta_Y));
        expand_dir_z_ = static_cast<int>(math::sign(delta_Z));
        if (expand_dir_x_ == 0 && expand_dir_y_ == 0 && expand_dir_z_ == 0) {
            return false;
        }

        double dis_x_over_t, dis_y_over_t, dis_z_over_t, t_max;
        dis_x_over_t = std::abs(end_x_d_ - start_x_d_);
        dis_y_over_t = std::abs(end_y_d_ - start_y_d_);
        dis_z_over_t = std::abs(end_z_d_ - start_z_d_);

        t_max = sqrt(dis_x_over_t * dis_x_over_t + dis_y_over_t * dis_y_over_t + dis_z_over_t * dis_z_over_t);

        dis_x_over_t /= t_max;
        dis_y_over_t /= t_max;
        dis_z_over_t /= t_max;

        t_when_step_x_ = expand_dir_x_ == 0 ? std::numeric_limits<double>::max() :
                            std::abs(resolution_ / dis_x_over_t);
        t_when_step_y_ = expand_dir_y_ == 0 ? std::numeric_limits<double>::max() :
                            std::abs(resolution_ / dis_y_over_t);
        t_when_step_z_ = expand_dir_z_ == 0 ? std::numeric_limits<double>::max() :
                            std::abs(resolution_ / dis_z_over_t);

        double start_grid_center_x, start_grid_center_y, start_grid_center_z;
        index2pos(start_x_i_, start_grid_center_x);
        index2pos(start_y_i_, start_grid_center_y);
        index2pos(start_z_i_, start_grid_center_z);
        double next_bound_x, next_bound_y, next_bound_z;
        next_bound_x = start_grid_center_x + expand_dir_x_ * resolution_ * 0.5;
        next_bound_y = start_grid_center_y + expand_dir_y_ * resolution_ * 0.5;
        next_bound_z = start_grid_center_z + expand_dir_z_ * resolution_ * 0.5;


        t_to_bound_x_ = expand_dir_x_ == 0 ? std::numeric_limits<double>::max() :
                        std::fabs(next_bound_x - start_x_d_) / dis_x_over_t;
        t_to_bound_y_ = expand_dir_y_ == 0 ? std::numeric_limits<double>::max() :
                        std::fabs(next_bound_y - start_y_d_) / dis_y_over_t;
        t_to_bound_z_ = expand_dir_z_ == 0 ? std::numeric_limits<double>::max() :
                        std::fabs(next_bound_z - start_z_d_) / dis_z_over_t;

        step_num_ = 0;
        first_point = true;
        return true;
    }

    inline constexpr void set_resolution(const double &resolution_in) noexcept {
        if(resolution_in < 0)
            throw std::runtime_error("resolution should be positive");
        resolution_ = resolution_in;
    }

    inline bool next(Eigen::Vector3i &ray_pt) {
        step_num_++;
        ray_pt.x() = cur_ray_pt_id_x_;
        ray_pt.y() = cur_ray_pt_id_y_;
        ray_pt.z() = cur_ray_pt_id_z_;
        
        if (cur_ray_pt_id_x_ == end_x_i_ && cur_ray_pt_id_y_ == end_y_i_ && cur_ray_pt_id_z_ == end_z_i_) {
            return false;
        }

        if (t_to_bound_x_ < t_to_bound_y_) {
            if (t_to_bound_x_ < t_to_bound_z_) {
                cur_ray_pt_id_x_ += expand_dir_x_;
                t_to_bound_x_ += t_when_step_x_;
            } else {
                cur_ray_pt_id_z_ += expand_dir_z_;
                t_to_bound_z_ += t_when_step_z_;
            }
        } else {
            if (t_to_bound_y_ < t_to_bound_z_) {
                cur_ray_pt_id_y_ += expand_dir_y_;
                t_to_bound_y_ += t_when_step_y_;
            } else {
                cur_ray_pt_id_z_ += expand_dir_z_;
                t_to_bound_z_ += t_when_step_z_;

            }
        }
        return true;
    }


private:
    inline constexpr void pos2index(const double d, int &id) const noexcept {
        if constexpr (center == CenterPosition::center_in_cornor)
            id = std::floor((d / resolution_));
        else
            id = static_cast<int>((d / resolution_ + math::sign(d) * 0.5));
    }

    inline constexpr void index2pos(const int &id, double &d) const noexcept { 
        if constexpr (center == CenterPosition::center_in_cornor)
            d = (id + 0.5) * resolution_;
        else if constexpr (center == CenterPosition::center_in_center)
            d = id * resolution_;
    }

    double resolution_{-1};
    bool first_point{true};

    double start_x_d_, start_y_d_, start_z_d_;
    double end_x_d_, end_y_d_, end_z_d_;
    double t_to_bound_x_, t_to_bound_y_, t_to_bound_z_;
    int expand_dir_x_, expand_dir_y_, expand_dir_z_;
    int end_x_i_, end_y_i_, end_z_i_;
    int start_x_i_, start_y_i_, start_z_i_;
    int cur_ray_pt_id_x_, cur_ray_pt_id_y_, cur_ray_pt_id_z_;
    double t_when_step_x_, t_when_step_y_, t_when_step_z_;
    int step_num_{0};
};
}
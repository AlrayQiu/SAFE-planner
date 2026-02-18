
#pragma once

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <random>
#include <vector>
#include <algorithm>

#include <Eigen/Eigen>

#include "safe_planner/map/imap.hpp"
#include "../utils/ikd_tree/ikd_Tree.h"


namespace safe_planner::planner {

namespace rrt_star {
enum class SearchAs : bool{
    free = false,
    occupied = true
};
constexpr struct SearchConfig{
    const float goal_sample_rate = 0.1;
    const float gamma_rrg = 15;
    const float ita = 1;
    const int   n = 1000;
    const float d = 3.0f;
    const float inv_d = 1.0f / 3.0f;
    const SearchAs unknown = SearchAs::free;
    const SearchAs out_of_bound = SearchAs::free;
} config;
struct RRTNode
{
    float x,y,z;
    int id;
    Eigen::Vector3f position() const{return {x,y,z};}

    inline RRTNode() = default;
    inline RRTNode(float x,float y,float z, int id):
    x(x),
    y(y),
    z(z),
    id(id){};

    inline RRTNode(const RRTNode& node) = default;
    inline RRTNode& operator=(const RRTNode& other) = default;
};



typedef KD_TREE<rrt_star::RRTNode> IKdTree;

static inline float rand_float(float min, float max) {
    static thread_local std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<float> dist(min, max);
    return dist(gen);
}

static inline void sample_free(
    const Eigen::Vector3f& box_min, 
    const Eigen::Vector3f& box_max, 
    rrt_star::RRTNode& point)
{
point.x = rand_float(box_min.x(),box_max.x());
point.y = rand_float(box_min.y(),box_max.y());
point.z = rand_float(box_min.z(),box_max.z());
}

static inline void steer(
    const Eigen::Vector3f& x_near, 
    const Eigen::Vector3f& x_sample, 
    const Eigen::Vector3f& goal,
    rrt_star::RRTNode& point){
    auto rand = rand_float(0, 1);
    Eigen::Vector3f dir;
    if(rand > rrt_star::config.goal_sample_rate){
        dir = x_sample - x_near;
    }
    else{
        dir = goal - x_near;
    }
    auto len = dir.norm();
    while(len == 0)
    {
        dir = Eigen::Vector3f::Random();
        len = dir.norm();
    }
    auto step = std::min(rrt_star::config.ita, len) / len;
    auto pos =  dir * step + x_near;
    point.x = pos.x();
    point.y = pos.y();
    point.z = pos.z();
}

static inline float radius(const float card_v){
    return std::max(rrt_star::config.gamma_rrg * std::pow(std::log(card_v)/card_v, rrt_star::config.inv_d), rrt_star::config.ita);
}

static inline float cost(const rrt_star::RRTNode& x, const rrt_star::RRTNode& y)
{
    return (y.position() - x.position()).norm();
}

}

template <class T>
requires std::derived_from<T, IMap>
static inline void rrt_star_search(
    const Eigen::Vector3f& start,
    const Eigen::Vector3f& goal,
    const T& map,
    std::vector<Eigen::Vector3f>&   path_out) 
{  
    path_out.clear();
    Eigen::Vector3f box_min,box_max;
    map.get_map_bound_d(box_min, box_max);
    rrt_star::RRTNode node{start.x(), start.y(), start.z(), 0};

    static std::vector<int>                parent{};
    static std::vector<float>              costs{};
    static std::vector<rrt_star::RRTNode>  nodes{};

    parent.reserve(rrt_star::config.n);
    costs.reserve(rrt_star::config.n);
    nodes.reserve(rrt_star::config.n);
    parent.clear();
    costs.clear();
    nodes.clear();
    
    std::unique_ptr<rrt_star::IKdTree> ikd_tree = std::make_unique<rrt_star::IKdTree>();
    ikd_tree->InitializeKDTree();
    ikd_tree->Build({node});
    parent.emplace_back(0);
    costs.emplace_back(0);
    nodes.emplace_back(node);
    int iter = 0;
    float c_min;
    rrt_star::RRTNode x_sample,x_nearest,x_new;
    int x_min_id;   
    vector<rrt_star::RRTNode,Eigen::aligned_allocator<rrt_star::RRTNode>> X_near{};
    vector<float> dis{};
    while(++iter <= rrt_star::config.n){
        sample_free(box_min, box_max, x_sample);
        ikd_tree->Nearest_Search(x_sample, x_nearest);
        const auto p_x_nearest = x_nearest.position();
        const auto p_x_sample = x_sample.position();
        steer(p_x_nearest,p_x_sample,goal,x_new);
        const auto p_x_new = x_new.position();
        if(map.check_line_d(p_x_nearest, p_x_new) == IMap::State::Unsafe) continue;
        ikd_tree->Radius_Search(x_new, rrt_star::radius(static_cast<float>(ikd_tree->size())), X_near);
        x_min_id = x_nearest.id;
        c_min = costs[x_min_id] + cost(x_new, x_nearest);
        X_near.erase(
            std::remove_if(X_near.begin(), X_near.end(), 
            [&](const auto& xn){ return map.check_line_d(xn.position(), p_x_new) == IMap::State::Unsafe; }), 
            X_near.end());
        for(const auto& x_near : X_near){
            auto c = costs[x_near.id] + cost(x_new, x_near);
            if(c_min <= c) continue;
            c_min = c;
            x_min_id = x_near.id;
        }
        parent.emplace_back(x_min_id);
        costs.emplace_back(c_min);
        x_new.id = nodes.size();
        nodes.emplace_back(x_new);
        auto x_new_cost = c_min;
        for(const auto& x_near : X_near){
            auto c = x_new_cost + cost(x_new, x_near);
            if(costs[x_near.id] <= c) continue;
            parent[x_near.id] = x_new.id;
            costs[x_near.id] = c;
        }
        ikd_tree->Add_Point(x_new);
    }
    rrt_star::RRTNode g{goal.x(),goal.y(),goal.z(), static_cast<int>(nodes.size()) };
    nodes.emplace_back(g);
    ikd_tree->Radius_Search(g, rrt_star::radius(static_cast<float>(ikd_tree->size())), X_near);
    c_min = std::numeric_limits<float>::infinity();
    int p_g = 0;
    for(const auto& x_near : X_near){
        if (map.check_line_d(x_near.position(), goal) == IMap::State::Unsafe) continue;
        auto c = costs[x_near.id] + cost(g, x_near);
        if(c_min < c) continue;
        c_min = c;
        p_g = x_near.id;
    }
    rrt_star::RRTNode curr = g;
    parent.emplace_back(p_g);
    while(curr.id != 0){
        path_out.emplace_back(curr.x,curr.y,curr.z);
        curr = nodes[parent[curr.id]];
    }
    path_out.emplace_back(start);
}
}
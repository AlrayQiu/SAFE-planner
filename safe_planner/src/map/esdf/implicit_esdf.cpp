#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/map/imap.hpp"
#include "safe_planner/map/rog_map.hpp"
#include "../utils/raycaster.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <array>
#include <limits>

#define CLASS(x)  template<class TMap> \
               requires std::derived_from<TMap, safe_planner::IMap> \
               x safe_planner::esdf::ImplicitESDF<TMap>
               
CLASS(class)::Impl{
public:
Impl(const TMap& map, const float resolution)
: map_(map)
, resolution(resolution)
, ray_caster_(resolution)
{ };

std::tuple<float,Eigen::Vector3f> get_esdf(const Eigen::Vector3f& pos, const float , const Eigen::Vector3f&){
    const int r = N / 2;
    Eigen::Vector3i index,g;
    Eigen::Vector3f p;
    map_.pos_to_index(pos, index);
    float min = std::numeric_limits<float>::max();
    for(int i = 0;i < N;++i)
    for(int j = 0;j < N;++j)
    for(int k = 0;k < N;++k)
    {
        const Eigen::Vector3i err{i - r,j - r,k - r};
        const auto n = err + index;
        if(map_.check_point_i(n) != IMap::State::Safe) continue;
        double d = distance_matrix[i][j][k];
        if(min <= d) continue;
        min = d;
        g = n;
    }
        ;
    if(min == std::numeric_limits<float>::max())
        return{0,Eigen::Vector3f::Zero()};
    map_.index_to_pos(g,p);
    p = pos - p;
    if(p.norm() < resolution)
        return{0,Eigen::Vector3f::Zero()};
    p = (p - p.normalized() * resolution);
    return {p.norm(),p.normalized()};
}
private:

static inline const int N = 11;

// 定义一个 N×N 的矩阵，存储到中心的距离
static inline const std::array<std::array<std::array<float, N>, N>, N> 
distance_matrix = [] {
    std::array<std::array<std::array<float, N>, N>, N> mat{};
    const int cx = N / 2;
    const int cy = N / 2;
    const int cz = N / 2;
    for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
    for (int k = 0; k < N; ++k) {
        mat[i][j][k] = std::sqrt((i - cx) * (i - cx) + (j - cy) * (j - cy) + (k - cz) * (k - cz));
    }
    }
    }
    return mat;
}();

const TMap& map_;
const float resolution;
map::utils::raycaster::Raycaster<center_position> ray_caster_;
};
CLASS()::ImplicitESDF(const TMap& map, const float resolution){
    impl_ = std::make_unique<Impl>(map, resolution);
}
CLASS()::~ImplicitESDF() = default;
typedef std::tuple<float,Eigen::Vector3f> ReturnType;
CLASS(ReturnType)::get_esdf(const Eigen::Vector3f& pos, const float init_guess_radius, const Eigen::Vector3f& ref_pos) const{
    return impl_->get_esdf(pos, init_guess_radius,ref_pos);
}

template class safe_planner::esdf::ImplicitESDF<safe_planner::map::ROGMap>;
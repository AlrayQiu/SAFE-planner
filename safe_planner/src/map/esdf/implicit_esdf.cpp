#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/map/rog_map.hpp"
#include <numbers>

#define CLASS(x)  template<class TMap> \
               requires std::derived_from<TMap, safe_planner::IMap> \
               x safe_planner::esdf::ImplicitESDF<TMap>
               
CLASS(class)::Impl{
public:
Impl(const TMap& map, const float resolution)
: map_(map)
, resolution(resolution)
{};

std::tuple<float,Eigen::Vector3f> get_esdf(const Eigen::Vector3f& pos, const float init_guess_radius){
    Eigen::Vector3f g{};
    search_for_negative(pos, init_guess_radius, g);
    return {g.norm(), g.normalized()};
}
private:

void search_for_negative(const Eigen::Vector3f& pos,const float init_guess_radius, Eigen::Vector3f& g){
    auto r = init_guess_radius == 0.0f ? resolution : init_guess_radius;
    auto b = 0.0f;
    bool flag = true;
    do{
        sample_top_.clear();
        sample_sphere(pos, r, resolution, sample_top_);
        for(const auto& p : sample_top_){
            if(map_.check_point_d(p) == IMap::State::Safe){
                flag = false;
                break;
            }
        }
        if(flag){
            b = r;
            r *= 2;
        }
    }while(flag);

    do{
        auto ra = (r + b) / 2;
        sample_.clear();
        sample_sphere(pos, ra, resolution, sample_);
        flag = true;
        for(const auto& p : sample_){
            if(map_.check_point_d(p) == IMap::State::Safe){
                flag = false;
                break;
            }
        }
        if(flag) b = ra;
        else    r = ra;
    }while(std::abs(r - b) > resolution * 0.1);

    sample_.clear();
    g.setZero();
    Eigen::Vector3f one_point{}; 
    sample_sphere(pos, r, resolution, sample_);
    for(const auto& p : sample_){
        if(map_.check_point_d(p) == IMap::State::Safe){
            continue;
        }
        g = p - pos;
        return;
    }
}
static void sample_sphere(
    const Eigen::Vector3f& center, 
    const float radius, const float resolution,
    std::vector<Eigen::Vector3f>& samples){
    for(float p = - radius; p <= radius; p += resolution){
        float r = std::sqrt(radius * radius - p * p);
        sample_circle(center, p, r, resolution, samples);
    }
}
static void sample_circle(const Eigen::Vector3f& center,const float p,const float r,const float resolution,std::vector<Eigen::Vector3f>& samples){
    for(float theta = -std::numbers::pi_v<float>; theta < std::numbers::pi_v<float>; theta += resolution / r){
        samples.emplace_back(center.x() + r * std::cos(theta), center.y() + r * std::sin(theta), center.z() + p);
    }
}
std::vector<Eigen::Vector3f> sample_top_;
std::vector<Eigen::Vector3f> sample_bottom_;
std::vector<Eigen::Vector3f> sample_;

const TMap& map_;
const float resolution;
};
CLASS()::ImplicitESDF(const TMap& map, const float resolution){
    impl_ = std::make_unique<Impl>(map, resolution);
}
CLASS()::~ImplicitESDF() = default;
typedef std::tuple<float,Eigen::Vector3f> ReturnType;
CLASS(ReturnType)::get_esdf(const Eigen::Vector3f& pos, const float init_guess_radius) const{
    return impl_->get_esdf(pos, init_guess_radius);
}

template class safe_planner::esdf::ImplicitESDF<safe_planner::map::ROGMap>;
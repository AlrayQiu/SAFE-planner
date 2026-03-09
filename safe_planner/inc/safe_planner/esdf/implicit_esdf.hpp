#pragma once

#include "safe_planner/map/imap.hpp"
#include <memory>

namespace safe_planner::esdf{
template <typename TMap>
requires std::derived_from<TMap, IMap>
class ImplicitESDF{
public:
    ImplicitESDF() = delete;
    ImplicitESDF(const TMap& map, const float resolution);
    ~ImplicitESDF();
    std::tuple<float,Eigen::Vector3f> get_esdf(const Eigen::Vector3f& , const float ,const Eigen::Vector3f& ref_pos = Eigen::Vector3f::Zero()) const;
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
}
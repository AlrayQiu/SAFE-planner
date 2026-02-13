#include "rrt_star.hpp"
#include "safe_planner/map/rog_map.hpp"
#include "safe_planner/planner/front_end.hpp"

#include <concepts>
#include <memory>
#include <format>
#include <stdexcept>
namespace safe_planner::planner {

template <class TMap>
requires std::derived_from<TMap, IMap>
class FrontEnd<TMap>::FrontEndImpl{
public:
    FrontEndImpl(front_end::Config&& config)
    : config_(std::move(config)) {};

    FrontEndImpl() = delete;
    ~FrontEndImpl() = default;

    void search(
        const Eigen::Vector3f& from,
        const Eigen::Vector3f& to,
        const TMap& map,
        std::vector<Eigen::Vector3f>& path)
    {
        if(config_.algo == front_end::Config::Algo::RRT) 
        {
            rrt_star_search<TMap>(from, to, map, path);
        }
        else throw std::runtime_error(std::format("config.algo no define with code {}", static_cast<int>(config_.algo))); 
    }
private:

    const front_end::Config config_;
};

template <class TMap>
requires std::derived_from<TMap, IMap>
FrontEnd<TMap>::FrontEnd(const TMap& map,front_end::Config&& config)
    :map_(map)
    ,impl_(std::make_unique<FrontEndImpl>(std::move(config))){};

template <class TMap>
requires std::derived_from<TMap, IMap>
FrontEnd<TMap>::~FrontEnd() = default;

template <class TMap>
requires std::derived_from<TMap, IMap>
void FrontEnd<TMap>::search(
    const Eigen::Vector3f& from,
    const Eigen::Vector3f& to,
    std::vector<Eigen::Vector3f>& path)
    {
        impl_->search(from,to,map_,path);
    };

template class FrontEnd<map::ROGMap>;
}
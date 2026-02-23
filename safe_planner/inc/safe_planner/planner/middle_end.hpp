#pragma once

#include "safe_planner/map/imap.hpp"
#include "safe_planner/trajectory/trajectroies.hpp"
#include <concepts>
#include <memory>
#include <vector>

namespace safe_planner::planner
{
namespace middle_end{
enum class Strategy{

};
}
template<class TMap, class TTraj>
requires std::derived_from<TMap, IMap> && std::derived_from<TTraj, ITrajectory>
class MiddleEnd{
public:
    MiddleEnd(const TMap&);
    MiddleEnd() = delete;
    ~MiddleEnd();

    void optimize(
        const std::vector<Eigen::Vector3f>& ref_path,
        const Eigen::Matrix<float,2,3>& begin_va,
        const Eigen::Matrix<float,2,3>& end_va,
        TTraj& traj);
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
}
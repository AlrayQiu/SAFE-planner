#pragma once

#include "../../front_end/jps.hpp"
#include "safe_planner/map/imap.hpp"
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>
#include <Eigen/Eigen>

namespace safe_planner::planner::middle_end::optimizer  {
    
class PVPairs{
public:
inline PVPairs()
:sf2(sf * sf)
,sf3(sf * sf2)
{
    pv_.reserve(10);
}
inline void cost(const Eigen::Vector3d& p, Eigen::Vector3d& gradp,double& costp)
{
    gradp.setZero();
    costp = 0;
    for(const auto& [q,v] : pv_){
        const auto d = (p - q).dot(v);
        const auto c = sf - d;
        if(c <= 0) continue;
        const auto c2 = c * c;
        const auto c3 = c * c2;
        const auto c4 = c2 * c2;
        const auto dcdp = -v;
        costp += c4;
        gradp += 4 * c3 * dcdp;
        continue;
        
    }
}
const double sf = 0.2;
const double sf2;
const double sf3;
std::vector<std::tuple<Eigen::Vector3d,Eigen::Vector3d>> pv_;
};
struct Segment{
int begin;
int end;
};

template<class TMap>
class Ego{
public:
inline Ego(size_t size,const TMap& map)
: n(size)
, map_(map) {
    segs_.reserve(n);
    pvs_.resize(n);
    is_new_.resize(n,false);
    flag.resize(n,false);
}

inline void collision_cost(const Eigen::Matrix3Xd& p, Eigen::Matrix3Xd& gradp, double& costp){
    int from = 0,to;
    seging = false;
    segs_.clear();
    
    for(int i = 1;i < n;i++){
        if(map_.check_point_d(p.col(i).cast<float>()) == IMap::State::Unsafe){
            seging = true;
            continue;
        }
        
        if(!seging){
            from = i;
            continue;
        }
        seging = false;
        to = i;
        segs_.emplace_back(from,to);
        from = to;
    }
    for(const auto& seg : segs_){
        bool has_new = false;
        for(int i = seg.begin + 1;i < seg.end - 1;++i){
            is_new_[i] = is_new_obs(i, p.col(i));
            if(!is_new_[i]) continue;
            has_new = true;
        }
        if(!has_new) continue;
        std::vector<Eigen::Vector3f> jps_pathf; 
        std::vector<Eigen::Vector3d> jps_path; 
        jps_search(p.col(seg.begin).cast<float>(), p.col(seg.end).cast<float>(), map_, jps_pathf);
        if(jps_pathf.size() == 2) continue;
        jps_path.reserve(jps_pathf.size());
        for(const auto& p : jps_pathf)
            jps_path.emplace_back(p.cast<double>());
        int got_intersection_id = -1;
        
        for (int i = seg.begin; i <= seg.end; ++i)
            flag[i] = false;

        for(int i = seg.begin + 1;i < seg.end - 1;++i){
            if(!is_new_[i]) continue;

            int jps_id = jps_path.size() / 2, jps_path_size = jps_path.size(),last_jps_id;
            Eigen::Vector3d p_law = (p.col(i + 1) - p.col(i-1)),intersection_point;
            float val = (jps_path[jps_id] - p.col(i)).dot(p_law);
            float last_val = val;
            while(true){
                last_jps_id = jps_id;
                if(val < 0)
                    ++jps_id;
                else
                    --jps_id;
                if(jps_id < 0 || jps_id >= jps_path_size) break;
                val = (jps_path[jps_id] - p.col(i)).dot(p_law);
                if(val * last_val > 0 || (abs(val) == 0 && abs(last_val) == 0)) continue;
                intersection_point =
                    jps_path[jps_id] +
                    ((jps_path[jps_id] - jps_path[last_jps_id]) *
                    (p_law.dot(p.col(i) - jps_path[jps_id]) / p_law.dot(jps_path[jps_id] - jps_path[last_jps_id]))
                    );
                got_intersection_id = i;
                break;
            }
            while(got_intersection_id >= 0)
            {
                double length = (intersection_point - p.col(i)).norm();
                if(length <= 1e-5){
                    got_intersection_id = -1;
                    break;
                }
                flag[i] = true;
                for (double a = length; a >= 0.0; a -= map_.resolution())
                {
                    
                    bool occ = map_.check_point_d(((a / length) * intersection_point + (1 - a / length) * p.col(i)).cast<float>()) == IMap::State::Unsafe;
                    
                    if (!occ && a > map_.resolution()) continue;
                    got_intersection_id = -1;
                    if (!occ) break;
                    got_intersection_id = i;
                    a += map_.resolution();
                    Eigen::Vector3d p_ = ((a / length) * intersection_point + (1 - a / length) * p.col(i));
                    Eigen::Vector3d v_ = (intersection_point - p.col(i)).normalized();
                    pvs_[i].pv_.emplace_back(p_, v_);
                    break;
                }
                break;
            }
        }
        if (got_intersection_id >= 0)
        {
          for (int j = got_intersection_id + 1; j <= seg.end; ++j)
            if (!flag[j]&& is_new_[j])
            {
                pvs_[j].pv_.emplace_back(pvs_[j - 1].pv_.back());
            }

          for (int j = got_intersection_id - 1; j >= seg.begin; --j)
            if (!flag[j] && is_new_[j])
            {
                pvs_[j].pv_.emplace_back(pvs_[j + 1].pv_.back());
            }
        }
    }
    double cost;
    Eigen::Vector3d grad;
    for(int i = 0;i < n;++i){
        pvs_[i].cost(p.col(i), grad, cost);
        gradp.col(i) << grad;
        costp += cost;
    }
}

private:
inline bool is_new_obs(const int i,const Eigen::Vector3d& q){
    for(const auto& [p,v] : pvs_[i].pv_)
    {
        const auto d = (p - q).dot(v);
        if(d >= 0) return false;
    }
    return true;
}
std::vector<Segment> segs_;
std::vector<PVPairs> pvs_;
std::vector<bool> is_new_;
std::vector<bool> flag;
const int n;
const TMap& map_;
bool new_loop = true;
bool last = false;
bool seging = false;
};
}
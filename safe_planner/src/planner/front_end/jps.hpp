/*url: https://github.com/KumarRobotics/jps3d/ */
#pragma once
#include "safe_planner/map/imap.hpp"
#include <Eigen/Eigen>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <unordered_map>
#include <vector>

namespace safe_planner::planner {

namespace jps {
static inline const int jump_size = 10;
struct JpsNode{
JpsNode(const Eigen::Vector3i& p,const int i):
point(p),
dir(0,0,0),
parent_id(-1),
h(0),
g(std::numeric_limits<double>::max()),
open(false),
close(false),
id(i)
{}

Eigen::Vector3i point;
Eigen::Vector3i dir{0,0,0};
long long parent_id = -1;
double h = 0;
double g = std::numeric_limits<double>::max();
bool open = false;
bool close = false;
int id;
};

struct JPS3DNeib {
// for each (dx,dy,dz) these contain:
//    ns: neighbors that are always added
//    f1: forced neighbors to check
//    f2: neighbors to add if f1 is forced
int ns[27][3][26];
int f1[27][3][12];
int f2[27][3][12];
// nsz contains the number of neighbors for the four different types of moves:
// no move (norm 0):        26 neighbors always added
//                          0 forced neighbors to check (never happens)
//                          0 neighbors to add if forced (never happens)
// straight (norm 1):       1 neighbor always added
//                          8 forced neighbors to check
//                          8 neighbors to add if forced
// diagonal (norm sqrt(2)): 3 neighbors always added
//                          8 forced neighbors to check
//                          12 neighbors to add if forced
// diagonal (norm sqrt(3)): 7 neighbors always added
//                          6 forced neighbors to check
//                          12 neighbors to add if forced
static inline constexpr int nsz[4][2] = {{26, 0}, {1, 8}, {3, 12}, {7, 12}};
JPS3DNeib(){
    int id = 0;
    for(int dz = -1; dz <= 1; ++dz) {
        for(int dy = -1; dy <= 1; ++dy) {
        for(int dx = -1; dx <= 1; ++dx) {
            int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
            for(int dev = 0; dev < nsz[norm1][0]; ++dev)
                neib(dx,dy,dz,norm1,dev,
                    ns[id][0][dev], ns[id][1][dev], ns[id][2][dev]);
            for(int dev = 0; dev < nsz[norm1][1]; ++dev)
            {
                force_neib(dx,dy,dz,norm1,dev,
                    f1[id][0][dev],f1[id][1][dev], f1[id][2][dev],
                    f2[id][0][dev],f2[id][1][dev], f2[id][2][dev]);
            }
            id ++;
        }
        }
    }
}
private:
void neib(int dx, int dy, int dz, int norm1, int dev, int& tx, int& ty, int& tz)
{
     switch(norm1)
  {
    case 0:
      switch(dev)
      {
        case 0: tx=1; ty=0; tz=0; return;
        case 1: tx=-1; ty=0; tz=0; return;
        case 2: tx=0; ty=1; tz=0; return;
        case 3: tx=1; ty=1; tz=0; return;
        case 4: tx=-1; ty=1; tz=0; return;
        case 5: tx=0; ty=-1; tz=0; return;
        case 6: tx=1; ty=-1; tz=0; return;
        case 7: tx=-1; ty=-1; tz=0; return;
        case 8: tx=0; ty=0; tz=1; return;
        case 9: tx=1; ty=0; tz=1; return;
        case 10: tx=-1; ty=0; tz=1; return;
        case 11: tx=0; ty=1; tz=1; return;
        case 12: tx=1; ty=1; tz=1; return;
        case 13: tx=-1; ty=1; tz=1; return;
        case 14: tx=0; ty=-1; tz=1; return;
        case 15: tx=1; ty=-1; tz=1; return;
        case 16: tx=-1; ty=-1; tz=1; return;
        case 17: tx=0; ty=0; tz=-1; return;
        case 18: tx=1; ty=0; tz=-1; return;
        case 19: tx=-1; ty=0; tz=-1; return;
        case 20: tx=0; ty=1; tz=-1; return;
        case 21: tx=1; ty=1; tz=-1; return;
        case 22: tx=-1; ty=1; tz=-1; return;
        case 23: tx=0; ty=-1; tz=-1; return;
        case 24: tx=1; ty=-1; tz=-1; return;
        case 25: tx=-1; ty=-1; tz=-1; return;
      }
      return;
    case 1:
      tx = dx; ty = dy; tz = dz; return;
    case 2:
      switch(dev)
      {
        case 0:
          if(dz == 0){
            tx = 0; ty = dy; tz = 0; return;
          }else{
            tx = 0; ty = 0; tz = dz; return;
          }
        case 1:
          if(dx == 0){
            tx = 0; ty = dy; tz = 0; return;
          }else{
            tx = dx; ty = 0; tz = 0; return;
          }
        case 2:
          tx = dx; ty = dy; tz = dz; return;
      }return;
    case 3:
      switch(dev)
      {
        case 0: tx = dx; ty =  0; tz =  0; return;
        case 1: tx =  0; ty = dy; tz =  0; return;
        case 2: tx =  0; ty =  0; tz = dz; return;
        case 3: tx = dx; ty = dy; tz =  0; return;
        case 4: tx = dx; ty =  0; tz = dz; return;
        case 5: tx =  0; ty = dy; tz = dz; return;
        case 6: tx = dx; ty = dy; tz = dz; return;
      }return;
  }
}
void force_neib( int dx, int dy, int dz, int norm1, int dev,
    int& fx, int& fy, int& fz,
    int& nx, int& ny, int& nz){
        switch(norm1)
  {
    case 1:
      switch(dev)
      {
        case 0: fx= 0; fy= 1; fz = 0; break;
        case 1: fx= 0; fy=-1; fz = 0; break;
        case 2: fx= 1; fy= 0; fz = 0; break;
        case 3: fx= 1; fy= 1; fz = 0; break;
        case 4: fx= 1; fy=-1; fz = 0; break;
        case 5: fx=-1; fy= 0; fz = 0; break;
        case 6: fx=-1; fy= 1; fz = 0; break;
        case 7: fx=-1; fy=-1; fz = 0; break;
        
      }
      nx = fx; ny = fy; nz = dz;
      // switch order if different direction
      if(dx != 0){
        fz = fx; fx = 0;
        nz = fz; nx = dx;
      }if(dy != 0){
        fz = fy; fy = 0;
        nz = fz; ny = dy;
      }
      return;
    case 2:
      if(dx == 0){
        switch(dev)
        {
          case 0:
            fx = 0; fy = 0; fz = -dz;
            nx = 0; ny = dy; nz = -dz;
            return;
          case 1:
            fx = 0; fy = -dy; fz = 0;
            nx = 0; ny = -dy; nz = dz;
            return;
          case 2:
            fx = 1; fy = 0; fz = 0;
            nx = 1; ny = dy; nz = dz;
            return;
          case 3:
            fx = -1; fy = 0; fz = 0;
            nx = -1; ny = dy; nz = dz;
            return;
          case 4:
            fx = 1; fy = 0; fz = -dz;
            nx = 1; ny = dy; nz = -dz;
            return;
          case 5:
            fx = 1; fy = -dy; fz = 0;
            nx = 1; ny = -dy; nz = dz;
            return;
          case 6:
            fx = -1; fy = 0; fz = -dz;
            nx = -1; ny = dy; nz = -dz;
            return;
          case 7:
            fx = -1; fy = -dy; fz = 0;
            nx = -1; ny = -dy; nz = dz;
            return;
          // Extras
          case 8:
            fx = 1; fy = 0; fz = 0;
            nx = 1; ny = dy; nz = 0;
            return;
          case 9:
            fx = 1; fy = 0; fz = 0;
            nx = 1; ny = 0; nz = dz;
            return;
          case 10:
            fx = -1; fy = 0; fz = 0;
            nx = -1; ny = dy; nz = 0;
            return;
          case 11:
            fx = -1; fy = 0; fz = 0;
            nx = -1; ny = 0; nz = dz;
            return;
        }
      }else if(dy == 0){
        switch(dev)
        {
          case 0:
            fx = 0; fy = 0; fz = -dz;
            nx = dx; ny = 0; nz = -dz;
            return;
          case 1:
            fx = -dx; fy = 0; fz = 0;
            nx = -dx; ny = 0; nz = dz;
            return;
          case 2:
            fx = 0; fy = 1; fz = 0;
            nx = dx; ny = 1; nz = dz;
            return;
          case 3:
            fx = 0; fy = -1; fz = 0;
            nx = dx; ny = -1;nz = dz;
            return;
          case 4:
            fx = 0; fy = 1; fz = -dz;
            nx = dx; ny = 1; nz = -dz;
            return;
          case 5:
            fx = -dx; fy = 1; fz = 0;
            nx = -dx; ny = 1; nz = dz;
            return;
          case 6:
            fx = 0; fy = -1; fz = -dz;
            nx = dx; ny = -1; nz = -dz;
            return;
          case 7:
            fx = -dx; fy = -1; fz = 0;
            nx = -dx; ny = -1; nz = dz;
            return;
          // Extras
          case 8:
            fx = 0; fy = 1; fz = 0;
            nx = dx; ny = 1; nz = 0;
            return;
          case 9:
            fx = 0; fy = 1; fz = 0;
            nx = 0; ny = 1; nz = dz;
            return;
          case 10:
            fx = 0; fy = -1; fz = 0;
            nx = dx; ny = -1; nz = 0;
            return;
          case 11:
            fx = 0; fy = -1; fz = 0;
            nx = 0; ny = -1; nz = dz;
            return;
        }
      }else{// dz==0
        switch(dev)
        {
          case 0:
            fx = 0; fy = -dy; fz = 0;
            nx = dx; ny = -dy; nz = 0;
            return;
          case 1:
            fx = -dx; fy = 0; fz = 0;
            nx = -dx; ny = dy; nz = 0;
            return;
          case 2:
            fx =  0; fy = 0; fz = 1;
            nx = dx; ny = dy; nz = 1;
            return;
          case 3:
            fx =  0; fy = 0; fz = -1;
            nx = dx; ny = dy; nz = -1;
            return;
          case 4:
            fx = 0; fy = -dy; fz = 1;
            nx = dx; ny = -dy; nz = 1;
            return;
          case 5:
            fx = -dx; fy = 0; fz = 1;
            nx = -dx; ny = dy; nz = 1;
            return;
          case 6:
            fx = 0; fy = -dy; fz = -1;
            nx = dx; ny = -dy; nz = -1;
            return;
          case 7:
            fx = -dx; fy = 0; fz = -1;
            nx = -dx; ny = dy; nz = -1;
            return;
          // Extras
          case 8:
            fx =  0; fy = 0; fz = 1;
            nx = dx; ny = 0; nz = 1;
            return;
          case 9:
            fx = 0; fy = 0; fz = 1;
            nx = 0; ny = dy; nz = 1;
            return;
          case 10:
            fx =  0; fy = 0; fz = -1;
            nx = dx; ny = 0; nz = -1;
            return;
          case 11:
            fx = 0; fy = 0; fz = -1;
            nx = 0; ny = dy; nz = -1;
            return;
        }
      }
      return;
    case 3:
      switch(dev)
      {
        case 0:
          fx = -dx; fy = 0; fz = 0;
          nx = -dx; ny = dy; nz = dz;
          return;
        case 1:
          fx = 0; fy = -dy; fz = 0;
          nx = dx; ny = -dy; nz = dz;
          return;
        case 2:
          fx = 0; fy = 0; fz = -dz;
          nx = dx; ny = dy; nz = -dz;
          return;
        // Need to check up to here for forced!
        case 3:
          fx = 0; fy = -dy; fz = -dz;
          nx = dx; ny = -dy; nz = -dz;
          return;
        case 4:
          fx = -dx; fy = 0; fz = -dz;
          nx = -dx; ny = dy; nz = -dz;
          return;
        case 5:
          fx = -dx; fy = -dy; fz = 0;
          nx = -dx; ny = -dy; nz = dz;
          return;
        // Extras
        case 6:
          fx = -dx; fy = 0; fz = 0;
          nx = -dx; ny = 0; nz = dz;
          return;
        case 7:
          fx = -dx; fy = 0; fz = 0;
          nx = -dx; ny = dy; nz = 0;
          return;
        case 8:
          fx = 0; fy = -dy; fz = 0;
          nx = 0; ny = -dy; nz = dz;
          return;
        case 9:
          fx = 0; fy = -dy; fz = 0;
          nx = dx; ny = -dy; nz = 0;
          return;
        case 10:
          fx = 0; fy = 0; fz = -dz;
          nx = 0; ny = dy; nz = -dz;
          return;
        case 11:
          fx = 0; fy = 0; fz = -dz;
          nx = dx; ny = 0; nz = -dz;
          return;
      }
      return;
  }
    }
};

template <class T>
requires std::derived_from<T, IGridMap>
class GraphSearch{
public:
inline GraphSearch(const T& m):map(m){
}
static inline std::vector<JpsNode> nodes{};
static inline std::unordered_map<long long, size_t> hash_points;
static inline JPS3DNeib jn3d_{};
static inline int xGoal_,yGoal_,zGoal_;
const T& map;


inline bool hasForced(int x, int y, int z, int dx, int dy, int dz) {
  int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
  int id = (dx+1)+3*(dy+1)+9*(dz+1);
  switch(norm1)
  {
    case 1:
      // 1-d move, check 8 neighbors
      for( int fn = 0; fn < 8; ++fn )
      {
        int nx = x + jn3d_.f1[id][0][fn];
        int ny = y + jn3d_.f1[id][1][fn];
        int nz = z + jn3d_.f1[id][2][fn];
        Eigen::Vector3i t{nx,ny,nz};
        if( map.check_point_i(t) == IMap::State::Unsafe )
          return true;
      }
      return false;
    case 2:
      // 2-d move, check 8 neighbors
      for( int fn = 0; fn < 8; ++fn )
      {
        int nx = x + jn3d_.f1[id][0][fn];
        int ny = y + jn3d_.f1[id][1][fn];
        int nz = z + jn3d_.f1[id][2][fn];
        Eigen::Vector3i t{nx,ny,nz};
        if( map.check_point_i(t) == IMap::State::Unsafe )
          return true;
      }
      return false;
    case 3:
      // 3-d move, check 6 neighbors
      for( int fn = 0; fn < 6; ++fn )
      {
        int nx = x + jn3d_.f1[id][0][fn];
        int ny = y + jn3d_.f1[id][1][fn];
        int nz = z + jn3d_.f1[id][2][fn];
        Eigen::Vector3i t{nx,ny,nz};
        if( map.check_point_i(t) == IMap::State::Unsafe )
          return true;
      }
      return false;
    default:
      return false;
  }
}


inline bool jump(int x, int y, int z, int dx, int dy, int dz, int& new_x, int& new_y, int& new_z,int dim) {
  new_x = x + dx;
  new_y = y + dy;
  new_z = z + dz;
  Eigen::Vector3i t{new_x, new_y, new_z};
  if(dim > jump_size)
    return true;
  if (map.check_point_i(t) == IMap::State::Unsafe)
    return false;

  if (new_x ==  xGoal_ && new_y == yGoal_ && new_z == zGoal_)
    return true;

  if (hasForced(new_x, new_y, new_z, dx, dy, dz))
    return true;

  const int id = (dx+1)+3*(dy+1)+9*(dz+1);
  const int norm1 = std::abs(dx) + std::abs(dy) +std::abs(dz);
  int num_neib = jn3d_.nsz[norm1][0];
  for( int k = 0; k < num_neib-1; ++k )
  {
    int new_new_x, new_new_y, new_new_z;
    if(jump(new_x,new_y,new_z,
          jn3d_.ns[id][0][k], jn3d_.ns[id][1][k], jn3d_.ns[id][2][k],
        new_new_x, new_new_y, new_new_z,dim+1)) return true;
  }


  return jump(new_x, new_y, new_z, dx, dy, dz, new_x, new_y, new_z,dim+1);
}
inline void get_jps_succ(
    JpsNode* curr, 
    std::vector<Eigen::Vector3i>& succ_ids, 
    std::vector<double>& succ_costs
    ) {
    const int norm1 = std::abs(curr->dir(0))+std::abs(curr->dir(1))+std::abs(curr->dir(2));
    int num_neib = jn3d_.nsz[norm1][0];
    int num_fneib = jn3d_.nsz[norm1][1];
    int id = (curr->dir(0)+1)+3*(curr->dir(1)+1)+9*(curr->dir(2)+1);

    for( int dev = 0; dev < num_neib+num_fneib; ++dev) {
        int new_x, new_y, new_z;
        int dx, dy, dz;
        if(dev < num_neib) {
        // std::cerr << id << ' ' << dev<< std::endl;
        dx = jn3d_.ns[id][0][dev];
        dy = jn3d_.ns[id][1][dev];
        dz = jn3d_.ns[id][2][dev];
        if(!jump(curr->point.x(), curr->point.y(), curr->point.z(),
                dx, dy, dz, new_x, new_y, new_z,0)) continue;
        }
        else {
        int nx = curr->point(0) + jn3d_.f1[id][0][dev-num_neib];
        int ny = curr->point(1) + jn3d_.f1[id][1][dev-num_neib];
        int nz = curr->point(2) + jn3d_.f1[id][2][dev-num_neib];
        Eigen::Vector3i t{nx,ny,nz};
        if(map.check_point_i(t) == IMap::State::Unsafe) {
            dx = jn3d_.f2[id][0][dev-num_neib];
            dy = jn3d_.f2[id][1][dev-num_neib];
            dz = jn3d_.f2[id][2][dev-num_neib];
            if(!jump(curr->point.x(), curr->point.y(), curr->point.z(),
                dx, dy, dz, new_x, new_y, new_z,0)) continue;
        }
        else
            continue;
        }
        succ_ids.emplace_back(new_x,new_y,new_z);
        succ_costs.push_back(std::sqrt(static_cast<double>((succ_ids.back() - curr->point).squaredNorm())));
    }

}

struct NodeCompare {
    bool operator()(const JpsNode* a1, const JpsNode* a2) const {
      double f1 = a1->g + a1->h;
      double f2 = a2->g + a2->h;
      if( ( f1 >= f2 - 0.000001) && (f1 <= f2 +0.000001) )
        return a1->g < a2->g; // if equal compare gvals
      return f1 > f2;
    }
};

inline void graph_search(
    const Eigen::Vector3i& start,
    const Eigen::Vector3i& goal,
    std::vector<Eigen::Vector3f>&   path_out)
{
    int max_expand = 0;
    
    std::priority_queue<JpsNode*,std::vector<JpsNode*>,NodeCompare> pq;
    nodes.clear();
    hash_points.clear();
    Eigen::Vector3i map_lower,map_upper;
    map.get_map_bound_i(map_lower,map_upper);
    const auto dim = map_upper - map_lower;
    const auto hash = [&](const Eigen::Vector3i& p) -> long long
    {
        return  p(0) + p(1) * dim(0) + p(2) * dim(1) * dim(2);
    };
	nodes.reserve(dim.prod());
    const auto start_id = hash(start);
    const auto goal_id = hash(goal);
    xGoal_ = goal.x();
    yGoal_ = goal.y();
    zGoal_ = goal.z();
    nodes.emplace_back(start,start_id);
    nodes.back().open = true;
    nodes.back().g = 0;
    nodes.back().dir = {0,0,0};
    hash_points.insert_or_assign(start_id,nodes.size() - 1);
    pq.emplace(&nodes.back());
    static std::vector<Eigen::Vector3i> succ_ids{};
    static std::vector<double> succ_costs{};
    while(!pq.empty())
    {
        max_expand++;
        auto curr = pq.top(); pq.pop();
        if(curr->id == goal_id) {
            break;
        }
        if(curr->close) continue;
        curr->close = true;
        succ_costs.clear();
        succ_ids.clear();
        get_jps_succ(curr, succ_ids, succ_costs);
        for( int s = 0, l = succ_ids.size(); s < l; ++s )
        {
            auto id = hash(succ_ids[s]);
            if(!hash_points.contains(id))
            {
                nodes.emplace_back(succ_ids[s],id);
                hash_points.insert_or_assign(id,nodes.size() - 1);
                nodes.back().h = std::sqrt(static_cast<double>((goal - succ_ids[s]).squaredNorm()));
            }
            
            JpsNode* child_ptr = &nodes[hash_points[id]];
            double tentative_gval = curr->g + succ_costs[s];
			if(child_ptr->open && child_ptr->close) continue;
            
            if(tentative_gval >= child_ptr->g ) continue;
            
			child_ptr->parent_id = curr->id;  
			child_ptr->g = tentative_gval;    
			
			child_ptr->dir = (child_ptr->point - curr->point);
			if(child_ptr->dir(0) != 0)
				child_ptr->dir(0) /= std::abs(child_ptr->dir(0));
			if(child_ptr->dir(1) != 0)
				child_ptr->dir(1) /= std::abs(child_ptr->dir(1));
			if(child_ptr->dir(2) != 0)
				child_ptr->dir(2) /= std::abs(child_ptr->dir(2));
			
			pq.emplace(child_ptr);
			child_ptr->open = true;
            
        }
        if(max_expand > 10000)break;
    }
	if(max_expand > 10000 || !hash_points.contains(goal_id)){
		Eigen::Vector3f p;
		map.index_to_pos(start,p);
		path_out.emplace_back(p);
		map.index_to_pos(goal,p);
		path_out.emplace_back(p);
		return;
	}
    JpsNode* node = &nodes[hash_points[goal_id]];
    Eigen::Vector3f np;
	map.index_to_pos(node->point,np);
	path_out.emplace_back(np);
    while(hash_points.contains(node->parent_id) && node->parent_id != start_id) {
		node = &nodes[hash_points[node->parent_id]];
		map.index_to_pos(node->point,np);
		path_out.emplace_back(np);
    }
	node = &nodes[hash_points[node->parent_id]];
	map.index_to_pos(node->point,np);
	path_out.emplace_back(np);
}
};
}


template <class T>
requires std::derived_from<T, IGridMap>
static inline std::vector<Eigen::Vector3f> remove_corner_pts(const std::vector<Eigen::Vector3f> &path,const T& map) {
  if (path.size() < 2)
    return path;

  // cut zigzag segment
  std::vector<Eigen::Vector3f> optimized_path;
  Eigen::Vector3f pose1 = path[0];
  Eigen::Vector3f pose2 = path[1];
  Eigen::Vector3f prev_pose = pose1;
  optimized_path.push_back(pose1);
  double cost1, cost2, cost3;

  if (map.check_line_d(pose1, pose2) == IMap::State::Safe)
    cost1 = (pose1 - pose2).norm();
  else
    cost1 = std::numeric_limits<double>::infinity();

  for (unsigned int i = 1; i < path.size() - 1; i++) {
    pose1 = path[i];
    pose2 = path[i + 1];
    if (map.check_line_d(pose1, pose2) == IMap::State::Safe)
      cost2 = (pose1 - pose2).norm();
    else
      cost2 = std::numeric_limits<double>::infinity();

    if (map.check_line_d(prev_pose, pose2) == IMap::State::Safe)
      cost3 = (prev_pose - pose2).norm();
    else
      cost3 = std::numeric_limits<double>::infinity();

    if (cost3 < cost1 + cost2)
      cost1 = cost3;
    else {
      optimized_path.push_back(path[i]);
      cost1 = (pose1 - pose2).norm();
      prev_pose = pose1;
    }
  }

  optimized_path.push_back(path.back());
  return optimized_path;
}



static inline std::vector<Eigen::Vector3f> remove_line_pts(const std::vector<Eigen::Vector3f> &path) {
  if (path.size() < 3)
    return path;

  std::vector<Eigen::Vector3f> new_path;
  new_path.push_back(path.front());
  for (unsigned int i = 1; i < path.size() - 1; i++) {
    Eigen::Vector3f p = (path[i + 1] - path[i]) - (path[i] - path[i - 1]);
    
    if (fabs(p(0)) + fabs(p(1)) > 1e-2)
    new_path.push_back(path[i]);
    
  }
  new_path.push_back(path.back());
  return new_path;
}



template <class T>
requires std::derived_from<T, IGridMap>
static inline void jps_search(
    const Eigen::Vector3f& start,
    const Eigen::Vector3f& goal,
    const T& map,
    std::vector<Eigen::Vector3f>&   path_out){
    jps::GraphSearch graph_search{map};
    path_out.clear();
    Eigen::Vector3i start_index, goal_index;
    map.pos_to_index(start,start_index);
    map.pos_to_index(goal,goal_index);
    if(map.check_point_i(start_index) != IMap::State::Safe)
        return ;
    if(map.check_point_i(goal_index) != IMap::State::Safe)
        return ;
    graph_search.graph_search(start_index, goal_index, path_out);
    path_out = remove_corner_pts(path_out, map);
    std::reverse(std::begin(path_out), std::end(path_out));
    path_out = remove_corner_pts(path_out, map);
    path_out = remove_line_pts(path_out);        
}

}
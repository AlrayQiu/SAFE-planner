#pragma once

#include <cstdint>
namespace safe_planner::planner::trajectory::utils{



template<int l, typename VecT>
inline void power_base(const float& v,const uint8_t& p,VecT& out){

    for(int i = 0; i < l; ++i){
    if(i < p){
        out(i) = 0;
    }
    else if(i == p){
        out(i) = 1;
        for(int j = 2;j <= i;++j)
        out(i) *= j;
    } else {
        out(i) = out(i - 1) * v / (i - p) * i;
    }
}
} 

}
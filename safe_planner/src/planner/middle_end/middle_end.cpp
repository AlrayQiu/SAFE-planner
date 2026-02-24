#include "safe_planner/planner/middle_end.hpp"
#include "optimizer/minco.hpp"
#include "safe_planner/map/imap.hpp"
#include "safe_planner/map/rog_map.hpp"
#include "safe_planner/trajectory/minco.hpp"
#include "optimizer/lbfgs.hpp"
#include <vector>

namespace safe_planner::planner{

#define CLASS(x)  template<class TMap, class TTraj> \
               requires std::derived_from<TMap, IMap> && std::derived_from<TTraj, ITrajectory> \
               x MiddleEnd<TMap, TTraj>


CLASS(class)::Impl
{
public:
Impl(const TMap& map)
: map_(map)
{}
void optimize(
    const std::vector<Eigen::Vector3f>& ref_path,
    const Eigen::Matrix<float,2,3>& begin_va,
    const Eigen::Matrix<float,2,3>& end_va,
    TTraj& traj){
        std::vector<trajectory::Minco::seconds> control_seconds(ref_path.size());
        for(size_t i = 0; i < control_seconds.size(); ++i){
            control_seconds[i] = trajectory::Minco::seconds{i};
        }
        trajectory::Minco minco(ref_path,control_seconds,begin_va,end_va);
        optimizable::MincoOpt opt(minco);
        Eigen::VectorXf x;
        opt.init_x(x);
        float fx;
        /* Set the minimization parameters */
        lbfgs::lbfgs_parameter_t params;
        params.g_epsilon = 1.0e-8;
        params.past = 5;
        params.delta = 1.0e-8;
        params.max_iterations = 50;
        optdata data{this, &opt, ref_path};
        /* Start minimization */
        auto ret = lbfgs::lbfgs_optimize(x,
                            fx,
                            cost_function,
                            nullptr,
                            monitor_progress,
                            &data,
                            params);
        
        if(ret >= 0) opt.set_x(x);
        std::cerr << lbfgs::lbfgs_strerror(ret) << std::endl;
        if(ret == lbfgs::LBFGSERR_INVALID_FUNCVAL) 
            std::cerr << x.transpose() << std::endl;
        traj = minco;
    }

private:


struct optdata
{
    MiddleEnd<TMap, TTraj>::Impl* instance;
    void* opt;
    const std::vector<Eigen::Vector3f>& ref_path;
};

static int monitor_progress(void *,
                               const Eigen::VectorXf &x,
                               const Eigen::VectorXf &g,
                               const float fx,
                               const float,
                               const int k,
                               const int )
{
    std::cout << std::setprecision(4)
                << "================================" << std::endl
                << "Iteration: " << k << std::endl
                << "Function Value: " << fx << std::endl
                << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
                << "Variables: " << std::endl
                << x.transpose() << std::endl;
    return 0;
}
static float cost_function(void *instance,
                               const Eigen::VectorXf &x,
                               Eigen::VectorXf &grad)
{
    optdata* data = static_cast<optdata*>(instance);
    auto opt = static_cast<optimizable::MincoOpt*>(data->opt);
    opt->set_x(x);
    opt->clear_j_g();
    opt->enum_p(0.05,[&](const int i, const float ratio, const Eigen::Vector3f& p){
        data->instance->check_point(i, ratio,p,data->ref_path,opt);
    });
    opt->get_g(x, grad);
    float j;
    opt->get_j(j);
    return j;
        
}
void check_point(
    const int i, const float& ratio,
    const Eigen::Vector3f& p, 
    const std::vector<Eigen::Vector3f>& ref_path,
    optimizable::MincoOpt* opt){

    if(map_.check_point_d(p) != IMap::State::Unsafe) return;
    Eigen::Vector3f g{};
    
    float l;
    point_line_dis_squared(p, ref_path[i], ref_path[i + 1], l, g);
    opt->add_j(l * distance_weight_) ;
    opt->add_g(ratio, i,g * distance_weight_);
}
void point_line_dis_squared(const Eigen::Vector3f& p, const Eigen::Vector3f& a, const Eigen::Vector3f& b, float& l, Eigen::Vector3f& g){
    Eigen::Vector3f u = b - a; 
    Eigen::Vector3f u_norm = u.normalized(); 
    Eigen::Vector3f ap = p - a; 
    float proj = ap.dot(u_norm);
    Eigen::Vector3f n = proj * u_norm - ap;
    l = n.squaredNorm();
    g = 2 * n.normalized();
}
const TMap& map_;
const float distance_weight_ = 100.0f;
};

CLASS()::MiddleEnd(const TMap& map){
    impl_ = std::make_unique<Impl>(map);
}

CLASS()::~MiddleEnd() = default;

CLASS(void)::optimize(
    const std::vector<Eigen::Vector3f>& ref_path,
    const Eigen::Matrix<float,2,3>& begin_va,
    const Eigen::Matrix<float,2,3>& end_va,
    TTraj& traj){
        impl_->optimize(ref_path, begin_va, end_va, traj);
}

template class MiddleEnd<safe_planner::map::ROGMap, safe_planner::planner::trajectory::Minco>;
}
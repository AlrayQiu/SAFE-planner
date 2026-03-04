#include "safe_planner/planner/middle_end.hpp"
#include "optimizer/lbfgs.hpp"
#include "optimizer/parameterized_urbs_optmizer.hpp"
#include "safe_planner/esdf/implicit_esdf.hpp"
#include "safe_planner/map/imap.hpp"
#include "safe_planner/map/rog_map.hpp"
#include "safe_planner/trajectory/bspline.hpp"
#include "trajectory/uniform_bspline_impl.hpp"
#include <memory>
#include <vector>

namespace safe_planner::planner{

#define CLASS(x)  template<class TMap> \
               requires std::derived_from<TMap, IMap>\
               x MiddleEnd<TMap>


CLASS(class)::Impl
{
public:
Impl(const TMap& map, const esdf::ImplicitESDF<TMap>& esdf)
: map_(map)
, esdf_(esdf)
{}
void optimize(
    const std::vector<Eigen::Vector3d>& ref_path,
    const Eigen::Matrix<double,3,2>& begin_va,
    const Eigen::Matrix<double,3,2>& end_va,
    trajectory::UniformBSpline& traj){
        trajectory::impl::UniformBSplineImpl urbs{ref_path,begin_va,end_va};
        middle_end::optimizer::ParameterizedUrbsOpt<TMap> opt(urbs,esdf_,map_);
        Eigen::VectorXd x;
        double fx;
        opt.init_x(x);
        lbfgs::lbfgs_parameter_t params;
        params.mem_size = 256;
        params.g_epsilon = 1.0e-8;
        params.past = 5;
        params.delta = 1.0e-6;
        params.max_iterations = 200;
        optdata data{this, &opt};
        // /* Start minimization */
        if(x.size() <= 1){
            urbs.build(traj);
            return;
        }
        auto ret = lbfgs::lbfgs_optimize(x,
                            fx,
                            cost_function,
                            nullptr,
                            nullptr,
                            &data,
                            params);
        opt.set_x(x);
        urbs.build(traj);
        std::cerr << lbfgs::lbfgs_strerror(ret) << std::endl;
    }
void test_optimizer(
    trajectory::MultiPoly& ,
    std::vector<Eigen::Vector3d>& ){
}
private:


struct optdata
{
    MiddleEnd<TMap>::Impl* instance;
    void* opt;
};

static int monitor_progress(void *,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double,
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

static double cost_function(void *instance,
                               const Eigen::VectorXd &x,
                               Eigen::VectorXd &grad)
{
    auto data = static_cast<optdata*>(instance);
    auto opt = static_cast<middle_end::optimizer::ParameterizedUrbsOpt<TMap>*>(data->opt);

    double j;
    opt->set_x(x);
    opt->get_g(grad);
    opt->get_j(j);
    return j;
}
const TMap& map_;
const esdf::ImplicitESDF<TMap>& esdf_;
};

CLASS()::MiddleEnd(const TMap& map, const esdf::ImplicitESDF<TMap>& esdf){
    impl_ = std::make_unique<Impl>(map, esdf);
}

CLASS()::~MiddleEnd() = default;

CLASS(void)::optimize(
    const std::vector<Eigen::Vector3d>& ref_path,
    const Eigen::Matrix<double,3,2>& begin_va,
    const Eigen::Matrix<double,3,2>& end_va,
    trajectory::UniformBSpline& traj){
        impl_->optimize(ref_path, begin_va, end_va, traj);
}
CLASS(void)::test_optimizer(
        trajectory::MultiPoly& traj,
        std::vector<Eigen::Vector3d>& gradients){
            impl_->test_optimizer(traj, gradients);
        }

template class MiddleEnd<safe_planner::map::ROGMap>;
// template class MiddleEnd<safe_planner::map::ROGMap, safe_planner::planner::trajectory::MincoTest>;
}
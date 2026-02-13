
#include "safe_planner/map/imap.hpp"
#include <memory>

namespace safe_planner::planner
{
namespace front_end{

struct Config{
enum class Algo{
    RRT = 0
};
Algo algo;
};

}
template<class TMap>
requires std::derived_from<TMap, IMap>
class FrontEnd{
public:
    FrontEnd(const TMap&,front_end::Config&& config);
    FrontEnd() = delete;
    ~FrontEnd();
    void search(
        const Eigen::Vector3f& from,
        const Eigen::Vector3f& to,
        std::vector<Eigen::Vector3f>& path);
private:
    const TMap& map_;
    class FrontEndImpl;
    std::unique_ptr<FrontEndImpl> impl_;
};
}
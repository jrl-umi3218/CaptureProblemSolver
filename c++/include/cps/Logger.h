#include <Eigen/Core>

#include <map>
#include <string>
#include <vector>

#include <cps/cps_api.h>

namespace cps
{
  class CPS_DLLAPI Logger
  {
  public:
    Logger(const std::string& filename);
  };
}

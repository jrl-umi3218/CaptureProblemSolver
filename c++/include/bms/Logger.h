#include <Eigen/Core>

#include <map>
#include <string>
#include <vector>

#include <bms/bms_api.h>

namespace cps
{
  class BMS_DLLAPI Logger
  {
  public:
    Logger(const std::string& filename);
  };
}

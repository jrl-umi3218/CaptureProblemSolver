#include <Eigen/Core>

#include <map>
#include <string>
#include <vector>

#include <bms_api.h>

namespace bms
{
  class BMS_DLLAPI Logger
  {
  public:
    Logger(const std::string& filename);
  };
}
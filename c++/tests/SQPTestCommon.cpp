#include "SQPTestCommon.h"

namespace bms
{
  RawProblem resampleProblem(const RawProblem& raw, int n)
  {
    auto raw2 = raw;
    VectorXd s = VectorXd::LinSpaced(n + 1, 0, 1);
    raw2.delta.resize(n);
    raw2.delta.array() = s.tail(n).array()*s.tail(n).array() - s.head(n).array()*s.head(n).array();

    //we invalidate the solution
    raw2.Phi_.resize(0);

    return raw2;
  }

  RawProblem resampleProblem(const std::string& filepath, int n)
  {
    RawProblem raw;
    std::string base = TESTS_DIR;
    raw.read(base + "/" + filepath);

    return resampleProblem(raw, n);
  }
}
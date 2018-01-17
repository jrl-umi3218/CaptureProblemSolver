#include <cps/Problem.h>
#include <cps/SQP.h>

using namespace Eigen;

namespace cps
{
  RawProblem resampleProblem(const RawProblem& raw, int n);
  RawProblem resampleProblem(const std::string& filepath, int n);
}

#include <cps/QRAlgorithms.h>

using namespace Eigen;

namespace cps
{
  SpecialQR::SpecialQR(DenseIndex kmax)
    : e_(kmax + 1)
    , c1_(kmax)
    , c2_(kmax)
    , d_(kmax)
    , l_(kmax)
    , dprecomp_(-VectorXd::LinSpaced(static_cast<int>(kmax),1, static_cast<int>(kmax)))
    , lprecomp_((ArrayXd::LinSpaced(static_cast<int>(kmax), 1, static_cast<int>(kmax)).sqrt()*ArrayXd::LinSpaced(static_cast<int>(kmax), 2, static_cast<int>(kmax) + 1).sqrt()).inverse().matrix())
  {
    G_.reserve(kmax);
    for (DenseIndex i = 0; i < kmax; ++i)
    {
      double c = 1 / sqrt(double(i) + 2.);
      G_.emplace_back(i, c, sqrt(double(i) + 1.)*c);
    }
  }
}

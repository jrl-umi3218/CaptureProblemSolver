#include "QRAlgorithms.h"

using namespace Eigen;

namespace bms
{
  SpecialQR::SpecialQR(int kmax)
    : e_(kmax + 1)
    , c1_(kmax)
    , c2_(kmax)
    , d_(kmax)
    , l_(kmax)
    , dprecomp_(-VectorXd::LinSpaced(kmax,1,kmax))
    , lprecomp_((ArrayXd::LinSpaced(kmax, 1, kmax).sqrt()*ArrayXd::LinSpaced(kmax, 2, kmax + 1).sqrt()).inverse().matrix())
  {
    G_.reserve(kmax);
    for (int i = 0; i < kmax; ++i)
    {
      double c = 1 / sqrt(i + 2.);
      G_.emplace_back(i, c, sqrt(i + 1.)*c);
    }
  }
}
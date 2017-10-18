#include "GivensSequence.h"

namespace bms
{
  void GivensSequence::extend(int incr)
  {
    for (auto& G : *this)
      G.extend(incr);
  }

  Eigen::MatrixXd GivensSequence::matrix(int n)
  {
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(n,n);
    for (const auto& G : *this)
      G.applyOnTheRightTo(Q);

    return Q;
  }
}
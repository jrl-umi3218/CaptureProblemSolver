#include "CondensedOrthogonalMatrix.h"

using namespace Eigen;

namespace bms
{
  CondensedOrthogonalMatrix::CondensedOrthogonalMatrix(int n, int kmax, int pmax)
    : n_(n)
    , sequences_(kmax)
    , transpositions_(VectorXi::LinSpaced(n,0,n-1))
  {
    for (auto& s : sequences_)
      s.reserve(pmax);
  }

  void CondensedOrthogonalMatrix::reset()
  {
    for (auto& s : sequences_)
      s.clear();
    transpositions_.indices() = VectorXi::LinSpaced(n_, 0, static_cast<int>(n_) - 1);
  }

  Eigen::MatrixXd CondensedOrthogonalMatrix::matrix()
  {
    MatrixXd M = MatrixXd::Identity(n_, n_);
    for (const auto& s : sequences_)
      s.applyOnTheRightTo(M);
    M = M*transpositions_;

    return M;
  }
}
#include "CondensedOrthogonalMatrix.h"

using namespace Eigen;

namespace bms
{
  CondensedOrthogonalMatrix::CondensedOrthogonalMatrix(int n, int kmax, int pmax, bool Ptranspose)
    : ptranspose_(Ptranspose)
    , n_(n)
    , sequences_(kmax)
    , transpositions_(VectorXi::LinSpaced(n,0,n-1))
  {
    for (auto& s : sequences_)
      s.reserve(pmax);
    Qh_.reserve(pmax);
  }

  void CondensedOrthogonalMatrix::reset(bool Ptranspose)
  {
    for (auto& s : sequences_)
      s.clear();
    Qh_.clear();
    transpositions_.indices() = VectorXi::LinSpaced(n_, 0, static_cast<int>(n_) - 1);
    ptranspose_ = Ptranspose;
  }

  Eigen::MatrixXd CondensedOrthogonalMatrix::matrix()
  {
    MatrixXd M = MatrixXd::Identity(n_, n_);
    this->applyOnTheRightTo(M);

    return M;
  }
}
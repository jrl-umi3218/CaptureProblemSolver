#include <cps/CondensedOrthogonalMatrix.h>

using namespace Eigen;

namespace cps
{
  CondensedOrthogonalMatrix::CondensedOrthogonalMatrix(int n, int kmax, int pmax, bool Ptranspose)
    : ptranspose_(Ptranspose)
    , n_(n)
    , sequences_(kmax)
    , transpositions_(VectorXi::LinSpaced(n,0,n-1))
    , identityTransposition_(VectorXi::LinSpaced(n, 0, n - 1))
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
    transpositions_.indices() = identityTransposition_;
    ptranspose_ = Ptranspose;
  }

  void CondensedOrthogonalMatrix::resize(int n, int kmax, int pmax)
  {
    n_ = n;
    sequences_.resize(kmax);
    for (auto& s : sequences_)
      s.reserve(pmax);
    Qh_.reserve(pmax);
    transpositions_.indices() = VectorXi::LinSpaced(n, 0, n - 1);
    identityTransposition_ = VectorXi::LinSpaced(n, 0, n - 1);
  }

  Eigen::MatrixXd CondensedOrthogonalMatrix::matrix()
  {
    MatrixXd M = MatrixXd::Identity(n_, n_);
    this->applyOnTheRightTo(M);

    return M;
  }
}

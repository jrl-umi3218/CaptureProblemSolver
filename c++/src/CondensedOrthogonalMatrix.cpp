/* Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
 *
 * This file is part of CPS.
 *
 * CPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CPS.  If not, see <http://www.gnu.org/licenses/>.
 */

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

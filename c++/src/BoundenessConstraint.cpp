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

#include <cps/BoundenessConstraint.h>

using namespace Eigen;

namespace cps
{
  BoundenessConstraint::BoundenessConstraint(const Eigen::VectorXd& delta, double alpha, double b)
    : n_(delta.size()), alpha_(alpha), b_(b), delta_(delta), y_(delta.size()+1)
  {
    y_[0] = 0;
  }

  void BoundenessConstraint::compute(double& val, const VectorConstRef& x) const
  {
    y_.tail(n_) = x.cwiseSqrt();
    val = -alpha_*y_[n_] - b_;
    for (DenseIndex i = 0; i<n_; ++i)
      val += delta_[i] / (y_[i + 1] + y_[i]);
  }

  void BoundenessConstraint::compute(double& val, VectorRef grad, const VectorConstRef& x) const
  {
    compute(val, x);
    grad[0] = -delta_[0] / (2 * y_[1] * y_[1] * y_[1]);
    for (DenseIndex i = 1; i < n_; ++i)
    {
      double d = y_[i + 1] + y_[i];
      double d2 = d*d;
      grad[i-1] -= delta_[i] / (2 * y_[i] * d2);
      grad[i] = -delta_[i] / (2 * y_[i+1] * d2);
    }
    grad[n_ - 1] -= alpha_ / (2 * y_[n_]);
  }

  void BoundenessConstraint::setAlpha(double alpha)
  {
    alpha_ = alpha;
  }

  void BoundenessConstraint::setb(double b)
  {
    b_ = b;
  }
}

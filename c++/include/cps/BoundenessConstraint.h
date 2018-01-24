#pragma once
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

#include <Eigen/Core>

#include <cps/cps_api.h>
#include <cps/defs.h>

namespace cps
{
  /** The function sum(delta_j/(sqrt(x(j+1))+sqrt(x(j))) - alpha sqrt(x(n)) - b*/
  class CPS_DLLAPI BoundenessConstraint
  {
  public:
    BoundenessConstraint(const Eigen::VectorXd& delta, double alpha, double b);

    void compute(double& val, const VectorConstRef& x) const;
    void compute(double& val, VectorRef grad, const VectorConstRef& x) const;

    // alternative version
    double compute(const VectorConstRef& x) const;

    void setAlpha(double alpha);
    void setb(double b);

  private:
    Eigen::DenseIndex n_;
    double alpha_;
    double b_;
    Eigen::VectorXd delta_;

    //intermediate data
    mutable Eigen::VectorXd y_;
  };
}

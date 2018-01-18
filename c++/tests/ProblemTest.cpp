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

#include <iostream>

#include <cps/Problem.h>
#include <cps/BoundenessConstraint.h>
#include <cps/SQP.h>

// boost
#define BOOST_TEST_MODULE ProblemTest
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace cps;

BOOST_AUTO_TEST_CASE(BoundenessConstraintTest)
{
  VectorXd delta = VectorXd::LinSpaced(10, 0.01, 0.19);
  BoundenessConstraint bc(delta, 1, 1);

  VectorXd x = VectorXd::Ones(10);
  VectorXd grad0(10);
  grad0 << -0.00875, -0.01, -0.015, -0.02, -0.025, -0.03, -0.035, -0.04, -0.045, -0.52375;

  double y;
  VectorXd grad(10);
  bc.compute(y, grad, x);

  BOOST_CHECK(std::abs(y + 1.495) <= 1e-15);
  BOOST_CHECK(grad.isApprox(grad0, 1e-15));
}

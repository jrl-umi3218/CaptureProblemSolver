#include <iostream>

#include <bms/Problem.h>
#include <bms/BoundenessConstraint.h>
#include <bms/SQP.h>

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

#include <Problem.h>
#include <BoundenessConstraint.h>
#include <SQP.h>

#include <iostream>

// boost
#define BOOST_TEST_MODULE ProblemTest
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace bms;

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

  BOOST_CHECK(abs(y + 1.495) <= 1e-15);
  BOOST_CHECK(grad.isApprox(grad0, 1e-15));
}

BOOST_AUTO_TEST_CASE(testFeasibilitySQP)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/data/Problem01.txt");

  Problem pb(raw);
  SQP sqp(static_cast<int>(raw.delta.size()));

  sqp.solveFeasibility(pb);

  double c;
  pb.nonLinearConstraint().compute(c, sqp.x());
  BOOST_CHECK(pb.linearConstraints().checkPrimal(sqp.x()));
  BOOST_CHECK(std::abs(c) <= 1e-7);
}
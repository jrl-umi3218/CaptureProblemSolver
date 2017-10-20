#include <LinearConstraints.h>
#include <ProblemMatrices.h>

#include<Eigen/SVD>

#include <iostream>

// boost
#define BOOST_TEST_MODULE LinearConstraintsTest
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace bms;

BOOST_AUTO_TEST_CASE(ActivationTest)
{
  VectorXd l = -VectorXd::Random(10).cwiseAbs();
  VectorXd u = VectorXd::Random(10).cwiseAbs();
  u[3] = l[3];
  u[6] = l[6];
  u[7] = l[7];

  LinearConstraints lc(l, u, -1, 1);

  BOOST_CHECK(lc.numberOfActiveConstraints() == 3);
  BOOST_CHECK(lc.activationStatus(0) == Activation::None);
  BOOST_CHECK(lc.activationStatus(1) == Activation::None);
  BOOST_CHECK(lc.activationStatus(2) == Activation::None);
  BOOST_CHECK(lc.activationStatus(3) == Activation::Equal);
  BOOST_CHECK(lc.activationStatus(4) == Activation::None);
  BOOST_CHECK(lc.activationStatus(5) == Activation::None);
  BOOST_CHECK(lc.activationStatus(6) == Activation::Equal);
  BOOST_CHECK(lc.activationStatus(7) == Activation::Equal);
  BOOST_CHECK(lc.activationStatus(8) == Activation::None);
  BOOST_CHECK(lc.activationStatus(9) == Activation::None);
  BOOST_CHECK(lc.activationStatus(10) == Activation::None);
  
  auto act = lc.activeSet();
  BOOST_CHECK(act.size() == 11);
  BOOST_CHECK(!act[0]);
  BOOST_CHECK(!act[1]);
  BOOST_CHECK(!act[2]);
  BOOST_CHECK(act[3]);
  BOOST_CHECK(!act[4]);
  BOOST_CHECK(!act[5]);
  BOOST_CHECK(act[6]);
  BOOST_CHECK(act[7]);
  BOOST_CHECK(!act[8]);
  BOOST_CHECK(!act[9]);
  BOOST_CHECK(!act[10]);

  auto actIdx = lc.activeSetIdx();
  BOOST_CHECK(actIdx[0] == 3);
  BOOST_CHECK(actIdx[1] == 6);
  BOOST_CHECK(actIdx[2] == 7);

  //activating a non activated constraint
  lc.activate(1, Activation::Lower);
  BOOST_CHECK(lc.activationStatus(1) == Activation::Lower);
  BOOST_CHECK(lc.activeSet()[1]);
  BOOST_CHECK(lc.numberOfActiveConstraints() == 4);

  //activating an activated constraint
  lc.activate(1, Activation::Upper);
  BOOST_CHECK(lc.activationStatus(1) == Activation::Upper);
  BOOST_CHECK(lc.activeSet()[1]);
  BOOST_CHECK(lc.numberOfActiveConstraints() == 4);

  //deactivating a constraint with activate
  lc.activate(3, Activation::None);
  BOOST_CHECK(lc.activationStatus(3) == Activation::None);
  BOOST_CHECK(!lc.activeSet()[3]);
  BOOST_CHECK(lc.numberOfActiveConstraints() == 3);

  //deactivating an activated constraint
  lc.deactivate(7);
  BOOST_CHECK(lc.activationStatus(7) == Activation::None);
  BOOST_CHECK(!lc.activeSet()[7]);
  BOOST_CHECK(lc.numberOfActiveConstraints() == 2);

  //deactivating a deactivated constraint
  lc.deactivate(8);
  BOOST_CHECK(lc.activationStatus(8) == Activation::None);
  BOOST_CHECK(!lc.activeSet()[8]);
  BOOST_CHECK(lc.numberOfActiveConstraints() == 2);

  //check that applyNullSpaceOnTheRight correctly computes the active set index
  MatrixXd Y(1, 8);
  lc.applyNullSpaceOnTheRight(Y, MatrixXd(1, 10));
  actIdx = lc.activeSetIdx();
  BOOST_CHECK(actIdx[0] == 1);
  BOOST_CHECK(actIdx[1] == 6);

  VectorXd y(11);
  VectorXd x = VectorXd::Random(2);
  lc.expandActive(y, x);
  BOOST_CHECK(y[1] == x[0]);
  BOOST_CHECK(y[6] == x[1]);
  BOOST_CHECK(y[0] == 0);
  BOOST_CHECK((y.array().segment(2, 4) == 0).all());
  BOOST_CHECK((y.array().tail(4) == 0).all());
}

BOOST_AUTO_TEST_CASE(NullspaceTest)
{
  VectorXd l = -VectorXd::Random(20).cwiseAbs();
  VectorXd u =  VectorXd::Random(20).cwiseAbs();
  LinearConstraints lc(l,u,-1,1);

  for (size_t i = 0; i < 21; ++i)
  {
    if (rand() > RAND_MAX / 2)
      lc.activate(i, Activation::Lower);
  }

  MatrixXd Ca = lc.matrixAct();
  MatrixXd N(20,20-lc.numberOfActiveConstraints());
  lc.applyNullSpaceOnTheRight(N, MatrixXd::Identity(20, 20));

  MatrixXd CaN = Ca*N;
  MatrixXd NtN = N.transpose()*N;
  BOOST_CHECK(CaN.isZero(1e-15));
  BOOST_CHECK(NtN.isDiagonal());

  VectorXd z = VectorXd::Random(20 - lc.numberOfActiveConstraints());
  VectorXd y(20);
  lc.applyNullSpaceOnTheLeft(y, z);
  VectorXd y2 = N*z;
  BOOST_CHECK(y.isApprox(y2, 1e-15));
}

BOOST_AUTO_TEST_CASE(MultiplicationTest)
{
  VectorXd l = -VectorXd::Random(20).cwiseAbs();
  VectorXd u = VectorXd::Random(20).cwiseAbs();
  LinearConstraints lc(l, u, -1, 1);

  MatrixXd C = lc.matrix();
  VectorXd x = VectorXd::Random(20);
  VectorXd y(21);
  lc.mult(y, x);
  VectorXd z = C*x;
  BOOST_CHECK(y.isApprox(z, 1e-15));

  VectorXd q = VectorXd::Random(21);
  VectorXd r(20);
  lc.transposeMult(r, q);
  VectorXd s = C.transpose()*q;
  BOOST_CHECK(r.isApprox(s, 1e-15));
}

BOOST_AUTO_TEST_CASE(PseudoInverseTest)
{
  VectorXd l = -VectorXd::Random(20).cwiseAbs();
  VectorXd u = VectorXd::Random(20).cwiseAbs();
  LinearConstraints lc(l, u, -1, 1);

  for (size_t i = 1; i < 20; ++i)
  {
    if (rand() > RAND_MAX / 2)
      lc.activate(i, Activation::Lower);
  }

  {
    MatrixXd Ca = lc.matrixAct();
    VectorXd b = VectorXd::Random(Ca.cols());
    VectorXd x(Ca.rows());

    lc.pinvTransposeMultAct(x, b);
    VectorXd y = Ca.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    BOOST_CHECK(x.isApprox(y, 1e-15));
  }

  {
    lc.activate(0, Activation::Lower);
    MatrixXd Ca = lc.matrixAct();
    VectorXd b = VectorXd::Random(Ca.cols());
    VectorXd x(Ca.rows());

    lc.pinvTransposeMultAct(x, b);
    VectorXd y = Ca.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    BOOST_CHECK(x.isApprox(y, 1e-15));
  }

  {
    lc.activate(20, Activation::Lower);
    MatrixXd Ca = lc.matrixAct();
    VectorXd b = VectorXd::Random(Ca.cols());
    VectorXd x(Ca.rows());

    lc.pinvTransposeMultAct(x, b);
    VectorXd y = Ca.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    BOOST_CHECK(x.isApprox(y, 1e-15));
  }

  {
    for (size_t i = 0; i < 20; ++i)
      lc.activate(i, Activation::Lower);
    lc.deactivate(20);

    MatrixXd Ca = lc.matrixAct();
    VectorXd b = VectorXd::Random(Ca.cols());
    VectorXd x(Ca.rows());

    lc.pinvTransposeMultAct(x, b);
    VectorXd y = Ca.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    BOOST_CHECK(x.isApprox(y, 1e-14));
  }

  {
    for (size_t i = 1; i < 21; ++i)
      lc.activate(i, Activation::Lower);
    lc.deactivate(0);

    MatrixXd Ca = lc.matrixAct();
    VectorXd b = VectorXd::Random(Ca.cols());
    VectorXd x(Ca.rows());

    lc.pinvTransposeMultAct(x, b);
    VectorXd y = Ca.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    BOOST_CHECK(x.isApprox(y, 1e-14));
  }
}

BOOST_AUTO_TEST_CASE(InitialPoint)
{
  // A point exists (0 is feasible)
  {
    VectorXd l = -VectorXd::Random(5).cwiseAbs();
    VectorXd u = VectorXd::Random(5).cwiseAbs();
    LinearConstraints lc(l, u, -1, 1);

    auto init = lc.initialPoint();
    BOOST_CHECK(init.first == FeasiblePointInfo::Found);
    BOOST_CHECK(lc.checkPrimal(init.second));
  }

  // No feasible point
  {
    VectorXd l = -VectorXd::Random(5).cwiseAbs();
    VectorXd u = VectorXd::Random(5).cwiseAbs();
    LinearConstraints lc(l, u, 1000, 1010);

    auto init = lc.initialPoint();
    BOOST_CHECK(init.first == FeasiblePointInfo::Infeasible);
  }

  // numerical issue with lower bound
  {
    VectorXd l = -VectorXd::Random(5).cwiseAbs();
    VectorXd u = VectorXd::Random(5).cwiseAbs();
    LinearConstraints lc(l, u, l.sum(), l.sum()+1);

    auto init = lc.initialPoint();
    BOOST_CHECK(init.first == FeasiblePointInfo::NumericalWarningLower);
    BOOST_CHECK(lc.checkPrimal(init.second));
  }

  // numerical issue with upper bound
  {
    VectorXd l = -VectorXd::Random(5).cwiseAbs();
    VectorXd u = VectorXd::Random(5).cwiseAbs();
    LinearConstraints lc(l, u, u.sum() - 1, u.sum()-1e-12);

    auto init = lc.initialPoint();
    BOOST_CHECK(init.first == FeasiblePointInfo::NumericalWarningUpper);
    BOOST_CHECK(lc.checkPrimal(init.second));
  }

  // numerical issue with both bounds
  {
    VectorXd l = -VectorXd::Random(5).cwiseAbs();
    VectorXd u = VectorXd::Random(5).cwiseAbs();
    LinearConstraints lc(l, u, l.sum(), u.sum() - 1e-12);

    auto init = lc.initialPoint();
    BOOST_CHECK(init.first == FeasiblePointInfo::NumericalWarningBoth);
    BOOST_CHECK(lc.checkPrimal(init.second));
  }

  // numerical issue with zonotope
  {
    VectorXd l = -1e-12* VectorXd::Random(5).cwiseAbs();
    VectorXd u = 1e-12*VectorXd::Random(5).cwiseAbs();
    LinearConstraints lc(l, u, -1, 1);

    auto init = lc.initialPoint();
    BOOST_CHECK(init.first == FeasiblePointInfo::TooSmallZonotope);
    BOOST_CHECK(lc.checkPrimal(init.second));
  }
}
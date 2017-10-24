#include <LinearConstraints.h>
#include <QuadraticObjective.h>

#include <iostream>
#include <vector>

// boost
#define BOOST_TEST_MODULE ObjectiveTest
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace bms;

//binary form of i
std::vector<bool> toVec(int i, int nbits)
{
  std::vector<bool> b;
  for (int j = nbits-1; j >= 0; --j)
  {
    int a = static_cast<int>(i / (1 << j));
    b.push_back(a!=0);
    i -= a*(1 << j);
  }
  return b;
}

void setActiveSet(LinearConstraints& lc, const std::vector<bool>& act)
{
  for (size_t i = 0; i < act.size(); ++i)
    lc.activate(i, act[i] ? Activation::Equal : Activation::None);
}

int nact(const std::vector<bool>& act)
{
  int n = 0;
  for (size_t i = 0; i < act.size(); ++i)
    n += act[i];
  return n;
}


/** Check if the matrix is a band matrix with lower bandwidth p and upper bandwidth q*/
bool isBandMatrix(const MatrixXd& M, DenseIndex p, DenseIndex q, double prec = 1e-15)
{
  auto m = M.rows();
  auto n = M.cols();

  bool b = true;
  for (DenseIndex i = p + 1; i < m && b; ++i)
    b = M.diagonal(-i).isZero(prec);
  for (DenseIndex i = q + 1; i < m && b; ++i)
    b = M.diagonal(i).isZero(prec);

  return b;
}

BOOST_AUTO_TEST_CASE(valueTest)
{
  VectorXd x = VectorXd::Random(10);
  VectorXd delta = VectorXd::LinSpaced(10, 0.01, 0.19);
  LeastSquareObjective obj(delta);

  MatrixXd J = obj.matrix();
  VectorXd Jx1 = J*x;
  VectorXd Jx2(9);
  obj.applyJToTheLeft(Jx2, x);
  BOOST_CHECK(Jx1.isApprox(Jx2));

  VectorXd JtJx1 = J.transpose()*Jx1;
  VectorXd JtJx2(10);
  obj.applyJTransposeToTheLeft(JtJx2, Jx2);
  BOOST_CHECK(JtJx1.isApprox(JtJx2));
}

BOOST_AUTO_TEST_CASE(projectedMatrixTest)
{
  int N = 12;
  VectorXd delta = VectorXd::LinSpaced(N, 0.01, 0.23);
  LeastSquareObjective obj(delta);
  LinearConstraints lc(N);

  for (int i = 0; i < (1 << (N+1)) - 1; ++i)
  {
    auto act = toVec(i, N+1);
    setActiveSet(lc, act);
    auto na = lc.numberOfActiveConstraints();
    auto J = obj.matrix();
    MatrixXd JN0(N - 1, N - na);
    lc.applyNullSpaceOnTheRight(JN0, J);
    auto JN = obj.projectedMatrix(na, act);

    BOOST_CHECK(JN.isApprox(JN0, 1e-15));
  }
}

BOOST_AUTO_TEST_CASE(qrTest)
{
  int n = 12;
  VectorXd delta = VectorXd::LinSpaced(n, 0.01, 0.02*n-0.01);
  LeastSquareObjective obj(delta);

  for (int i=0; i < (1 << (n + 1)) - 1; ++i)
  {
    auto act = toVec(i, n + 1);
    MatrixXd R(n - 1, n-nact(act));
    CondensedOrthogonalMatrix Q(n-1, n-1, 2*n);
    obj.qr(R, Q, act);

    MatrixXd J = obj.projectedMatrix(act);
    Q.applyTo(J);
    BOOST_CHECK(isBandMatrix(R, 0, 2, 1e-14));
    BOOST_CHECK(R.isApprox(J, 1e-15));
  }
}
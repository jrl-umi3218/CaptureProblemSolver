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

BOOST_AUTO_TEST_CASE(projectedMatrix)
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
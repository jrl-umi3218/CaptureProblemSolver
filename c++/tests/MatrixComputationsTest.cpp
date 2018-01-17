#include <iostream>

#include <cps/CondensedOrthogonalMatrix.h>
#include <cps/Givens.h>
#include <cps/GivensSequence.h>
#include <cps/ProblemMatrices.h>
#include <cps/QRAlgorithms.h>
#include <cps/QuadraticObjective.h>

// boost
#define BOOST_TEST_MODULE MatrixComputationsTest
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace cps;

BOOST_AUTO_TEST_CASE(GivensTest)
{
  MatrixXd M = MatrixXd::Random(5, 5);
  
  //zeroing M(4,0) with row 3
  Givens G1(M, 3, 4, 0);
  G1.applyTo(M);
  BOOST_CHECK(std::abs(M(4,0)) <= 1e-15);

  //zeroing M(3,0) with row 2
  Givens G2(M, 2, 0);
  G2.applyTo(M);
  BOOST_CHECK(std::abs(M(3, 0)) <= 1e-15);
  BOOST_CHECK(std::abs(M(4, 0)) <= 1e-15);

  //hand computed value of c and s to zero M(2,0), using row 0
  double a = M(0, 0);
  double b = M(2, 0);
  double c = a / sqrt(a*a + b*b);
  double s = -b / sqrt(a*a + b*b);
  Givens G3(0, 2, c, s);
  G3.applyTo(M);
  BOOST_CHECK(std::abs(M(2, 0)) <= 1e-15);
  BOOST_CHECK(std::abs(M(3, 0)) <= 1e-15);
  BOOST_CHECK(std::abs(M(4, 0)) <= 1e-15);

  //hand computed value of c and s to zero M(1,0)
  a = M(0, 0);
  b = M(1, 0);
  c = a / sqrt(a*a + b*b);
  s = -b / sqrt(a*a + b*b);
  Givens G4(0, 1, c, s);
  G4.applyTo(M);
  BOOST_CHECK(std::abs(M(1, 0)) <= 1e-15);
  BOOST_CHECK(std::abs(M(2, 0)) <= 1e-15);
  BOOST_CHECK(std::abs(M(3, 0)) <= 1e-15);
  BOOST_CHECK(std::abs(M(4, 0)) <= 1e-15);

  //test on a block
  Givens G5(M.bottomRightCorner(2,4), 0, 0);
  G5.extend(3);
  G5.applyTo(M.rightCols(4));
  BOOST_CHECK(std::abs(M(4, 1)) <= 1e-15);
}


BOOST_AUTO_TEST_CASE(GivensSequenceTest)
{
  MatrixXd M = MatrixXd::Random(5, 5);
  MatrixXd M2 = M;

  Givens G1(M, 3, 0);
  G1.applyTo(M);
  Givens G2(M, 2, 0);
  G2.applyTo(M);
  Givens G3(M, 1, 0);
  G3.applyTo(M);
  Givens G4(M, 0, 0);
  G4.applyTo(M);

  GivensSequence seq = { G1, G2, G3 };
  seq.push_back(G4);
  seq.applyTo(M2);

  BOOST_CHECK(M.isApprox(M2, 1e-15));
}

BOOST_AUTO_TEST_CASE(CondensedOrthogonalMatrixTest)
{
  MatrixXd M = MatrixXd::Random(5, 5);
  MatrixXd M2 = M;

  //Let's make a weird Givens QR algorithm
  Givens G1(M, 0, 3, 0);
  G1.applyTo(M);
  Givens G2(M, 4, 0, 0);
  G2.applyTo(M);
  Givens G3(M, 1, 4, 0);
  G3.applyTo(M);
  Givens G4(M, 2, 1, 0);
  G4.applyTo(M);
  Givens G5(M, 0, 3, 1);
  G5.applyTo(M);
  Givens G6(M, 4, 0, 1);
  G6.applyTo(M);
  Givens G7(M, 1, 4, 1);
  G7.applyTo(M);
  Givens G8(M, 0, 3, 2);
  G8.applyTo(M);
  Givens G9(M, 4, 0, 2);
  G9.applyTo(M);
  Givens G10(M, 0, 3, 3);
  G10.applyTo(M);
  VectorXi idx(5); idx << 0, 1, 0, 0, 3;
  Transpositions<Dynamic> P(idx);
  M = P.transpose()*M;

  CondensedOrthogonalMatrix Q(5, 4, 4);
  Q.Q(0).push_back(G1);
  Q.Q(0).push_back(G2);
  Q.Q(0).push_back(G3);
  Q.Q(0).push_back(G4);
  Q.Q(1).push_back(G5);
  Q.Q(1).push_back(G6);
  Q.Q(1).push_back(G7);
  Q.Q(2).push_back(G8);
  Q.Q(2).push_back(G9);
  Q.Q(3).push_back(G10);
  Q.P().indices() = idx;
  Q.applyTo(M2);

  BOOST_CHECK(M.isApprox(M2, 1e-14));
}

bool isUpperTriangular(const MatrixXd& M)
{
  return M.triangularView<Eigen::StrictlyLower>().toDenseMatrix().isZero(1e-15);
}

/** Check if the matrix is a band matrix with lower bandwidth p and upper bandwidth q*/
bool isBandMatrix(const MatrixXd& M, DenseIndex p, DenseIndex q)
{
  auto m = M.rows();
  auto n = M.cols();

  bool b = true;
  for (DenseIndex i = p+1; i < m && b; ++i)
    b = M.diagonal(-i).isZero(1e-15);
  for (DenseIndex i = q+1; i < m && b; ++i)
    b = M.diagonal(i).isZero(1e-15);

  return b;
}


BOOST_AUTO_TEST_CASE(HessenbergQRTest)
{
  MatrixXd M = MatrixXd::Random(6, 7);
  M.bottomLeftCorner(5, 5).triangularView<Eigen::StrictlyLower>().setZero();  //make the matrix Hessenberg

  MatrixXd M2 = M;
  GivensSequence seq;

  bool b = hessenbergQR(M, seq, false, 1e-7);
  BOOST_CHECK(isUpperTriangular(M));
  BOOST_CHECK(b);

  MatrixXd Q = seq.matrix(6);
  MatrixXd QR = Q*M;
  BOOST_CHECK(QR.isApprox(M2, 1e-15));
}

BOOST_AUTO_TEST_CASE(TridiagonalQRTest)
{
  MatrixXd M = MatrixXd::Random(6, 6);
  M.bottomLeftCorner(5, 5).triangularView<Eigen::StrictlyLower>().setZero();  //make the matrix tridiagonal (set lower part to 0)
  M.topRightCorner(5, 5).triangularView<Eigen::StrictlyUpper>().setZero();  //make the matrix tridiagonal (set upper part to 0)

  MatrixXd M2 = M;
  GivensSequence seq;

  bool b = tridiagonalQR(M, seq, false, 1e-7);
  BOOST_CHECK(isBandMatrix(M,0,2));
  BOOST_CHECK(b);

  MatrixXd Q = seq.matrix(6);
  MatrixXd QR = Q*M;
  BOOST_CHECK(QR.isApprox(M2, 1e-15));
}

BOOST_AUTO_TEST_CASE(SpecialQRTest)
{
  SpecialQR qr(10);

  VectorXd e = VectorXd::Random(8);
  MatrixXd J1 = buildJj(e);
  MatrixXd J1back = J1;

  GivensSequence Q;
  bool b = qr.compute(J1, Q, false);

  BOOST_CHECK(isBandMatrix(J1, 0, 2));
  BOOST_CHECK(!b);
  BOOST_CHECK((Q.matrix(9)*J1).isApprox(J1back, 1e-15));


  //check reentry
  e = VectorXd::Random(9);
  MatrixXd J2 = buildJj(e);
  MatrixXd J2back = J2;
  b = qr.compute(J2, Q, false);
  BOOST_CHECK(isBandMatrix(J2, 0, 2));
  BOOST_CHECK(!b);
  BOOST_CHECK((Q.matrix(10)*J2).isApprox(J2back, 1e-15));

  //variant case
  MatrixXd J3 = buildJpm1(e);
  MatrixXd J3back = J3;
  b = qr.compute(J3, Q, true);
  BOOST_CHECK(isBandMatrix(J3, 0, 2));
  BOOST_CHECK(b);
  BOOST_CHECK((Q.matrix(9)*J3).isApprox(J3back, 1e-15));
}

BOOST_AUTO_TEST_CASE(SpecialQRTestExtended)
{
  int n = 8;
  VectorXd d = VectorXd::Random(n);
  VectorXd e = d.cwiseInverse();

  SpecialQR qr(n+1);
  LeastSquareObjective ls(d);

  MatrixXd J1(n + 1, n);
  MatrixXd R1(n + 1, n);
  GivensSequence Q1;
  Q1.reserve(n);
  ls.buildJj(J1, 0, n-1, StartType::Case3, EndType::Case1);
  bool b1 = qr.compute(R1, e, Q1, EndType::Case1);
  BOOST_CHECK(isBandMatrix(R1, 0, 2));
  BOOST_CHECK(b1);
  BOOST_CHECK((Q1.matrix(n+1)*R1).isApprox(J1, 1e-15));

  MatrixXd J2(n, n);
  MatrixXd R2(n, n);
  GivensSequence Q2;
  Q2.reserve(n-1);
  ls.buildJj(J2, 0, n - 1, StartType::Case3, EndType::Case2);
  bool b2 = qr.compute(R2, e, Q2, EndType::Case2);
  BOOST_CHECK(isBandMatrix(R2, 0, 2));
  BOOST_CHECK(b2);
  BOOST_CHECK((Q2.matrix(n)*R2).isApprox(J2, 1e-15));

  MatrixXd J3(n, n + 1);
  MatrixXd R3(n, n + 1);
  GivensSequence Q3;
  Q3.reserve(n-1);
  ls.buildJj(J3, 0, n - 1, StartType::Case3, EndType::Case3);
  bool b3 = qr.compute(R3, e, Q3, EndType::Case3);
  BOOST_CHECK(isBandMatrix(R3, 0, 2));
  BOOST_CHECK(b3);
  BOOST_CHECK((Q3.matrix(n)*R3).isApprox(J3, 1e-15));

  MatrixXd J4(n + 1, n + 1);
  MatrixXd R4(n + 1, n + 1);
  GivensSequence Q4;
  Q4.reserve(n);
  ls.buildJj(J4, 0, n - 1, StartType::Case3, EndType::Case4);
  bool b4 = qr.compute(R4, e, Q4, EndType::Case4);
  BOOST_CHECK(isBandMatrix(R4, 0, 2));
  BOOST_CHECK(!b4);
  BOOST_CHECK((Q4.matrix(n + 1)*R4).isApprox(J4, 1e-15));
}

#include <Givens.h>
#include <GivensSequence.h>
#include <ProblemMatrices.h>
#include <QRAlgorithms.h>

#include <iostream>

// boost
#define BOOST_TEST_MODULE MatrixComputationsTest
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace bms;

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

  bool b = hessenbergQR(M, seq);
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

  bool b = tridiagonalQR(M, seq);
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
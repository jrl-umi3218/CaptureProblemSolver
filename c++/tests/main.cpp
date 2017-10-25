#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <Eigen/QR>

#include <iostream>
#include <chrono>

#include "LeastSquare.h"
#include "Problem.h"
#include "ProblemMatrices.h"
#include "QRAlgorithms.h"
#include "QuadraticObjective.h"
#include "SQP.h"

using namespace Eigen;
using namespace bms;

void exampleMatrices(int n)
{
  VectorXd d = VectorXd::Random(n).cwiseAbs();

  std::cout << "J=\n" << buildJ(d) << std::endl;
  std::cout << "J0=\n" << buildJ0(d) << std::endl;
  std::cout << "Jj=\n" << buildJj(d) << std::endl;
  std::cout << "J{p-1}=\n" << buildJpm1(d) << std::endl;
  std::cout << "Cz=\n" << buildCZ(n) << std::endl;
  std::cout << "C=\n" << buildC(n) << std::endl;
}

/** n: size of matrices, N: number of tests*/
void QRPerformances(int n, const int N)
{
  VectorXd d = VectorXd::Random(n).cwiseAbs();
  MatrixXd J = buildJ(d);

  //dummy accumulator
  double acc = 0;

  std::cout << "On matrix J" << std::endl;
  //Eigen QR
  {
    //preallocations
    HouseholderQR<MatrixXd> qr(n - 1, n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      qr.compute(J);
      acc += qr.matrixQR()(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - Eigen::HouseholderQR: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count())/N << " microseconds" << std::endl;
  }

  //hessenbergQR
  {
    //preallocation
    MatrixXd Jcopy(n - 1, n);
    GivensSequence Q;
    Q.reserve(n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = J;
      hessenbergQR(Jcopy, Q);
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - bms::hessenbergQR " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  //tridiagonalQR
  {
    //preallocation
    MatrixXd Jcopy(n - 1, n);
    GivensSequence Q;
    Q.reserve(n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = J;
      tridiagonalQR(Jcopy, Q);
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - bms::tridiagonalQR " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  //copy overhead
  {
    //preallocation
    MatrixXd Jcopy(n - 1, n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = J;
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - copy overhead " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  std::cout << "\nOn matrix Jj" << std::endl;
  MatrixXd Jj = buildJj(d.head(n-1));
  //tridiagonalQR
  {
    //preallocation
    MatrixXd Jcopy(n, n);
    GivensSequence Q;
    Q.reserve(n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = Jj;
      tridiagonalQR(Jcopy, Q);
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - bms::tridiagonalQR " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  //specialQR
  {
    //preallocation
    MatrixXd Jcopy(n, n);
    SpecialQR qr(n);
    GivensSequence Q;
    Q.reserve(n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = Jj;
      qr.compute(Jcopy, Q,false);
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - bms::SpecialQR " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  //copy overhead
  {
    //preallocation
    MatrixXd Jcopy(n, n);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = Jj;
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - copy overhead " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }


  std::cout << "\ndummy: " << acc << std::endl;
  
}

void LSPerformance(int n, const int N)
{
  std::cout << "LSPerformance, size = " << n << std::endl;

  double d = 0;

  VectorXd l = -VectorXd::Random(n).cwiseAbs();
  VectorXd u = VectorXd::Random(n).cwiseAbs();
  LinearConstraints lc(l, u, -1, 1);

  VectorXd j = 100*VectorXd::Random(n);
  double c = -200;

  VectorXd delta = VectorXd::LinSpaced(n, 0.01, 0.02*n-0.01);
  LeastSquareObjective obj(delta);
  VectorXd Jx0 = VectorXd::Zero(n - 1);

  LeastSquare ls(n);
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      lc.resetActivation();
      ls.solveFeasibility(j,c,lc);
      d += ls.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << "LSFeasibility: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      lc.resetActivation();
      ls.solve(obj, Jx0, j, c, lc);
      d += ls.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << "LS: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }
  std::cout << d << std::endl;
}

void QRJAPerformance(int n, const int N)
{
  std::cout << "QRJAPreformance, size = " << n << std::endl;
  VectorXd d = VectorXd::Random(n).cwiseAbs();
  LeastSquareObjective obj(d);
  LinearConstraints lc(n);
  for (size_t i = 0; i < static_cast<size_t>(n); ++i)
  {
    if (rand() > RAND_MAX / 2)
      lc.activate(i, Activation::Lower);
  }
  auto na = lc.numberOfActiveConstraints();
  MatrixXd J = obj.matrix();
  MatrixXd JA(n - 1, n - na);
  lc.applyNullSpaceOnTheRight(JA, J);

  double acc = 0;
  {
    HouseholderQR<MatrixXd> qr(n - 1, n - na);
    MatrixXd Jcopy(n - 1, n - na);
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = JA;
      qr.compute(Jcopy);
      acc += qr.matrixQR()(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " HouseholderQR " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  {
    CondensedOrthogonalMatrix Q(n, n, 2 * n);
    MatrixXd Jcopy(n - 1, n - na);
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = JA;
      obj.qr(Jcopy, Q, na, lc.activeSet());
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " dedicated QR " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  {
    //preallocation
    MatrixXd Jcopy(n - 1, n - na);

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      Jcopy = JA;
      acc += Jcopy(0, 0);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - copy overhead " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  std::cout << acc << std::endl;
}

void readTest(const std::string& filepath)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);
}

void testSQPFeas(const std::string& filepath)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);
  
  Problem pb(raw);
  SQP sqp(static_cast<int>(raw.delta.size()));

  sqp.solveFeasibility(pb);
}

void testSQP(const std::string& filepath)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);

  Problem pb(raw);
  if (raw.delta.size() <= 15)
    pb.objective().precompute(1);
  SQP sqp(static_cast<int>(raw.delta.size()));

  sqp.solve(pb);
  std::cout << "  x = " << sqp.x().transpose() << std::endl;
  std::cout << "Phi = " << raw.Phi_.tail(raw.Phi_.size()-1).transpose() << std::endl;
  
  double fx;
  double fphi;
  pb.nonLinearConstraint().compute(fx, sqp.x());
  pb.nonLinearConstraint().compute(fphi, raw.Phi_.tail(raw.Phi_.size() - 1));
  std::cout << "  f(x) = " << fx << std::endl;
  std::cout << "f(Phi) = " << fphi << std::endl;
}

void SQPPerformance(const std::string& filepath, const int N)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);

  Problem pb(raw);
  if (raw.delta.size() <= 15)
    pb.objective().precompute(1);
  SQP sqp(static_cast<int>(raw.delta.size()));

  double d = 0;
  //{
  //  auto start_time = std::chrono::high_resolution_clock::now();
  //  for (int i = 0; i < N; ++i)
  //  {
  //    sqp.solveFeasibility(pb);
  //    d += sqp.x()[0];
  //  }
  //  auto end_time = std::chrono::high_resolution_clock::now();
  //  auto time = end_time - start_time;
  //  std::cout << "SQPFeasibility: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  //}
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      sqp.solve(pb);
      d += sqp.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << "SQP: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }
  std::cout << d << std::endl;
}

int main()
{
  //exampleMatrices(10);
  //QRPerformances(10, 10000);
  //QRPerformances(20, 10000);
  //QRPerformances(100, 1000);
  //LSPerformance(10, 1000);
  //LSPerformance(10, 10000);
  //LSPerformance(100, 1000);
  //LSPerformance(200, 1000);
  //LSPerformance(500, 1000);

  //readTest("data/Problem01.txt");

  //testSQP("data/Problem01.txt"); 
  SQPPerformance("data/Problem01.txt", 25000);
  SQPPerformance("data/Problem02.txt", 25000);
  SQPPerformance("data/Problem03.txt", 25000);
  //QRJAPerformance(10, 10000);
  //QRJAPerformance(20, 10000);
  //QRJAPerformance(50, 10000);
  //QRJAPerformance(100, 10000);
  //QRJAPerformance(200, 1000);
  //QRJAPerformance(500, 1000);

#ifdef WIN32
  system("pause");
#endif
}
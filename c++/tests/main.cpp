#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <Eigen/QR>

#include <fstream>
#include <iostream>
#include <chrono>

#include "LeastSquare.h"
#include "Problem.h"
#include "ProblemMatrices.h"
#include "QRAlgorithms.h"
#include "QuadraticObjective.h"
#include "SQP.h"
#include "SQPTestCommon.h"
#include "toMatlab.h"

using namespace Eigen;
using namespace bms;

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

void SQPPerformance(const std::string& filepath, int n, const int N)
{
  const int precompMax = 20;

  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);
  if (n > 0)
    raw = resampleProblem(raw, n);

  Problem pb(raw);
  Problem pb2(raw);
  if (raw.delta.size() <= precompMax)
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    pb.objective().precompute(1);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - precomputation: " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(time).count()) << " ms" << std::endl;
  }

  SQP sqp(static_cast<int>(raw.delta.size()));

  std::cout << "test SQP for n = " << raw.delta.size() << std::endl;
  double d = 0;
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      sqp.solveFeasibility(pb);
      d += sqp.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - SQPFeasibility: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      sqp.solve(pb);
      d += sqp.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - SQP: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  if (raw.delta.size() <= precompMax)
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      sqp.solve(pb2);
      d += sqp.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - SQP, no precomp: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }

  {
    std::vector<Activation> act = pb.linearConstraints().activationStatus();
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      pb.linearConstraints().setActivationStatus(act);
      sqp.solve(pb);
      act = pb.linearConstraints().activationStatus();
      d += sqp.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - SQP, warm start: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }
  std::cout << d << std::endl;
}

void mapFeasibleInputs()
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/data/Problem01.txt");

  Problem pb(raw);
  SQP sqp(static_cast<int>(raw.delta.size()));
  std::vector<Vector3d> points;
  points.reserve(500000);

  for (double zf = 0.0; zf <= 2; zf += 0.04)
  {
    for (double zi = 0.0; zi <= 2; zi += 0.04)
    {
      std::cout << zi << std::endl;
      for (double dzi = -5; dzi <= 5; dzi += 0.2)
      {
        pb.set_zi(zi);
        pb.set_dzi(dzi);
        pb.set_zf(zf);
        auto s = sqp.solveFeasibility(pb);
        if (s == SolverStatus::Converge)
        {
          double v;
          pb.nonLinearConstraint().compute(v, sqp.x());
          if (std::abs(v) < 1e-5)
            points.push_back({ zi,dzi,zf });
        }
      }
    }
  }

  Eigen::MatrixXd P(3, static_cast<DenseIndex>(points.size()));
  for (DenseIndex i = 0; i < P.cols(); ++i)
    P.col(i) = points[static_cast<size_t>(i)];

  std::ofstream aof("points.m");
  aof << "P = " << (toMatlab)P << ";" << std::endl;
}

int main()
{
  //QRPerformances(10, 10000);
  //QRPerformances(20, 10000);
  //QRPerformances(100, 1000);
  //LSPerformance(10, 1000);
  //LSPerformance(10, 10000);
  //LSPerformance(100, 1000);
  //LSPerformance(200, 1000);
  //LSPerformance(500, 1000);

  SQPPerformance("data/Problem01.txt", -1, 5000);
  //SQPPerformance("data/Problem01.txt", 20, 5000);
  //SQPPerformance("data/Problem01.txt", 50, 5000);
  //SQPPerformance("data/Problem01.txt", 100, 5000);
  //SQPPerformance("data/Problem01.txt", 200, 1000);
  //SQPPerformance("data/Problem01.txt", 500, 1000);
  //SQPPerformance("data/Problem01.txt", 1000, 1000);
  //QRJAPerformance(10, 10000);
  //QRJAPerformance(20, 10000);
  //QRJAPerformance(50, 10000);
  //QRJAPerformance(100, 10000);
  //QRJAPerformance(200, 1000);
  //QRJAPerformance(500, 1000);

  //mapFeasibleInputs();

#ifdef WIN32
  system("pause");
#endif
}

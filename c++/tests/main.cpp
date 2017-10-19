#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <Eigen/QR>

#include <iostream>
#include <chrono>

#include "LeastSquare.h"
#include "ProblemMatrices.h"
#include "QRAlgorithms.h"

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
  double d = 0;

  VectorXd l = -VectorXd::Random(n).cwiseAbs();
  VectorXd u = VectorXd::Random(n).cwiseAbs();
  LinearConstraints lc(l, u, -1, 1);

  VectorXd j = VectorXd::Random(n);
  double c = -10;

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
}

int main()
{
  //exampleMatrices(10);
  //QRPerformances(10, 10000);
  //QRPerformances(20, 10000);
  //QRPerformances(100, 1000);
  LSPerformance(10, 1000);
  LSPerformance(10, 10000);
  LSPerformance(100, 1000);
  LSPerformance(200, 1000);
  LSPerformance(500, 1000);

#ifdef WIN32
  system("pause");
#endif
}
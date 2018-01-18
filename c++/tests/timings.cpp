/* Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
 *
 * This file is part of CPS.
 *
 * CPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CPS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "timings.h"

using namespace Eigen;
using namespace cps;

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
    std::cout << " - Eigen::HouseholderQR: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
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
  MatrixXd Jj = buildJj(d.head(n - 1));
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
      qr.compute(Jcopy, Q, false);
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

  VectorXd j = 100 * VectorXd::Random(n);
  double c = -200;

  VectorXd delta = VectorXd::LinSpaced(n, 0.01, 0.02*n - 0.01);
  LeastSquareObjective obj(delta);
  VectorXd Jx0 = VectorXd::Zero(n - 1);

  LeastSquare ls(n);
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      lc.resetActivation();
      ls.solveFeasibility(j, c, lc);
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

#ifdef LSSOL_OPTION
  auto opt = ls.parameters();
  ls.parameters(opt.optimizedSolver(false));
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
    std::cout << "LS-lssol: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }
#endif

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

#ifdef LSSOL_OPTION
  {
    auto opt = sqp.LSParameters();
    sqp.LSParameters(opt.optimizedSolver(false));
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i)
    {
      sqp.solve(pb);
      d += sqp.x()[0];
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << " - SQP, lssol: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / N << " microseconds" << std::endl;
  }
#endif
  std::cout << d << std::endl;
}

void SQPTimings(const std::vector<std::string>& filepath, const std::vector<int>& n, const int N)
{
  const int precompMax = 20;

  std::vector<RawProblem> raw;
  std::string base = TESTS_DIR;
  for (auto s : filepath)
  {
    RawProblem r;
    r.read(base + "/" + s);
    raw.push_back(r);
  }

  //dummy
  double d = 0;

  for (auto ni : n)
  {
    std::vector<Problem> precompPb;
    std::vector<Problem> pb;

    bool precomp;
    if (ni > 0)
    {
      for (auto r : raw)
      {
        r = resampleProblem(r, ni);
        pb.emplace_back(r);
        if (r.delta.size() <= precompMax)
          precompPb.emplace_back(r);
      }
    }
    else
    {
      ni = static_cast<int>(raw[0].delta.size());
      for (auto r : raw)
      {
        if (ni != static_cast<int>(r.delta.size()))
        {
          throw std::runtime_error("All problem must have the same size. Consider resampling if it is not the case.");
        }
        pb.emplace_back(r);
        if (r.delta.size() <= precompMax)
          precompPb.emplace_back(r);
      }
    }

    std::cout << "Problem size: " << ni << std::endl;
    if (ni <= precompMax)
    {
      auto start_time = std::chrono::high_resolution_clock::now();
      for (auto& pbi: precompPb)
        pbi.objective().precompute(1);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << " - precomputation: " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(time).count())/precompPb.size() << " ms" << std::endl;
    }

    SQP sqp(ni);
    {
      auto start_time = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i)
      {
        for (const auto& pbi : pb)
        {
          sqp.solve(pbi);
          d += sqp.x()[0];
        }
      }
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << " - SQP: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / (N*pb.size()) << " microseconds" << std::endl;
    }

    if (ni <= precompMax)
    {
      auto start_time = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i)
      {
        for (const auto& pbi : precompPb)
        {
          sqp.solve(pbi);
          d += sqp.x()[0];
        }
      }
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << " - SQP, precomp: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / (N*precompPb.size()) << " microseconds" << std::endl;
    }

#ifdef LSSOL_OPTION
    {
      auto opt = sqp.LSParameters();
      sqp.LSParameters(opt.optimizedSolver(false));
      auto start_time = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i)
      {
        for (const auto& pbi : pb)
        {
          sqp.solve(pbi);
          d += sqp.x()[0];
        }
      }
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << " - SQP, lssol: " << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time).count()) / (N*pb.size()) << " microseconds" << std::endl;
    }
#endif
  }

  std::cout << d;
}


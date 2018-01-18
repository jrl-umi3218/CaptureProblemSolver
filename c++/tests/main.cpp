#include <cps/Problem.h>
#include <cps/SQP.h>
#include <cps/toMatlab.h>

#include "SQPTestCommon.h"
#include "timings.h"

using namespace Eigen;
using namespace cps;


void mapFeasibleInputs()
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/data/Problem01.txt");

  Problem pb(raw);
  SQP sqp(static_cast<int>(raw.delta.size()));
  std::vector<Vector3d> points;
  points.reserve(500000);

  for (double target_height = 0.0; target_height <= 2; target_height += 0.04)
  {
    for (double init_zbar = 0.0; init_zbar <= 2; init_zbar += 0.04)
    {
      std::cout << init_zbar << std::endl;
      for (double init_zbar_deriv = -5; init_zbar_deriv <= 5; init_zbar_deriv += 0.2)
      {
        pb.set_init_zbar(init_zbar);
        pb.set_init_zbar_deriv(init_zbar_deriv);
        pb.set_target_height(target_height);
        auto s = sqp.solveFeasibility(pb);
        if (s == SolverStatus::Converge)
        {
          double v;
          pb.nonLinearConstraint().compute(v, sqp.x());
          if (std::abs(v) < 1e-5)
            points.push_back({init_zbar, init_zbar_deriv, target_height});
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

void SQPSolveTest(const std::string& filepath)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);

  Problem pb(raw);
  pb.objective().precompute(1);
  SQP sqp(static_cast<int>(raw.delta.size()));
  auto s = sqp.solve(pb);
  std::cout << static_cast<int>(s) << std::endl;
  std::cout << sqp.x().transpose() << std::endl;
  std::cout << raw.Phi_.transpose() << std::endl;

  auto statistics = sqp.statistics();
  std::cout << "(iter, act, deact, #actCstr)" << std::endl;
  for (size_t i = 0; i < statistics.lsStats.size(); ++i)
  {
    auto si = statistics.lsStats[i];
    std::cout << si.iter << ", " << si.activation << ", " << si.deactivation << ", " << si.activeConstraints << std::endl;
  }
}

int main()
{
  //QRPerformances(10, 10000);
  //QRPerformances(20, 10000);
  //QRPerformances(100, 1000);
  //LSPerformance(10, 5000);
  //LSPerformance(20, 5000);
  //LSPerformance(100, 1000);
  //LSPerformance(200, 1000);
  //LSPerformance(500, 1000);

  //SQPPerformance("data/Problem01.txt", -1, 5000);
  //SQPPerformance("data/Problem02.txt", -1, 5000);
  //SQPPerformance("data/Problem03.txt", -1, 5000);
  //SQPPerformance("data/Problem01.txt", 20, 5000);
  //SQPPerformance("data/Problem02.txt", 20, 5000);
  //SQPPerformance("data/Problem03.txt", 20, 5000);
  //SQPPerformance("data/Problem01.txt", 50, 5000);
  //SQPPerformance("data/Problem02.txt", 50, 5000);
  //SQPPerformance("data/Problem03.txt", 50, 5000);
  //SQPPerformance("data/Problem01.txt", 100, 5000);
  //SQPPerformance("data/Problem02.txt", 100, 5000);
  //SQPPerformance("data/Problem03.txt", 100, 5000);
  //SQPPerformance("data/Problem01.txt", 200, 1000);
  //SQPPerformance("data/Problem02.txt", 200, 1000);
  //SQPPerformance("data/Problem03.txt", 200, 1000);
  //SQPPerformance("data/Problem01.txt", 500, 1000);
  //SQPPerformance("data/Problem01.txt", 1000, 1000);
  //QRJAPerformance(10, 10000);
  //QRJAPerformance(20, 10000);
  //QRJAPerformance(50, 10000);
  //QRJAPerformance(100, 10000);
  //QRJAPerformance(200, 1000);
  //QRJAPerformance(500, 1000);

  //mapFeasibleInputs();
  //SQPSolveTest("data/Problem02.txt");

  //SQPTimings({ "data/Problem01.txt", "data/Problem02.txt" , "data/Problem03.txt" , "data/Problem04.txt" }, { -1, 15, 20, 50, 100 }, 1000);
  SQPTimings({ "data/Problem01.txt" }, { -1, 15, 20, 50, 100 }, 1000);

#ifdef WIN32
  system("pause");
#endif
}

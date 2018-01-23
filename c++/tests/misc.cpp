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

#include "misc.h"

#include <cps/Problem.h>
#include <cps/SQP.h>
#include <cps/toMatlab.h>

#include <fstream>
#include <vector>

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
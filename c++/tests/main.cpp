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

#include <cps/Problem.h>
#include <cps/SQP.h>

#include "misc.h"
#include "SQPTestCommon.h"
#include "timings.h"

using namespace Eigen;
using namespace cps;

/** Solve the problem specified by filepath (relatively to TEST_DIR) and display some statistics.*/
void SQPSolveTest(const std::string& filepath)
{
  //reading a problem from file
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/" + filepath);

  //creating Problem and SQP objects
  Problem pb(raw);
  SQP sqp(static_cast<int>(raw.delta.size()));

  //solving and displaying
  std::cout << "Solving problem " << filepath << "\n" <<std::endl;
  auto s = sqp.solve(pb);
  std::cout << "  solver status: " << static_cast<int>(s) << std::endl;
  std::cout << "  sqp solution:  " << sqp.x().transpose() << std::endl;
  std::cout << "  raw solution:  " << raw.Phi_.tail(raw.Phi_.size()-1).transpose() << "\n" << std::endl;

#ifdef USE_STATS
  auto statistics = sqp.statistics();
  std::cout << " (iter,\tact,\tdeact,\t#actCstr)" << std::endl;
  for (size_t i = 0; i < statistics.lsStats.size(); ++i)
  {
    auto si = statistics.lsStats[i];
    std::cout << "  " << si.iter << ",\t" << si.activation << ",\t" << si.deactivation << ",\t" << si.activeConstraints << std::endl;
  }
  std::cout << std::endl;
#endif //USE_STATS
}

// Additional arguments are paths of problem files, relative to TESTS_DIR
int main(int argc, char* argv[])
{
  if (argc > 1)
  {
    for (int i = 1; i < argc; ++i)
    {
      SQPSolveTest(argv[i]);
    }
  }
  else
  {
    SQPSolveTest("data/Problem01.txt");
  }

#ifdef WIN32
  system("pause");
#endif
}

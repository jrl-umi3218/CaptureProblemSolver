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

#include <iostream>

#include <cps/SQP.h>

#include "SQPTestCommon.h"

// boost
#define BOOST_TEST_MODULE SQPTests
#include <boost/test/unit_test.hpp>

using namespace Eigen;
using namespace cps;

BOOST_AUTO_TEST_CASE(testFeasibilitySQP)
{
  RawProblem raw;
  std::string base = TESTS_DIR;
  raw.read(base + "/data/Problem01.txt");

  Problem pb(raw);
  SQP sqp(static_cast<int>(raw.delta.size()));

  auto s = sqp.solveFeasibility(pb);

  double c;
  pb.nonLinearConstraint().compute(c, sqp.x());
  BOOST_CHECK(s == SolverStatus::Converge);
  BOOST_CHECK(pb.linearConstraints().checkPrimal(sqp.x()));
  BOOST_CHECK(std::abs(c) <= 1e-7);
}


BOOST_AUTO_TEST_CASE(BasicSQPTest)
{
  {
    RawProblem raw;
    std::string base = TESTS_DIR;
    raw.read(base + "/data/Problem01.txt");

    Problem pb(raw);
    if (raw.delta.size() <= 15)
      pb.objective().precompute(1);
    SQP sqp(static_cast<int>(raw.delta.size()));

    auto s = sqp.solve(pb);
    BOOST_CHECK(s == SolverStatus::Converge);
    BOOST_CHECK(sqp.x().isApprox(raw.Phi_.tail(raw.Phi_.size()-1), 1e-8));
  }
  {
    RawProblem raw;
    std::string base = TESTS_DIR;
    raw.read(base + "/data/Problem02.txt");

    Problem pb(raw);
    if (raw.delta.size() <= 15)
      pb.objective().precompute(1);
    SQP sqp(static_cast<int>(raw.delta.size()));

    auto s = sqp.solve(pb);
    BOOST_CHECK(s == SolverStatus::Converge);
    BOOST_CHECK(sqp.x().isApprox(raw.Phi_.tail(raw.Phi_.size() - 1), 1e-6));
  }
  {
    RawProblem raw;
    std::string base = TESTS_DIR;
    raw.read(base + "/data/Problem03.txt");

    Problem pb(raw);
    if (raw.delta.size() <= 15)
      pb.objective().precompute(1);
    SQP sqp(static_cast<int>(raw.delta.size()));

    auto s = sqp.solve(pb);
    BOOST_CHECK(s == SolverStatus::Converge);
    BOOST_CHECK(sqp.x().isApprox(raw.Phi_.tail(raw.Phi_.size() - 1), 1e-7));
  }
}

BOOST_AUTO_TEST_CASE(ResampledProblems)
{
  auto raw = resampleProblem("data/Problem01.txt", 20);
  SQP sqp(static_cast<int>(raw.delta.size()));

  auto s = sqp.solve(raw);

  BOOST_CHECK(s == SolverStatus::Converge);
}

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

#include "SQPTestCommon.h"

namespace cps
{
  RawProblem resampleProblem(const RawProblem& raw, int n)
  {
    auto raw2 = raw;
    VectorXd s = VectorXd::LinSpaced(n + 1, 0, 1);
    raw2.delta.resize(n);
    raw2.delta.array() = s.tail(n).array()*s.tail(n).array() - s.head(n).array()*s.head(n).array();

    //we invalidate the solution
    raw2.Phi_.resize(0);

    return raw2;
  }

  RawProblem resampleProblem(const std::string& filepath, int n)
  {
    RawProblem raw;
    std::string base = TESTS_DIR;
    raw.read(base + "/" + filepath);

    return resampleProblem(raw, n);
  }
}
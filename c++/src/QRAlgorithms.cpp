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

#include <cps/QRAlgorithms.h>

using namespace Eigen;

namespace cps
{
  SpecialQR::SpecialQR(DenseIndex kmax)
    : e_(kmax + 1)
    , c1_(kmax)
    , c2_(kmax)
    , d_(kmax)
    , l_(kmax)
    , dprecomp_(-VectorXd::LinSpaced(static_cast<int>(kmax),1, static_cast<int>(kmax)))
    , lprecomp_((ArrayXd::LinSpaced(static_cast<int>(kmax), 1, static_cast<int>(kmax)).sqrt()*ArrayXd::LinSpaced(static_cast<int>(kmax), 2, static_cast<int>(kmax) + 1).sqrt()).inverse().matrix())
  {
    G_.reserve(kmax);
    for (DenseIndex i = 0; i < kmax; ++i)
    {
      double c = 1 / sqrt(double(i) + 2.);
      G_.emplace_back(i, c, sqrt(double(i) + 1.)*c);
    }
  }
}

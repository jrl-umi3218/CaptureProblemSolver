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

#include <cps/Givens.h>

namespace cps
{
  Givens::Givens()
    : Givens(0,0,1,0)
  {
  }

  Givens::Givens(Index i, Index j, double c, double s)
    : i_(i), j_(j), Jt_(c,-s)
  {
  }

  Givens::Givens(Index i, double c, double s)
    : Givens(i,i+1,c,s)
  {
  }

  void Givens::extend(Index incr)
  {
    i_ += incr;
    j_ += incr;
  }
}

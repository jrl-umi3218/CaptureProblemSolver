#pragma once
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

/** Compute a discrete view of the feasible space in (z_i, dz_i, z_f), i.e. the
  * set of inputs in this space for which the problem is feasible.
  * This is output in points.m.
  */
void mapFeasibleInputs();
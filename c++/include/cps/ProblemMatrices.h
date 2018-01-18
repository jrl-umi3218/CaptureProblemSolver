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

#include <Eigen/Core>

#include <cps/cps_api.h>

namespace cps
{
  /** A collection of functions to generate the problem matrices.
    *
    * Mostly for debugging purposes
    */

  /** Build J such that the objective is ||Jx||^2 where d = 1/delta */
  CPS_DLLAPI Eigen::MatrixXd buildJ(const Eigen::VectorXd& d);

  /** Build J_0 without the initial 0 rows (see doc)*/
  CPS_DLLAPI Eigen::MatrixXd buildJ0(const Eigen::VectorXd& e);
  
  /** Build J_j (see doc)*/
  CPS_DLLAPI Eigen::MatrixXd buildJj(const Eigen::VectorXd& e);
  
  /** Build J_{p-1} (see doc)*/
  CPS_DLLAPI Eigen::MatrixXd buildJpm1(const Eigen::VectorXd& e);

  /** Build the matrix of zonotopique constraints l_i <= x_i - x_{i-1} <= u_i
    * n is the size of x
    */
  CPS_DLLAPI Eigen::MatrixXd buildCZ(int n);

  /** [C_Z;e_n^T] where e_n is the n-th column of the nxn identity matrix.*/
  CPS_DLLAPI Eigen::MatrixXd buildC(int n);
}

#pragma once

#include <Eigen/Core>

#include <bms_api.h>

namespace bms
{
  /** A collection of functions to generate the problem matrices.
    *
    * Mostly for debugging purposes
    */

  /** Build J such that the objective is ||Jx||^2 where d = 1/sqrt(delta) */
  BMS_DLLAPI Eigen::MatrixXd buildJ(const Eigen::VectorXd& d);

  /** Build J_0 without the initial 0 rows (see doc)*/
  BMS_DLLAPI Eigen::MatrixXd buildJ0(const Eigen::VectorXd& e);
  
  /** Build J_j (see doc)*/
  BMS_DLLAPI Eigen::MatrixXd buildJj(const Eigen::VectorXd& e);
  
  /** Build J_{p-1} (see doc)*/
  BMS_DLLAPI Eigen::MatrixXd buildJpm1(const Eigen::VectorXd& e);

  /** Build the matrix of zonotopique constraints l_i <= x_i - x_{i-1} <= u_i
    * n is the size of x
    */
  BMS_DLLAPI Eigen::MatrixXd buildCZ(int n);

  /** [C_Z;e_n^T] where e_n is the n-th column of the nxn identity matrix.*/
  BMS_DLLAPI Eigen::MatrixXd buildC(int n);
}
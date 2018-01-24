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

#include <cps/SQP.h>

namespace cps
{
  /** A solver for the capture problem.
    *
    * The class is a wrapper around SQP to handle the differences between the
    * original capture problem and the penalty-based relaxation solved by the
    * SQP. It offers two main functions:
    *  - solve
    *  - solveFeasibility
    */
  class CPS_DLLAPI CaptureSolver
  {
  public:
    /** Build a solver for problems of size n, with an accepted violation on 
      * the non-linear constraint.
      */
    CaptureSolver(int n, double violation = 1e-6);
    /** Set the accepted violation on the non-linear constraint.*/
    void violation(double v);
    /** Get the current accepted violation on the non-linear constraint.*/
    double violation() const;
    /** Set a new set of parameters for the SQP.*/
    void SQPParameters(const SQP::Parameters& p);
    /** Set a new set of parameters for the SQP's least-square solver.*/
    void LSParameters(const LeastSquare::Parameters& p);
    /** Get the current parameters of the SQP.*/
    const SQP::Parameters& SQPParameters() const;
    /** Get the current parameters of the SQP's least-square solver.*/
    const LeastSquare::Parameters& LSParameters() const;

    /** Solve the following problem: 
      *
      * min. 1/2 ||J x||^2
      * s.t. b(x) = 0
      *      l_i <= x_i-x_{i-1} <= u_i for i=0..n-1  (x_{-1} = 0)
      *      xln <= x_{n-1} <= xun
      *
      * where b(x) is the boundedness constraint violation described by pb.nonLinearConstraint(),
      * J is described by pb.objective() and the linear constraints by pb.linearConstraints().
      */
    SolverStatus solve(const Problem& pb);
    /** Solve the following problem: 
      *
      * min. 1/2 ||b(x)||^2
      * s.t. l_i <= x_i-x_{i-1} <= u_i for i=0..n-1  (x_{-1} = 0)
      *      xln <= x_{n-1} <= xun
      *
      * where b(x) is the boundedness constraint violation described by pb.nonLinearConstraint(),
      * and the constraints are described by pb.linearConstraints().
      */
    SolverStatus solveFeasibility(const Problem& pb);

    /** Retrieve the solution (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& x() const;
    /** Retrieve the Lagrange multipliers (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& lambda() const;
    /** Retrieve the current active set.*/
    const std::vector<Activation>& activeSet() const;
    /** Get the number of SQP iterations of the last run.*/
    int numberOfIterations() const;

    /** Get stats on the last run. Only meaningful if USE_STATS is defined*/
    const stats::SQPStats& statistics() const;

  private:
    double violation_;
    SQP sqp_;
    
  };
}
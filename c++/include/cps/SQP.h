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
#include <cps/LeastSquare.h>
#include <cps/Problem.h>
#include <cps/QuadraticObjective.h>
#include <cps/Statistics.h>

namespace cps
{
  /** A solver for the (modified) capture problem.
    *
    * The class in itself handle the memory and statistics useful for its two
    * main methods:
    *  - solve
    *  - solveFeasibility
    */
  class CPS_DLLAPI SQP
  {
  public:
    /** A set of parameters for controlling the execution of the resolution.*/
    class Parameters
    {
    public:
      int maxIter()           const { return maxIter_; }
      double mu()             const { return mu_; }
      double beta()           const { return beta_; }
      double c1()             const { return c1_; }
      double smallestLSStep() const { return smallestLSStep_; }
      double tau_p()          const { return tau_p_; }
      double tau_d()          const { return tau_d_; }
      double feasibilityEps() const { return feasibilityEps_; }

      Parameters& maxIter       (int i)    { maxIter_ = i;        return *this; }
      Parameters& mu            (double d) { mu_ = d;             return *this; }
      Parameters& beta          (double d) { beta_ = d;           return *this; }
      Parameters& c1            (double d) { c1_ = d;             return *this; }
      Parameters& smallestLSStep(double d) { smallestLSStep_ = d; return *this; }
      Parameters& tau_p         (double d) { tau_p_ = d;          return *this; }
      Parameters& tau_d         (double d) { tau_d_ = d;          return *this; }
      Parameters& feasibilityEps(double d) { feasibilityEps_ = d; return *this; }

    private:
      int maxIter_ = 100;            //maximum number of iterations
      double mu_ = 100000;           //penalty parameter
      double beta_ = 0.9;            //backtracking multiplier in line search
      double c1_ = 0.01;             //gradient coefficient in line search
      double smallestLSStep_ = 1e-8;
      double tau_p_ = 1e-6;          //precision parameter on primal condition
      double tau_d_ = 1e-6;          //precision parameter on dual condition
      double feasibilityEps_ = 1e-8; //margin of feasibility for which a potential numerical issue is detected
    };

    /** Build a solver for problems of size n.*/
    SQP(int n);
    /** Set a new set of parameters for the SQP.*/
    void SQPParameters(const Parameters& p);
    /** Set a new set of parameters for the underlying least-square solver.*/
    void LSParameters(const LeastSquare::Parameters& p);
    /** Get the current parameters of the SQP.*/
    const Parameters& SQPParameters() const;
    /** Get the current parameters of the underlying least-square solver.*/
    const LeastSquare::Parameters& LSParameters() const;

    /***/
    SolverStatus solve(const Problem& pb);
    /***/
    SolverStatus solveFeasibility(const Problem& pb);

    /** Retrieve the solution (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& x() const;
    /** Retrieve the Lagrange multipliers (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& lambda() const;
    /** Retrieve the current active set.*/
    const std::vector<Activation>& activeSet() const;
    /** Get the number of iterations of the last run.*/
    int numberOfIterations() const;

    /** Get stats on the last run. Only meaningful if USE_STATS is defined*/
    const stats::SQPStats& statistics() const;

  private:
    bool checkKKT(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, 
                  const LinearConstraints& lc, const LeastSquareObjective* const obj = nullptr) const;

    Eigen::DenseIndex n_;
    LeastSquare ls_;

    //optimization parameters
    Parameters params_;

    //optimization data
    LinearConstraints shiftedLC_;                 //constraints passed to the LS
    Eigen::VectorXd x_;                           //value of the iterate
    Eigen::VectorXd xa_;                          //value of the current point in the line search
    Eigen::VectorXd j_;                           //gradient of the non-linear function
    Eigen::VectorXd lambda_;                      //Lagrange multipliers
    mutable Eigen::VectorXd Cx_;                  //value of C*x
    mutable Eigen::VectorXd Cl_;                  //value of C^T*lambda
    mutable Eigen::VectorXd Jx_;                  //value of J*x
    mutable Eigen::VectorXd JtJx_;                //value of J^T * J *x
    mutable Eigen::VectorXd Jp_;                  //value of J*p
    mutable Eigen::VectorXd gradL_;               //gradient of the Lagrangian
    std::vector<Activation> currentActiveSet_;    //current active set
    std::vector<Activation> previousActiveSet_;   //active set backup
    int k_;                                       //iteration number

    //stats
    stats::SQPStats stats_;
  };
}

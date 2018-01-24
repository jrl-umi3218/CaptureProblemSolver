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

#include <cps/CaptureSolver.h>

using namespace Eigen;

namespace cps
{
  CaptureSolver::CaptureSolver(int n, double viol)
    : sqp_(n)
    , violation_(viol)
  {
  }

  void CaptureSolver::violation(double v)
  {
    assert(v > 0);
    violation_ = v;
  }

  double CaptureSolver::violation() const
  {
    return violation_;
  }

  void CaptureSolver::SQPParameters(const SQP::Parameters& p)
  {
    sqp_.SQPParameters(p);
  }

  void CaptureSolver::LSParameters(const LeastSquare::Parameters& p)
  {
    sqp_.LSParameters(p);
  }

  const SQP::Parameters & CaptureSolver::SQPParameters() const
  {
    return sqp_.SQPParameters();
  }

  const LeastSquare::Parameters & CaptureSolver::LSParameters() const
  {
    return sqp_.LSParameters();
  }

  SolverStatus CaptureSolver::solve(const Problem& pb)
  {
    //check input
    if (pb.init_omega_min() > pb.init_omega_max())
    {
      return SolverStatus::NoLinearlyFeasiblePoint;
    }

    auto s = sqp_.solve(pb);
    
    //check if the problem is feasible for the non-linear constraint
    if (s == SolverStatus::Converge)
    {
      if (std::abs(pb.nonLinearConstraint().compute(sqp_.x())) <= violation_)
      {
        return SolverStatus::Converge;
      }
      else
      {
        return SolverStatus::Unfeasible;
      }
    }
    else
    {
      return s;
    }
  }

  SolverStatus CaptureSolver::solveFeasibility(const Problem& pb)
  {
    return sqp_.solveFeasibility(pb);;
  }

  const Eigen::VectorXd & CaptureSolver::x() const
  {
    return sqp_.x();
  }

  const Eigen::VectorXd & CaptureSolver::lambda() const
  {
    return sqp_.lambda();
  }

  const std::vector<Activation>& CaptureSolver::activeSet() const
  {
    return sqp_.activeSet();
  }

  int CaptureSolver::numberOfIterations() const
  {
    return sqp_.numberOfIterations();;
  }

  const stats::SQPStats & CaptureSolver::statistics() const
  {
    return sqp_.statistics();
  }
}
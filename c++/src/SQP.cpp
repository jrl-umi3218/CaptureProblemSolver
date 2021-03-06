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

using namespace Eigen;

namespace cps
{
  SQP::SQP(int n)
    : n_(n)
    , ls_(n)
    , shiftedLC_(n)
    , x_(n)
    , xa_(n)
    , j_(n)
    , lambda_(n + 1)
    , Cx_(n + 1)
    , Cl_(n)
    , Jx_(n - 1)
    , JtJx_(n)
    , Jp_(n - 1)
    , gradL_(n)
    , currentActiveSet_(n + 1, Activation::None)
    , previousActiveSet_(n + 1, Activation::None)
  {
  }

  void SQP::SQPParameters(const Parameters & p)
  {
    params_ = p;
  }

  void SQP::LSParameters(const LeastSquare::Parameters & p)
  {
    ls_.parameters(p);
  }

  const SQP::Parameters& SQP::SQPParameters() const
  {
    return params_;
  }

  const LeastSquare::Parameters & SQP::LSParameters() const
  {
    return ls_.parameters();
  }


  SolverStatus SQP::solve(const Problem& pb)
  {
    k_ = 0;  // reset iteration number
    STATISTICS(stats_.reset());

    //shortcuts
    const auto& lc = pb.linearConstraints();
    const auto& nlc = pb.nonLinearConstraint();
    const auto& obj = pb.objective();

    //active sets
    currentActiveSet_ = lc.activationStatus();
    previousActiveSet_ = lc.activationStatus();

    //penality parameter
    double mu = params_.mu();
    double mu2 = mu*mu;

    //feasible (x,lambda)
    FeasiblePointInfo fpi;
    std::tie(fpi, x_) = lc.initialPoint(true, params_.feasibilityEps());
    if (fpi != FeasiblePointInfo::Found) //FIXME: we can be more precise here
    {
      //attempt to find a point without the active set
      std::tie(fpi, x_) = lc.initialPoint(false, params_.feasibilityEps());
      if (fpi != FeasiblePointInfo::Found)
        return SolverStatus::NoLinearlyFeasiblePoint;
      else
      {
        //reset the active sets
        for (size_t i = 0; i < currentActiveSet_.size(); ++i)
        {
          if (lc.activationStatus(i) != Activation::Equal)
          {
            currentActiveSet_[i] = Activation::None;
            previousActiveSet_[i] = Activation::None;
          }
        }

      }
    }
    lambda_.setZero();

    //main loop
    for (k_ = 0; k_ < params_.maxIter(); ++k_)
    {
      double f;
      nlc.compute(f, j_, x_);
      f *= mu;
      j_ *= mu;
      obj.applyJToTheLeft(Jx_, x_);
      double v = 0.5* Jx_.squaredNorm();

      if (checkKKT(x_, lambda_, f, j_, lc, &obj))
        return SolverStatus::Converge;

      //solve the least square
      shiftedLC_ = lc.shift(x_, true);
      previousActiveSet_ = currentActiveSet_;
      shiftedLC_.setActivationStatus(currentActiveSet_);
      if (ls_.solve(obj, Jx_, j_, f, shiftedLC_) != SolverStatus::Converge)
      {
        STATISTICS(stats_.lsStats.push_back(ls_.statistics()));
        return SolverStatus::Fail;
      }
      STATISTICS(stats_.lsStats.push_back(ls_.statistics()));
      currentActiveSet_ = shiftedLC_.activationStatus();
      const auto& p = ls_.x();
      const auto& pl = ls_.lambda();

      //line search
      double a = 1;
      obj.applyJToTheLeft(Jx_, x_);
      obj.applyJToTheLeft(Jp_, p);
      double cgp = params_.c1()*(f*j_.dot(p) + Jx_.dot(Jp_));

      xa_ = x_ + a*p;
      double fa;
      nlc.compute(fa, xa_);
      double va = obj.value(xa_);
      double val0 = 0.5*f*f + v;
      STATISTICS(int kls = 0);
      while (0.5*mu2*fa*fa + va > val0 + a*cgp)
      {
        STATISTICS(++kls);
        a = params_.beta() * a;
        if (a < params_.smallestLSStep())
        {
          STATISTICS(stats_.lineSearchSteps.push_back(kls));
          return SolverStatus::LineSearchFailed;
        }
        xa_ = x_ + a*p;
        nlc.compute(fa, xa_);
        va = obj.value(xa_);
      }
      STATISTICS(stats_.lineSearchSteps.push_back(kls));
      if (x_.isApprox(xa_))
      {
        ++k_; //we did full iteration, we need to count it.
        //x_ didn't change so we just need to update lambda but not the objective-related values
        lambda_ += a*(pl - lambda_);
        if (checkKKT(x_, lambda_, f, j_, lc, &obj))
          return SolverStatus::Converge;
        else
          return SolverStatus::NumericallyEquivalentIterates;
      }
      x_ = xa_;
      lambda_ += a*(pl - lambda_);
      if (a < 1)
        currentActiveSet_ = previousActiveSet_;
    }

    return SolverStatus::MaxIteration;
  }



  SolverStatus SQP::solveFeasibility(const Problem& pb)
  {
    k_ = 0;  // reset iteration number

    //shortcuts
    const auto& lc = pb.linearConstraints();
    const auto& nlc = pb.nonLinearConstraint();

    //feasible (x,lambda)
    FeasiblePointInfo fpi;
    std::tie(fpi, x_) = lc.initialPoint(false, params_.feasibilityEps());
    if (fpi != FeasiblePointInfo::Found) //FIXME: we can be more precise here
      return SolverStatus::NoLinearlyFeasiblePoint; 
    lambda_.setZero();

    //main loop
    for (int k = 0; k < params_.maxIter(); ++k)
    {
      double f;
      nlc.compute(f, j_, x_);
      //std::cout << "||c||^2 = " << f*f << std::endl;

      if (checkKKT(x_, lambda_, f, j_, lc))
        return SolverStatus::Converge;

      //solve the least square
      shiftedLC_ = lc.shift(x_, true);
      ls_.solveFeasibility(j_, f, shiftedLC_);
      const auto& p = ls_.x();
      const auto& pl = ls_.lambda();

      //line search
      double a = 1;
      double cgp = params_.c1()*f*j_.dot(p);

      xa_ = x_ + a*p;
      double fa;
      nlc.compute(fa, xa_);
      while (0.5*fa*fa > 0.5*f*f + a*cgp)
      {
        a = params_.beta() * a;
        if (a < params_.smallestLSStep())
          return SolverStatus::LineSearchFailed;
        xa_ = x_ + a*p;
        nlc.compute(fa, xa_);
      }
      if (x_.isApprox(xa_))
      {
        //x_ didn't change so we just need to update lambda but not the objective-related values
        lambda_ += a*(pl - lambda_);
        if (checkKKT(x_, lambda_, f, j_, lc))
          return SolverStatus::Converge;
        else
          return SolverStatus::NumericallyEquivalentIterates;
      }
      x_ = xa_;
      lambda_ += a*(pl - lambda_);
    }

    return SolverStatus::MaxIteration;
  }

  const Eigen::VectorXd & SQP::x() const
  {
    return x_;
  }
  const Eigen::VectorXd & SQP::lambda() const
  {
    return lambda_;
  }

  const std::vector<Activation>& SQP::activeSet() const
  {
    return currentActiveSet_;
  }

  int SQP::numberOfIterations() const
  {
    return k_;
  }

  const stats::SQPStats & SQP::statistics() const
  {
    return stats_;
  }

  bool SQP::checkKKT(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, 
                     const LinearConstraints& lc, const LeastSquareObjective* const obj) const
  {
    //see Stan's thesis $4.3.5
    double tx = params_.tau_p()*(1 + x.lpNorm<Infinity>());
    double tl = params_.tau_d()*(1 + lambda.lpNorm<Infinity>());

    lc.transposeMult(Cl_, lambda);
    gradL_ = f*g + Cl_;
    if (obj)
    {
      obj->applyJToTheLeft(Jx_, x);
      obj->applyJTransposeToTheLeft(JtJx_, Jx_);
      gradL_ += JtJx_;
    }

    if (gradL_.lpNorm<Infinity>() <= tl)
    {
      lc.mult(Cx_, x);
      bool b = true;
      const auto& l = lc.l();
      const auto& u = lc.u();
      for (DenseIndex i = 0; b && i <= n_; ++i)
      {
        b = (std::abs(Cx_[i] - l[i]) <= tx && lambda[i] <= tl)
          || (l[i] - tx <= Cx_[i] && Cx_[i] <= u[i] + tx && std::abs(lambda[i]) <= tl)
          || (std::abs(Cx_[i] - u[i]) <= tx && lambda[i] >= -tl);
      }
      return b;
    }
    else
      return false;
  }
}

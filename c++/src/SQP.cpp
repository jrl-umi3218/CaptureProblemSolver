#include <SQP.h>

#include <iostream>

using namespace Eigen;

namespace bms
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


  SolverStatus SQP::solve(const Problem& pb)
  {
    //shortcuts
    const auto& lc = pb.linearConstraints();
    const auto& nlc = pb.nonLinearConstraint();
    const auto& obj = pb.objective();

    //active sets
    currentActiveSet_ = lc.activationStatus();
    previousActiveSet_ = lc.activationStatus();

    //penality parameter
    double mu = 100000;
    double mu2 = mu*mu;

    //feasible (x,lambda)
    FeasiblePointInfo fpi;
    std::tie(fpi, x_) = lc.initialPoint(true);
    if (fpi != FeasiblePointInfo::Found) //FIXME: we can be more precise here
    {
      //attempt to find a point without the active set
      std::tie(fpi, x_) = lc.initialPoint(false);
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
    for (int k = 0; k < maxIter_; ++k)
    {
      double f;
      nlc.compute(f, j_, x_);
      f *= mu;
      j_ *= mu;
      obj.applyJToTheLeft(Jx_, x_);
      double v = 0.5* Jx_.squaredNorm();
      //std::cout << "1/2||f||^2" << 0.5*f*f << std::endl;
      //std::cout << "obj = " << 0.5*f*f + v << std::endl;

      if (checkKKT(x_, lambda_, f, j_, lc, &obj))
        return SolverStatus::Converge;

      //solve the least square
      shiftedLC_ = lc.shift(x_, true);
      previousActiveSet_ = currentActiveSet_;
      shiftedLC_.setActivationStatus(currentActiveSet_);
      ls_.solve(obj, Jx_, j_, f, shiftedLC_);
      currentActiveSet_ = shiftedLC_.activationStatus();
      const auto& p = ls_.x();
      const auto& pl = ls_.lambda();

      //line search
      double a = 1;
      obj.applyJToTheLeft(Jx_, x_);
      obj.applyJToTheLeft(Jp_, p);
      double cgp = c1*(f*j_.dot(p) + Jx_.dot(Jp_));

      xa_ = x_ + a*p;
      double fa;
      nlc.compute(fa, xa_);
      double va = obj.value(xa_);
      double val0 = 0.5*f*f + v;
      while (0.5*mu2*fa*fa + va > val0 + a*cgp)
      {
        a = beta * a;
        if (a < smallestLSStep)
          return SolverStatus::StepTooSmall;
        xa_ = x_ + a*p;
        nlc.compute(fa, xa_);
        va = obj.value(xa_);
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
    //shortcuts
    const auto& lc = pb.linearConstraints();
    const auto& nlc = pb.nonLinearConstraint();

    //feasible (x,lambda)
    FeasiblePointInfo fpi;
    std::tie(fpi, x_) = lc.initialPoint();
    if (fpi != FeasiblePointInfo::Found) //FIXME: we can be more precise here
      return SolverStatus::NoLinearlyFeasiblePoint; 
    lambda_.setZero();

    //main loop
    for (int k = 0; k < maxIter_; ++k)
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
      double cgp = c1*f*j_.dot(p);

      xa_ = x_ + a*p;
      double fa;
      nlc.compute(fa, xa_);
      while (0.5*fa*fa > 0.5*f*f + a*cgp)
      {
        a = beta * a;
        if (a < smallestLSStep)
          return SolverStatus::StepTooSmall;
        xa_ = x_ + a*p;
        nlc.compute(fa, xa_);
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

  bool SQP::checkKKT(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, 
                     const LinearConstraints& lc, const LeastSquareObjective* const obj) const
  {
    //see Stan's thesis $4.3.5
    double tx = tau_p*(1 + x.lpNorm<Infinity>());
    double tl = tau_d*(1 + lambda.lpNorm<Infinity>());

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

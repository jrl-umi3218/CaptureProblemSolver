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
    , lambda_(n+1)
    , Cx_(n+1)
    , Cl_(n)
    , gradL_(n)
  {
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

      if (checkKKTFeasibility(x_, lambda_, f, j_, lc))
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
      lambda_ += a*(lambda_ - pl);
    }

    return SolverStatus();
  }

  const Eigen::VectorXd & SQP::x() const
  {
    return x_;
  }
  const Eigen::VectorXd & SQP::lambda() const
  {
    return lambda_;
  }

  bool SQP::checkKKTFeasibility(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, const LinearConstraints& lc) const
  {
    //see Stan's thesis $4.3.5
    double tx = tau_p*(1 + x.lpNorm<Infinity>());
    double tl = tau_d*(1 + lambda.lpNorm<Infinity>());

    lc.transposeMult(Cl_, lambda);
    gradL_ = f*g + Cl_;

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

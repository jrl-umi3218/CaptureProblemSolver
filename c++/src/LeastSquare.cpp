#include "LeastSquare.h"

using namespace Eigen;

namespace bms
{
  LeastSquare::LeastSquare(int n)
    : n_(10)
    , maxIter_(10*n)
    , x_(n)
    , p_(n)
    , z_(n)
    , lambda_(n + 1)
    , lambdaAct_(n + 1)
    , jN_(n)
  {
  }

  LSStatus LeastSquare::solveFeasibility(const VectorConstRef& j, double c, LinearConstraints& lc)
  {
    assert(j.size() == n_);
    assert(lc.size() == n_);

    x_.setZero();
    assert(lc.checkPrimal(x_) && "This solver supposes that 0 is a feasible point");

    //main loop
    int k = 0;
    for (k = 0; k < maxIter_; ++k)
    {
      auto jN = jN_.head(lc.nullSpaceSize());
      auto z = z_.head(lc.nullSpaceSize());
      lc.applyNullSpaceOnTheRight(jN.transpose(), j.transpose());   //jN^T = j^T N
      z = -(c+j.dot(x_))/jN.squaredNorm() * jN;                     //-pinv(jN^T*(c+j^T*x)
      lc.applyNullSpaceOnTheLeft(p_, z);                            //p = Nz

      if (p_.lpNorm<Infinity>() < 1e-12)  //FIXME: hard coded threshold?
      {
        //compute Lagrange multipliers
        auto lambdaAct = lambdaAct_.head(lc.numberOfActiveConstraints());
        lc.pinvTransposeMultAct(lambdaAct, j);
        lambdaAct *= -(c + j.dot(x_));
        lc.expandActive(lambda_, lambdaAct);
        if (lc.checkDual(lambda_))
          return LSStatus::Converge;
        else
          lc.deactivateMaxLambda(lambda_);
      }
      else
      {
        x_ = lc.performQPstep(x_, p_);
      }
    }

    return LSStatus::MaxIteration;
  }

  const Eigen::VectorXd & LeastSquare::x() const
  {
    return x_;
  }

  const Eigen::VectorXd & LeastSquare::lambda() const
  {
    return lambda_;
  }
}
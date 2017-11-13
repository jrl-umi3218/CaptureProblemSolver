#include "LeastSquare.h"
#include "QRAlgorithms.h"

#include <Eigen/SVD>

using namespace Eigen;

namespace bms
{
  LeastSquare::LeastSquare(int n)
    : n_(n)
    , x_(n)
    , p_(n)
    , z_(n)
    , lambda_(n + 1)
    , lambdaAct_(n + 1)
    , jN_(n)
    , A_(n,n)
    , b_(n)
    , Jx_(n-1)
    , JtJx_(n)
    , tmp_(n)
    , Q_(n,n,2*n,true)
  {
    Qg_.reserve(n);
  }

  void LeastSquare::parameters(const LeastSquare::Parameters& param)
  {
    params_ = param;
  }

  const LeastSquare::Parameters& LeastSquare::parameters() const
  {
    return params_;
  }

  SolverStatus LeastSquare::solve(const LeastSquareObjective& obj, const VectorConstRef& Jx0, const VectorConstRef& j, double c, LinearConstraints& lc)
  {
    assert(j.size() == n_);
    assert(obj.size() == n_);
    assert(Jx0.size() == n_ - 1);
    assert(lc.size() == n_);

    x_.setZero();
    assert(lc.checkPrimal(x_) && "This solver supposes that 0 is a feasible point");

    STATISTICS(stats_.reset());

    //main loop
    int k = 0;
    if (params_.maxIter() <= 0) { params_.maxIter(10 * static_cast<int>(n_)); } //default value for maxIter if need be.
    bool skip = false;
    for (k = 0; k < params_.maxIter(); ++k)
    {
      auto nz = lc.nullSpaceSize();
      if (nz > 0 && !skip)
      {
        auto A = A_.leftCols(nz);
        auto jN = jN_.head(nz);
        auto z = z_.head(nz);
        //We form A = [j^T * N; R_A] where R_A is the obtain by the qr of J_A
        lc.applyNullSpaceOnTheRight(jN.transpose(), j.transpose());   //jN^T = j^T N
        A.row(0) = jN.transpose();
        obj.qr(A.bottomRows(n_ - 1), Q_, lc.numberOfActiveConstraints(), lc.activeSet(), 1);
        //Compute b
        obj.applyJToTheLeft(b_.tail(n_ - 1), x_);
        b_.tail(n_ - 1) += Jx0;
        b_[0] = c + j.dot(x_);
        //QR of A (A is upper Hessenberg
        bool fullRank = hessenbergQR(A.topLeftCorner(std::min(n_, nz + 1), nz), Qg_, false, params_.rankThreshold());
        //std::cout << "A= \n" << A << std::endl;
        //z = A^-1 b
        Q_.applyTo(b_);
        Qg_.applyTo(b_);
        //std::cout << "b = " << b_.transpose() << std::endl; 
        if (fullRank)
        {
          z = -b_.head(nz);
          A.topLeftCorner(nz, nz).triangularView<Upper>().solveInPlace(z);
        }
        else
        {
          z.head(nz - 1) = -b_.head(nz - 1);
          A.topLeftCorner(nz - 1, nz - 1).triangularView<Upper>().solveInPlace(z.head(nz - 1));
          z[nz - 1] = 0;
          STATISTICS(++stats_.rankLoss);
        }
        lc.applyNullSpaceOnTheLeft(p_, z);                            //p = Nz
      }
      else
        p_.setZero();

      if (p_.lpNorm<Infinity>() < params_.minNorm_p())  //FIXME: hard coded threshold?
      {
        //compute Lagrange multipliers for the active constraints: -C_A^+T(A^T(Ax+[c;Jx0]))
        obj.applyJToTheLeft(Jx_, x_);
        Jx_ += Jx0;
        obj.applyJTransposeToTheLeft(JtJx_, Jx_);
        tmp_ = -JtJx_ - (c + j.dot(x_))*j;
        auto lambdaAct = lambdaAct_.head(lc.numberOfActiveConstraints());
        lc.pinvTransposeMultAct(lambdaAct, tmp_);
        //std::cout << "lambda_act = " << lambdaAct.transpose() << std::endl;
        //std::cout << "lambda_pinv = " << lc.matrixAct().transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(tmp_) << std::endl;
        //std::cout << "err = " << (lc.matrixAct().transpose()*lambdaAct - tmp_).transpose() << std::endl;
        //std::cout << "optim = " << ((c + j.dot(x_))*j + obj.matrix().transpose()*(obj.matrix()*x_) + lc.matrixAct().transpose()*lambdaAct).transpose() << std::endl;
        lc.expandActive(lambda_, lambdaAct);
        if (lc.checkDual(lambda_, params_.dualEps()))
        {
          STATISTICS(stats_.iter = k + 1);
          STATISTICS(stats_.activeConstraints = (int)lc.numberOfActiveConstraints());
          return SolverStatus::Converge;
        }
        else
        {
          lc.deactivateMaxLambda(lambda_);
          STATISTICS(++stats_.deactivation);
        }
        skip = false;
      }
      else
      {
        lc.performQPstep(x_, p_, &skip);
        STATISTICS(if (!skip) { ++stats_.activation; });
      }
    }

    STATISTICS(stats_.iter = k);
    STATISTICS(stats_.activeConstraints = (int)lc.numberOfActiveConstraints());
    return SolverStatus::MaxIteration;
  }

  SolverStatus LeastSquare::solveFeasibility(const VectorConstRef& j, double c, LinearConstraints& lc)
  {
    assert(j.size() == n_);
    assert(lc.size() == n_);

    x_.setZero();
    assert(lc.checkPrimal(x_) && "This solver supposes that 0 is a feasible point");

    //main loop
    int k = 0;
    if (params_.maxIter() <= 0) { params_.maxIter(10 * static_cast<int>(n_)); } //default value for maxIter if need be.
    for (k = 0; k < params_.maxIter(); ++k)
    {
      auto jN = jN_.head(lc.nullSpaceSize());
      auto z = z_.head(lc.nullSpaceSize());
      lc.applyNullSpaceOnTheRight(jN.transpose(), j.transpose());   //jN^T = j^T N
      z = -(c+j.dot(x_))/jN.squaredNorm() * jN;                     //-pinv(jN^T*(c+j^T*x)
      lc.applyNullSpaceOnTheLeft(p_, z);                            //p = Nz

      if (p_.lpNorm<Infinity>() < params_.minNorm_p())  //FIXME: hard coded threshold?
      {
        //compute Lagrange multipliers
        auto lambdaAct = lambdaAct_.head(lc.numberOfActiveConstraints());
        lc.pinvTransposeMultAct(lambdaAct, j);
        lambdaAct *= -(c + j.dot(x_));
        lc.expandActive(lambda_, lambdaAct);
        if (lc.checkDual(lambda_, params_.dualEps()))
          return SolverStatus::Converge;
        else
          lc.deactivateMaxLambda(lambda_);
      }
      else
      {
        lc.performQPstep(x_, p_);
      }
    }

    return SolverStatus::MaxIteration;
  }

  const Eigen::VectorXd & LeastSquare::x() const
  {
    return x_;
  }

  const Eigen::VectorXd & LeastSquare::lambda() const
  {
    return lambda_;
  }

  const stats::LSStats & LeastSquare::statistics() const
  {
    return stats_;
  }
}
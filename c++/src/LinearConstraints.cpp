#include "LinearConstraints.h"

using namespace Eigen;


namespace
{
  void applyFirst(bms::MatrixRef Y, const bms::MatrixConstRef& X, DenseIndex k)
  {
    assert(Y.rows() == k && Y.cols() == X.cols());
    Y.setZero();
    Y.row(k - 1) = X.row(k - 1);
    for (auto i = k - 2; i >= 0; --i)
      Y.row(i) = Y.row(i + 1) + X.row(i);
  }

  void applyLast(bms::MatrixRef Y, const bms::MatrixConstRef& X, DenseIndex k)
  {
    assert(Y.rows() == k && Y.cols() == X.cols());
    if (k == 1)
      Y = X.bottomRows<1>();
    else
    {
      auto n = X.rows();
      Y.row(0) = -X.row(n - k);
      for (DenseIndex i = 1; i < k - 1; ++i)
        Y.row(i) = Y.row(i - 1) - X.row(n - k + i);
      Y.row(k - 1) = -Y.row(k - 2) + X.bottomRows<1>();
    }
  }

  void apply(bms::MatrixRef Y, const bms::MatrixConstRef& X, DenseIndex s, DenseIndex k)
  {
    assert(Y.rows() == k - 1 && Y.cols() == X.cols());
    double d = 1. / k;
    for (DenseIndex i = 0; i < k - 1; ++i)
    {
      double m = (i + 1 - k)*d;
      Y.row(i) = m * X.row(s);
      for (DenseIndex j = 1; j < i + 1; ++j)
        Y.row(i) += m * X.row(s + j);
      m = (i+1)*d;
      for (DenseIndex j = i + 1; j < k; ++j)
        Y.row(i) += m*X.row(s + j);
    }
  }
}

namespace bms
{
  LinearConstraints::LinearConstraints(const Eigen::VectorXd& l, const Eigen::VectorXd& u, double xln, double xun)
    : n_(l.size())
    , l_(l.size()+1)
    , u_(l.size()+1)
    , na_(0)
    , activationStatus_(l.size()+1, Activation::None)
    , activeSet_(l.size()+1, false)
    , validIdx_(false)
    , validActIdx_(false)
    , idx_(l.size(),0)
    , Cx_(l.size()+1)
    , Cp_(l.size()+1)
  {
    assert(l.size() == u.size());
    assert(((u - l).array() >= 0).all());

    l_ << l, xln;
    u_ << u, xun;

    actIdx_.reserve(n_);

    for (int i = 0; i < n_+1; ++i)
    {
      if (l_[i] == u_[i])
        activate(i, Activation::Equal);
    }
  }

  DenseIndex LinearConstraints::size() const
  {
    return n_;
  }

  DenseIndex LinearConstraints::nullSpaceSize() const
  {
    return n_ - na_;
  }

  void LinearConstraints::activate(size_t i, Activation a)
  {
    assert(0 <= i && i <= static_cast<size_t>(n_));
    activationStatus_[i] = a;
    bool b = a != Activation::None;
    if (activeSet_[i] && !b)
      --na_;
    else if (!activeSet_[i] && b)
      ++na_;
    activeSet_[i] = b;
    validIdx_ = false;
    validActIdx_ = false;
  }

  void LinearConstraints::deactivate(size_t i)
  {
    assert(0 <= i && i <= static_cast<size_t>(n_));
    activationStatus_[i] = Activation::None;
    if (activeSet_[i])
      --na_;
    activeSet_[i] = false;
    validIdx_ = false;
    validActIdx_ = false;
  }

  Activation LinearConstraints::activationStatus(size_t i) const
  {
    return activationStatus_[i];
  }

  const std::vector<bool>& LinearConstraints::activeSet() const
  {
    return activeSet_;
  }

  const std::vector<Eigen::DenseIndex>& LinearConstraints::activeSetIdx() const
  {
    if (!validActIdx_)
      computeActIdx();
    return actIdx_;
  }

  Eigen::DenseIndex LinearConstraints::numberOfActiveConstraints() const
  {
    return na_;
  }

  void LinearConstraints::applyNullSpaceOnTheRight(MatrixRef Y, const MatrixConstRef& X) const
  {
    assert(X.cols() == n_);
    assert(Y.rows() == X.rows() && Y.cols() == n_ - na_);

    DenseIndex c = -1;
    DenseIndex k = 0;
    std::fill(idx_.begin(), idx_.end(), -1);
    actIdx_.resize(na_);

    //skip the first group of activated constraints
    while (k < n_ && activeSet_[static_cast<size_t>(k)])
    {
      actIdx_[static_cast<size_t>(k)] = k;
      ++k;
    }

    size_t ca = static_cast<size_t>(k);
    auto i = k;
    for (; i < n_; ++i)
    {
      if (activeSet_[static_cast<size_t>(i)])
      {
        Y.col(c) += X.col(i);
        actIdx_[ca] = i;
        ++ca;
      }
      else
      {
        ++c;
        if (c >= n_ - na_)
          break;
        Y.col(c) = X.col(i);
      }
      idx_[static_cast<size_t>(i)] = c;
    }
    for (++i; i < n_; ++i)
    {
      actIdx_[ca] = i;
      ++ca;
    }
    if (activeSet_.back())
      actIdx_[ca] = n_;

    validIdx_ = true;
    validActIdx_ = true;
  }

  void LinearConstraints::applyNullSpaceOnTheLeft(MatrixRef Y, const MatrixConstRef& X) const
  {
    assert(X.rows() == n_ - na_);
    assert(Y.rows() == n_ && Y.cols() == X.cols());

    if (!validIdx_)
      computeIdx();

    for (DenseIndex i = 0; i < n_; ++i)
    {
      auto j = idx_[static_cast<size_t>(i)];
      if (j >= 0)
        Y.row(i) = X.row(j);
      else
        Y.row(i).setZero();
    }
  }

  void LinearConstraints::mult(MatrixRef Y, const MatrixConstRef& X) const
  {
    assert(X.rows() == n_);
    assert(Y.rows() == n_ + 1 && Y.cols() == X.cols());

    Y.topRows(n_) = X;
    Y.middleRows(1, n_ - 1) -= X.topRows(n_ - 1);
    Y.bottomRows<1>() = X.bottomRows<1>();
  }

  void LinearConstraints::transposeMult(MatrixRef Y, const MatrixConstRef& X) const
  {
    assert(X.rows() == n_ + 1);
    assert(Y.rows() == n_ && Y.cols() == X.cols());

    Y = X.topRows(n_);
    Y.topRows(n_ - 1) -= X.middleRows(1, n_ - 1);
    Y.bottomRows<1>() += X.bottomRows<1>();
  }

  void LinearConstraints::pinvTransposeMultAct(MatrixRef Y, const MatrixConstRef& X) const
  {
    //Fixme: we are re-scanning the active constraint here. Should we maintain some data to avoid making checks all the time ?
    DenseIndex k0, kn;
    if (activeSet_[0])
    {
      DenseIndex k = 1;
      while (k < na_ && activeSet_[static_cast<size_t>(k)]) { ++k; }
      applyFirst(Y.topRows(k), X, k);
      k0 = k;
    }
    else
      k0 = 0;

    if (activeSet_[n_])
    {
      DenseIndex k = 1;
      while (k < na_ && activeSet_[static_cast<size_t>(n_ - k)]) { ++k; }
      applyLast(Y.bottomRows(k), X, k);
      kn = n_ - k;
    }
    else
      kn = n_;

    DenseIndex k = 0;
    for (auto i = k0; i < kn; ++i)
    {
      if (activeSet_[static_cast<size_t>(i)])
        ++k;
      else
      {
        if (k > 0)
        {
          apply(Y.middleRows(k0, k), X, i - k - 1, k + 1);
          k0 += k;
          k = 0;
        }
      }
    }
    if (k>0)
      apply(Y.middleRows(k0, k), X, kn - k - 1, k + 1);
  }

  void LinearConstraints::expandActive(VectorRef y, const VectorConstRef& x) const
  {
    assert(y.size() == n_ + 1);
    assert(x.size() == na_);

    if (!validActIdx_)
      computeActIdx();

    y.setZero();
    for (DenseIndex i = 0; i < na_; ++i)
      y[actIdx_[static_cast<size_t>(i)]] = x[i];
  }

  Eigen::VectorXd LinearConstraints::performQPstep(const Eigen::VectorXd& x, const Eigen::VectorXd& p)
  {
    mult(Cx_, x);  //Cx = C*x
    mult(Cp_, p);  //Cp = C*p

    double a = std::numeric_limits<double>::max();
    size_t iact = -1;
    Activation type = Activation::None;

    //Looking for the maximum a such that x+a*p is still feasible
    //a must satisfy [l;xln] <= C(x+a*p) <= [u;xun]
    for (DenseIndex i = 0; i <= n_; ++i)
    {
      double ai;
      Activation typei;
      if (!activeSet_[static_cast<size_t>(i)] && Cp_[i] > 0)
      {
        assert(std::abs(Cp_[i]) > 1e-12);
        ai = (u_[i] - Cx_[i]) / Cp_[i];
        typei = Activation::Upper;
      }
      else if (!activeSet_[static_cast<size_t>(i)] && Cp_[i] < 0)
      {
        assert(std::abs(Cp_[i]) > 1e-12);
        ai = (l_[i] - Cx_[i]) / Cp_[i];
        typei = Activation::Lower;
      }
      if (ai < a)
      {
        a = ai;
        iact = i;
        type = typei;
      }
    }

    if (a < 0)
    {
      assert(a > -1e-14); //FIXME: relative value ?
      a = 0;
    }

    if (a < 1)
    {
      activate(iact, type);
      return x + a*p;
    }
    else
      return x+p;
  }

  void LinearConstraints::deactivateMaxLambda(const VectorConstRef& lambda)
  {
    if (!validActIdx_)
      computeActIdx();
    
    double lmax = 0;
    size_t imax;

    for (size_t i = 0; i < actIdx_.size(); ++i)
    {
      auto k = actIdx_[i];
      switch (activationStatus_[static_cast<size_t>(k)])
      {
      case Activation::Lower: if (lambda[k] > lmax) { lmax = lambda[k]; imax = i; } break;
      case Activation::Upper: if (-lambda[k] > lmax) { lmax = -lambda[k]; imax = i; } break;
      default: break;
      }
    }

    deactivate(imax);
  }


  Eigen::MatrixXd LinearConstraints::matrix() const
  {
    MatrixXd C = MatrixXd::Identity(n_ + 1, n_);
    C.diagonal<-1>().setConstant(-1);
    C(n_, n_ - 1) = 1;

    return C;
  }

  Eigen::MatrixXd LinearConstraints::matrixAct() const
  {
    MatrixXd Ca = MatrixXd::Zero(na_, n_);

    DenseIndex r = 0;
    if (activeSet_[0])
    {
      Ca(0, 0) = 1;
      ++r;
    }

    for (DenseIndex i = 1; i < n_; ++i)
    {
      if (activeSet_[static_cast<size_t>(i)])
      {
        Ca(r, i - 1) = -1;
        Ca(r, i) = 1;
        ++r;
      }
    }

    if (activeSet_[static_cast<size_t>(n_)])
      Ca(r, n_ - 1) = 1;

    return Ca;
  }

  bool LinearConstraints::checkPrimal(const VectorConstRef& x, double eps) const
  {
    assert(x.size() == n_);
    mult(Cx_, x);

    return ((Cx_.array() >= (l_.array() - eps)).all() && (Cx_.array() <= (u_.array() + eps)).all());
  }

  bool LinearConstraints::checkDual(const VectorConstRef& lambda, double eps) const
  {
    if (!validActIdx_)
      computeActIdx();

    bool b = true;
    for (size_t i = 0; i < actIdx_.size() && b; ++i)
    {
      auto k = actIdx_[i];
      switch (activationStatus_[static_cast<size_t>(k)])
      {
      case Activation::Lower: b = lambda[k] <=  eps; break;
      case Activation::Upper: b = lambda[k] >= -eps; break;
      default: break;
      }
    }

    return b;
  }


  void LinearConstraints::computeActIdx() const
  {
    actIdx_.resize(na_);
    size_t ca = 0;
    for (DenseIndex i=0; i <= n_; ++i)
    {
      if (activeSet_[static_cast<size_t>(i)])
      {
        actIdx_[ca] = i;
        ++ca;
      }
    }
    validActIdx_ = true;
  }

  void LinearConstraints::computeIdx() const
  {
    DenseIndex c = -1;
    DenseIndex k = 0;
    std::fill(idx_.begin(), idx_.end(), -1);

    //skip the first group of activated constraints
    while (k <n_ && activeSet_[static_cast<size_t>(k)])
      ++k;

    for (auto i = k; i < n_; ++i)
    {
      if (!activeSet_[static_cast<size_t>(i)])
      {
        ++c;
        if (c >= n_ - na_)
          break;
      }
      idx_[static_cast<size_t>(i)] = c;
    }
    validIdx_ = true;
  }
}
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
    , l_(l)
    , u_(u)
    , xln_(xln)
    , xun_(xun)
    , na_(0)
    , activationStatus_(n_+1,Activation::None)
    , activeSet_(n_+1,false)
    , validIdx_(false)
    , idx_(n_,0)
  {
    assert(l.size() == u.size());
    assert(((u - l).array() >= 0).all());

    for (int i = 0; i < n_; ++i)
    {
      if (l[i] == u[i])
        activate(i, Activation::Equal);
    }
    if (xln == xun)
      activate(n_, Activation::Equal);
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
  }

  void LinearConstraints::deactivate(size_t i)
  {
    assert(0 <= i && i <= static_cast<size_t>(n_));
    activationStatus_[i] = Activation::None;
    if (activeSet_[i])
      --na_;
    activeSet_[i] = false;
    validIdx_ = false;
  }

  Activation LinearConstraints::activationStatus(size_t i) const
  {
    return activationStatus_[i];
  }

  const std::vector<bool>& LinearConstraints::activeSet() const
  {
    return activeSet_;
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

    //skip the first group of activated constraints
    while (k <n_ && activeSet_[static_cast<size_t>(k)])
      ++k;

    for (auto i = k; k < n_; ++i)
    {
      if (activeSet_[static_cast<size_t>(i)])
        Y.col(c) += X.col(i);
      else
      {
        ++c;
        if (c >= n_ - na_)
          break;
        Y.col(c) = X.col(i);
      }
      idx_[static_cast<size_t>(i)] = c;
    }
    validIdx_ = true;
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

  void LinearConstraints::pinvTransposeMult(MatrixRef Y, const MatrixConstRef& X)
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

  Eigen::VectorXd LinearConstraints::performQPstep(const Eigen::VectorXd& x, const Eigen::VectorXd& p)
  {
    return Eigen::VectorXd();
  }

  Eigen::MatrixXd LinearConstraints::matrix() const
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


  void LinearConstraints::computeIdx() const
  {
    DenseIndex c = -1;
    DenseIndex k = 0;
    std::fill(idx_.begin(), idx_.end(), -1);

    //skip the first group of activated constraints
    while (k <n_ && activeSet_[static_cast<size_t>(k)])
      ++k;

    for (auto i = k; k < n_; ++i)
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
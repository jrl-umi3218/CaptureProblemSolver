#include <iostream>

#include <bms/LinearConstraints.h>

using namespace Eigen;

namespace
{
  template<typename Derived1, typename Derived2>
  void applyFirst(const MatrixBase<Derived1>& Y_, const MatrixBase<Derived2>& X, DenseIndex k)
  {
    Eigen::MatrixBase<Derived1>& Y = const_cast<Eigen::MatrixBase<Derived1>&>(Y_);
    assert(Y.rows() == k && Y.cols() == X.cols());
    Y.setZero();
    Y.row(k - 1) = X.row(k - 1);
    for (auto i = k - 2; i >= 0; --i)
      Y.row(i) = Y.row(i + 1) + X.row(i);
  }

  template<typename Derived1, typename Derived2>
  void applyLast(const MatrixBase<Derived1>& Y_, const MatrixBase<Derived2>& X, DenseIndex k)
  {
    Eigen::MatrixBase<Derived1>& Y = const_cast<Eigen::MatrixBase<Derived1>&>(Y_);
    assert(Y.rows() == k && Y.cols() == X.cols());
    if (k == 1)
      Y = X.template bottomRows<1>();
    else
    {
      auto n = X.rows();
      Y.row(0) = -X.row(n - k);
      for (DenseIndex i = 1; i < k - 1; ++i)
        Y.row(i) = Y.row(i - 1) - X.row(n - k + i);
      Y.row(k - 1) = -Y.row(k - 2) + X.template bottomRows<1>();
    }
  }

  template<typename Derived1, typename Derived2>
  void apply(const MatrixBase<Derived1>& Y_, const MatrixBase<Derived2>& X, DenseIndex s, DenseIndex k)
  {
    Eigen::MatrixBase<Derived1>& Y = const_cast<Eigen::MatrixBase<Derived1>&>(Y_);
    assert(Y.rows() == k - 1 && Y.cols() == X.cols());
    double d = 1. / double(k);
    for (DenseIndex i = 0; i < k - 1; ++i)
    {
      double m = double(i + 1 - k)*d;
      Y.row(i) = m * X.row(s);
      for (DenseIndex j = 1; j < i + 1; ++j)
        Y.row(i) += m * X.row(s + j);
      m = double(i+1)*d;
      for (DenseIndex j = i + 1; j < k; ++j)
        Y.row(i) += m*X.row(s + j);
    }
  }
}

namespace bms
{
  LinearConstraints::LinearConstraints(int n)
    : n_(n)
    , l_(n+ 1)
    , u_(n + 1)
    , na_(0)
    , activationStatus_(n + 1, Activation::None)
    , activeSet_(n + 1, false)
    , validIdx_(false)
    , validActIdx_(false)
    , idx_(n, 0)
    , Cx_(n + 1)
    , Cp_(n + 1)
    , actl_(n + 1)
    , actu_(n + 1)
  {
  }

  LinearConstraints::LinearConstraints(const Eigen::VectorXd& l, const Eigen::VectorXd& u, double xln, double xun)
    : LinearConstraints(static_cast<int>(l.size()))
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

  LinearConstraints& LinearConstraints::operator= (const LinearConstraints::Shift& s)
  {
    this->operator=(s.c);
    //auto xm = Ref<const VectorXd>(s.x);
    mult(Cx_, s.x);

    l_ -= Cx_;
    u_ -= Cx_;

    if (s.feasible)
    {
      for (DenseIndex i = 0; i <= n_; ++i)
      {
        if (activationStatus_[static_cast<size_t>(i)]==Activation::Equal)
        {
          l_[i] = 0.;
          u_[i] = 0.;
        }
        else
        {
          l_[i] = std::min(l_[i], 0.);
          u_[i] = std::max(u_[i], 0.);
        }
      }
    }

    return *this;
  }

  LinearConstraints::Shift LinearConstraints::shift(const Eigen::VectorXd& x, bool feasible) const
  {
    assert(!feasible || checkPrimal(x, 1e-14));
    return{ *this, x, feasible };
  }

  const Eigen::VectorXd & LinearConstraints::l() const
  {
    return l_;
  }

  const Eigen::VectorXd & LinearConstraints::u() const
  {
    return u_;
  }

  std::pair<FeasiblePointInfo, VectorXd> LinearConstraints::initialPoint(bool takeActivationIntoAccount, double eps) const
  {
    std::pair<FeasiblePointInfo, VectorXd> ret;
    ret.second.resize(n_);

    //take into account the activation status if applicable
    if (takeActivationIntoAccount)
    {
      for (DenseIndex i = 0; i <= n_; ++i)
      {
        switch (activationStatus_[static_cast<int>(i)])
        {
        case Activation::None:  actl_[i] = l_[i]; actu_[i] = u_[i]; break;
        case Activation::Lower: actl_[i] = l_[i]; actu_[i] = l_[i]; break;
        case Activation::Upper: actl_[i] = u_[i]; actu_[i] = u_[i]; break;
        case Activation::Equal: actl_[i] = l_[i]; actu_[i] = u_[i]; break;
        }
      }
    }
    else
    {
      actl_ = l_;
      actu_ = u_;
    }

    double s = 0;
    double t = 0;
    //s = sum(l), t = sum(u-l)
    for (DenseIndex i = 0; i < n_; ++i)
    {
      s += actl_[i];
      t += actu_[i];
    }
    t -= s;

    if (std::abs(t) > double(n_)*eps)
    {
      double al = (actl_[n_] - s) / t;
      double au = (actu_[n_] - s) / t;
      if (al <= 1 && au >= 0)
      {
        double f = (std::max(al, 0.) + std::min(au, 1.)) / 2;
        ret.second[0] = actl_[0] + f*(actu_[0] - actl_[0]);
        for (DenseIndex i = 1; i < n_; ++i)
          ret.second[i] = ret.second[i - 1] + actl_[i] + f*(actu_[i] - actl_[i]);

        if (al > 1 - eps || au < eps)
          ret.first = FeasiblePointInfo::TooSmallFeasibilityZone;
        else
        {
          if (std::abs(al) < eps)
          {
            if (std::abs(au - 1) < eps)
              ret.first = FeasiblePointInfo::NumericalWarningBoth;
            else
              ret.first = FeasiblePointInfo::NumericalWarningLower;
          }
          else
          {
            if (std::abs(au - 1) < eps)
              ret.first = FeasiblePointInfo::NumericalWarningUpper;
            else
              ret.first = FeasiblePointInfo::Found;
          }
        }
      }
      else
        ret.first = FeasiblePointInfo::Infeasible;
    }
    else
    {
      ret.second[0] = actl_[0];
      for (DenseIndex i = 1; i < n_; ++i)
        ret.second[i] = ret.second[i - 1] + actl_[i];
      if (l_[n_] <= ret.second[n_ - 1] && ret.second[n_ - 1] <= u_[n_])
        ret.first = FeasiblePointInfo::TooSmallZonotope;
      else
        ret.first = FeasiblePointInfo::Infeasible;
    }

    return ret;
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

  void LinearConstraints::resetActivation()
  {
    for (size_t i = 0; i < static_cast<size_t>(n_); ++i)
    {
      if (activationStatus_[i] != Activation::Equal)
        deactivate(i);
    }
  }

  Activation LinearConstraints::activationStatus(size_t i) const
  {
    return activationStatus_[i];
  }

  const std::vector<bool>& LinearConstraints::activeSet() const
  {
    return activeSet_;
  }

  const std::vector<Activation>& LinearConstraints::activationStatus() const
  {
    return activationStatus_;
  }

  void LinearConstraints::setActivationStatus(const std::vector<Activation>& act)
  {
    assert(act.size() == activationStatus_.size());
    for (size_t i = 0; i < act.size(); ++i)
    {
      if (activationStatus_[i] != Activation::Equal)
      {
        if (act[i] != Activation::Equal)
          activate(i, act[i]);
        else
          throw std::runtime_error("You can't activate constraint as equality. This is defined by the value of the bounds.");
      }
      else
      {
        if (act[i] != Activation::Equal)
          throw std::runtime_error("You can't change the status of an equality constraint. This is defined by the value of the bounds.");
      }
    }
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

  void LinearConstraints::changeBounds(DenseIndex i, double l, double u)
  {
    assert(i >= 0 && i <= n_);
    l_[i] = l;
    u_[i] = u;
    if (l == u)
      activate(static_cast<int>(i), Activation::Equal);
    else
      deactivate(static_cast<int>(i));
  }

  void LinearConstraints::changeBounds(const VectorConstRef l, const VectorConstRef u)
  {
    assert(l.size() == n_ || l.size() == n_ + 1);
    assert(l.size() == u.size());
    for (DenseIndex i = 0; i < l.size(); ++i)
      changeBounds(i, l[i], u[i]);
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

  void LinearConstraints::performQPstep(VectorXd& x, const VectorXd& p, bool* fullstep)
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
      if (!activeSet_[static_cast<size_t>(i)])
      {
        assert(std::abs(Cp_[i]) > 1e-12);
        double ai;
        Activation typei;
        if (Cp_[i] > 0)
        {
          ai = (u_[i] - Cx_[i]) / Cp_[i];
          typei = Activation::Upper;
        }
        else
        {
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
    }

    if (a < 0)
    {
      assert(a > -1e-14); //FIXME: relative value ?
      a = 0;
    }

    if (a < 1)
    {
      //std::cout << "activate " << iact << std::endl;
      if (fullstep)
        *fullstep = false;
      activate(iact, type);
      x += a*p;
    }
    else
    {
      //std::cout << "full step of norm " << p.lpNorm<Infinity>() << std::endl;
      if (fullstep)
        *fullstep = true;
      x += p;
    }
  }

  void LinearConstraints::deactivateMaxLambda(const VectorConstRef& lambda)
  {
    if (!validActIdx_)
      computeActIdx();
    
    double lmax = 0;
    DenseIndex kmax;

    for (size_t i = 0; i < actIdx_.size(); ++i)
    {
      auto k = actIdx_[i];
      switch (activationStatus_[static_cast<size_t>(k)])
      {
      case Activation::Lower: if (lambda[k] > lmax) { lmax = lambda[k]; kmax = k; } break;
      case Activation::Upper: if (-lambda[k] > lmax) { lmax = -lambda[k]; kmax = k; } break;
      default: break;
      }
    }

    //std::cout << "deactivate " << lmax << std::endl;
    deactivate(static_cast<size_t>(kmax));
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

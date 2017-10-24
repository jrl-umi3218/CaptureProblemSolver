#include "ProblemMatrices.h"
#include "QuadraticObjective.h"

#include <iostream>

using namespace Eigen;

namespace bms
{
  LeastSquareObjective::LeastSquareObjective(const VectorXd& delta)
    : n_(delta.size())
    , delta_(delta)
    , d_(delta.cwiseInverse())
    , e_(delta.size()+1)
    , qr_(delta.size())
  {
  }

  LeastSquareObjective::Index LeastSquareObjective::size() const
  {
    return n_;
  }

  double LeastSquareObjective::value(const VectorConstRef& x) const
  {
    assert(x.size() == n_);
    double vi = d_[1] * x[1] - (d_[0] + d_[1])*x[0];
    double v = vi*vi;
    for (Index i = 1; i < n_ - 1; ++i)
    {
      vi = d_[i] * x[i - 1] - (d_[i] + d_[i + 1])*x[i] + d_[i + 1] * x[i + 1];
      v += vi*vi;
    }
    return v/2;
  }

  void LeastSquareObjective::applyJToTheLeft(MatrixRef Y, const MatrixConstRef& X) const
  {
    assert(X.rows() == n_);
    assert(Y.rows() == n_ - 1 && Y.cols() == X.cols());

    Y.row(0) = d_[1] * X.row(1) - (d_[0] + d_[1])*X.row(0);
    for (Index i = 1; i < n_ - 1; ++i)
      Y.row(i) = d_[i] * X.row(i - 1) - (d_[i] + d_[i + 1])*X.row(i) + d_[i + 1] * X.row(i + 1);
  }

  void LeastSquareObjective::applyJTransposeToTheLeft(MatrixRef Y, const MatrixConstRef& X) const
  {
    assert(X.rows() == n_ - 1);
    assert(Y.rows() == n_ && Y.cols() == X.cols());

    Y.row(0) = d_[1] * X.row(1) - (d_[0] + d_[1])*X.row(0);
    for (Index i = 1; i < n_ - 2; ++i)
      Y.row(i) = d_[i] * X.row(i - 1) - (d_[i] + d_[i + 1])*X.row(i) + d_[i + 1] * X.row(i + 1);
    Y.row(n_ - 2) = d_[n_ - 2] * X.row(n_ - 3) - (d_[n_ - 2] + d_[n_ - 1])*X.row(n_ - 2);
    Y.row(n_ - 1) = d_[n_ - 1] * X.row(n_ - 2);
  }

  void LeastSquareObjective::qr(MatrixRef R, CondensedOrthogonalMatrix & Q, const std::vector<bool>& act, Index shift) const
  {
    int na = 0;
    for (auto b : act)
      na += b;

    return qr(R, Q, na, act, shift);
  }

  void LeastSquareObjective::qr(MatrixRef R, CondensedOrthogonalMatrix& Q, Index nact, const std::vector<bool>& act, Index shift) const
  {
    assert(act.size() == n_ + 1);
    assert(R.rows() == n_ - 1 && R.cols() == n_ - nact);
    assert(Q.size() == n_ - 1 + shift);
    R.setZero();
    Q.reset(true);

    //permutation management
    Index up = 0;
    const int reduce = -1;

    if (nact == n_) return;

    //number of consecutive activated constraints at the begining
    Index a0 = 0;
    while (act[static_cast<size_t>(a0)]) ++a0;

    //number of consecutive activated constraints at the end
    Index ap = 0;
    while (act[static_cast<size_t>(n_ - ap)]) ++ap;

    StartType startType = StartType::Case2;
    Index startOffset = 0;
    Index addDim = 0;
    if (a0 > 0)
    {
      startType = StartType::Case1;
      startOffset = a0 - 1;
      addDim = 1;
    }

    //cases where all the active constraints are located at the begining and the end
    if (a0 + ap == nact)
    {
      switch (ap)
      {
      case 0:
        if (nact != n_ - 1)
        {
          qrJj(R.topRows(n_ - 1 - startOffset), Q.Q(0), shift + startOffset, a0, n_ - 1, startType, EndType::Case3);
          Q.P().indices().segment(shift, n_ - 1 - startOffset).setLinSpaced(static_cast<int>(shift + startOffset), static_cast<int>(shift + n_ - 2));
        }
        else
        {
          R(0, 0) = d_(n_ - 1);
          Q.P().indices()[0] = static_cast<int>(n_ - 2);
        }
        break;
      case 1:
        qrJj(R.topRows(n_ - 1 - startOffset), Q.Q(0), shift + startOffset, a0, n_ - 1, startType, EndType::Case2);
        Q.P().indices().segment(shift, n_ - 1 - startOffset).setLinSpaced(static_cast<int>(shift + startOffset), static_cast<int>(shift + n_ - 2));
        break;
      default:
        qrJj(R.topRows(n_ + 1 - ap - startOffset), Q.Q(0), shift + startOffset, a0, n_ - ap, startType, EndType::Case1);
        Q.P().indices().segment(shift, n_ + 1 - ap - startOffset).setLinSpaced(static_cast<int>(shift + startOffset), static_cast<int>(shift + n_ - ap));
        break;
      }
      return;
    }

    Index k = a0; //number of constraints visited (initial scan of the end excluded)
    Index i1 = 0; //size of the first block of inactive constraints
    size_t q = 0; //numbers of GivensSequence used
    up = startOffset;
    while (k + i1 <= n_ && !act[static_cast<size_t>(k + i1)]) ++i1;
    qrJj(R.block(startOffset - up, 0, i1 + addDim, i1), Q.Q(q), shift + startOffset, a0, a0 + i1 - 1, startType, EndType::Case4);
    Q.P().indices().segment(shift, i1 + addDim).setLinSpaced(static_cast<int>(shift + startOffset), static_cast<int>(shift + startOffset + i1 + addDim - 1)); //FIXME: we are doing one permutation too many if startType is Case1
    ++q;
    if (startType == StartType::Case1)  up += 1;

    k += i1;
    Index c = i1 - 1;
    Index a1 = 0;
    while (k + a1 <= n_ && act[static_cast<size_t>(k + a1)]) ++a1;
    k += a1;
    up += a1 - 1;

    Index startHessenberg = k - 2 - up;
    Index finalSize;
    while (true)
    {
      Index ii = 0;
      Index ai = 0;

      while (k + ii <= n_ && !act[static_cast<size_t>(k + ii)]) ++ii;
      while (k + ii + ai <= n_ && act[static_cast<size_t>(k + ii + ai)]) ++ai;

      // if we reach the end
      if (k + ii + ai == n_ + 1)
      {
        assert(ai == ap);
        switch (ai)
        {
        case 0: 
          qrJj(R.block(k - 1 - up, c, ii - 1, ii), Q.Q(q), shift + k - 1, n_ + 1 - ii, n_ - 1, StartType::Case3, EndType::Case3);
          Q.P().indices().segment(shift + k - 1 - up, ii - 1).setLinSpaced(static_cast<int>(shift + k - 1), static_cast<int>(shift + k - 3 + ii));
          finalSize = k - up + ii - 2;
          break;
        case 1: 
          qrJj(R.block(k - 1 - up, c, ii, ii), Q.Q(q), shift + k - 1, n_ - ii, n_ - 1, StartType::Case3, EndType::Case2);
          Q.P().indices().segment(shift + k - 1 - up, ii).setLinSpaced(static_cast<int>(shift + k - 1), static_cast<int>(shift + k - 2 + ii));
          finalSize = k - up + ii - 1;
          break;
        default: 
          qrJj(R.block(k - 1 - up, c, ii + 1, ii), Q.Q(q), shift + k - 1, n_ + 1 - ii - ai, n_ - ai, StartType::Case3, EndType::Case1);
          Q.P().indices().segment(shift + k - 1 - up, ii + 1).setLinSpaced(static_cast<int>(shift + k - 1), static_cast<int>(shift + k - 1 + ii));
          finalSize = k - up + ii - 1;
          break;
        }
        ++q;
        break;
      }
      else
      {
        qrJj(R.block(k - 1 - up, c, ii + 1, ii + 1), Q.Q(q), shift + k - 1, k, k + ii - 1, StartType::Case3, EndType::Case4);
        Q.P().indices().segment(shift + k - 1 - up, ii + 1).setLinSpaced(static_cast<int>(shift + k - 1), static_cast<int>(shift + k - 1 + ii));
        up += ai;
        ++q;
        c += ii;
        k += ii + ai;
      }
    }

    //at this point, we still have a (partially) Hessenberg matrix, so we need a last pass to get a triangular matrix
    if (finalSize>startHessenberg+1)
      tridiagonalQR(R.block(startHessenberg, startHessenberg, finalSize - startHessenberg, R.cols() - startHessenberg), Q.Qh());
    Q.Qh().extend(static_cast<int>(shift + startHessenberg));
  }

  Eigen::MatrixXd LeastSquareObjective::matrix() const
  {
    return buildJ(d_);
  }

  MatrixXd LeastSquareObjective::projectedMatrix(const std::vector<bool>& act) const
  {
    int na = 0;
    for (auto b : act)
      na += b;

    return projectedMatrix(na, act);
  }

  MatrixXd LeastSquareObjective::projectedMatrix(Index nact, const std::vector<bool>& act) const
  {
    assert(act.size() == n_ + 1);
    MatrixXd JN = MatrixXd::Zero(n_ - 1, n_ - nact);

    if (nact == n_) return JN;

    //number of consecutive activated constraints at the begining
    Index a0 = 0;
    while (act[static_cast<size_t>(a0)]) ++a0;

    //number of consecutive activated constraints at the end
    Index ap = 0;
    while (act[static_cast<size_t>(n_ - ap)]) ++ap;

    StartType startType = StartType::Case2;
    Index startOffset = 0;
    Index addDim = 0;
    if (a0 > 0)
    {
      startType = StartType::Case1;
      startOffset = a0 - 1;
      addDim = 1;
    }

    //cases where all the active constraints are located at the begining and the end
    if (a0 + ap == nact)
    {
      switch (ap)
      {
      case 0: 
        if (nact != n_ - 1) buildJj(JN.bottomRows(n_-1-startOffset), a0, n_ - 1, startType, EndType::Case3);
        else JN(n_ - 2, 0) = d_(n_ - 1);
          break;
      case 1: 
        buildJj(JN.bottomRows(n_-1 - startOffset), a0, n_ - 1, startType, EndType::Case2);
        break;
      default: 
        buildJj(JN.middleRows(startOffset, n_ + 1 - ap - startOffset), a0, n_ - ap, startType, EndType::Case1);
        break;
      }
      return JN;
    }

    Index k = a0;
    Index i1 = 0;
    while (k + i1 <= n_ && !act[static_cast<size_t>(k + i1)]) ++i1;
    buildJj(JN.block(startOffset, 0, i1 + addDim, i1), a0, a0 + i1 - 1, startType, EndType::Case4);

    k += i1;
    Index c = i1 - 1;
    Index a1 = 0;
    while (k + a1 <= n_ && act[static_cast<size_t>(k + a1)]) ++a1;
    k += a1;

    while (true)
    {
      Index ii = 0;
      Index ai = 0;

      while (k + ii <= n_ && !act[static_cast<size_t>(k + ii)]) ++ii;
      while (k + ii + ai <= n_ && act[static_cast<size_t>(k + ii + ai)]) ++ai;

      // if we reach the end
      if (k + ii + ai == n_ + 1)
      {
        assert(ai == ap);
        switch (ai)
        {
        case 0: buildJj(JN.block(k - 1, c, ii - 1, ii), n_ + 1 - ii, n_ - 1, StartType::Case3, EndType::Case3); break;
        case 1: buildJj(JN.block(k - 1, c, ii, ii), n_ - ii, n_ - 1, StartType::Case3, EndType::Case2); break;
        default: buildJj(JN.block(k - 1, c, ii + 1, ii), n_ + 1 - ii - ai, n_ - ai, StartType::Case3, EndType::Case1);  break;
        }
        break;
      }
      else
      {
        buildJj(JN.block(k - 1, c, ii + 1, ii + 1), k, k + ii - 1, StartType::Case3, EndType::Case4);
        c += ii;
        k += ii + ai;
      }
    }
    return JN;
  }

  void LeastSquareObjective::buildJj(MatrixRef Jj, Index dstart, Index dend, StartType startType, EndType endType) const
  {
    Index n = dend - dstart + 1;
    auto e = e_.head(n + 1);
    e.head(n) = d_.segment(dstart, n);
    e(n) = 0;

    Index s, ks;
    switch (startType)
    {
    case StartType::Case1: s = 1; ks = 1; break;
    case StartType::Case2: s = 0; ks = 1; break;
    case StartType::Case3: s = 0; ks = 0; break;
    }

    Index er, eb, ke;
    switch (endType)
    {
    case EndType::Case1: eb = 1; er = 0; ke = n - 2; break;
    case EndType::Case2: eb = 0; er = 0; ke = n - 2; break;
    case EndType::Case3: eb = 0; er = 1; ke = n - 2; break;
    case EndType::Case4: eb = 0; er = 0; ke = n - 1; break;
    }

    if (n == 0)
    {
      assert(Jj.rows() == 0 && Jj.cols() == 1);
      return;
    }
    else
    {
      Index k = ke - ks + 1;
      assert(Jj.rows() == s + k + 1 + eb && Jj.cols() == k + 1 + er);
      Jj.setZero();
      Jj.middleRows(s, k + 1).diagonal<-1>().head(k) = e.segment(ks, k);
      Jj.middleRows(s, k + 1).diagonal().tail(k) = -e.segment(ks, k) - e.segment(ks + 1, k);
      Jj.middleRows(s, k + 1).diagonal<+1>().head(k) = e.segment(ks, k);
    }
    //set S
    switch (startType)
    {
    case StartType::Case1: Jj(0, 0) = e[0]; Jj(1, 0) = -e[0] - e[1]; break;
    case StartType::Case2: Jj(0, 0) = -e[0] - e[1]; break;
    case StartType::Case3: Jj(0, 0) = -e[0]; break;
    }
    //set E
    if (endType == EndType::Case1 || endType == EndType::Case3)
      Jj.bottomRightCorner(1, 1)(0, 0) = e[n - 1];
  }

  void LeastSquareObjective::qrJj(MatrixRef R, GivensSequence& Q, Index extend, Index dstart, Index dend, StartType startType, EndType endType) const
  {
    auto n = dend - dstart + 1;
    if (n == 0)
    {
      Q.clear();
      assert(R.rows() == 0 && R.cols() == 1);
      return;
    }
    auto e = d_.segment(dstart, n);
    switch (startType)
    {
    case StartType::Case1: 
      buildJj(R, dstart, dend, startType, endType);
      tridiagonalQR(R.bottomRows(R.rows() - 1), Q, false);
      Q.extend(1);
      if (endType == EndType::Case1) //last line is zero
        tridiagonalQR(R.topRows(R.rows() - 1), Q, true);
      else
        tridiagonalQR(R, Q, true);
      break;
    case StartType::Case2: 
      buildJj(R, dstart, dend, startType, endType);
      tridiagonalQR(R, Q, true);
      break;
    case StartType::Case3: 
      qr_.compute(R, e, Q, endType); 
      break;
    }

    Q.extend(static_cast<int>(extend));
  }
}
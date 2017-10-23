#include "ProblemMatrices.h"
#include "QuadraticObjective.h"

using namespace Eigen;

namespace bms
{
  LeastSquareObjective::LeastSquareObjective(const VectorXd& delta)
    : n_(delta.size())
    , delta_(delta)
    , d_(delta.cwiseInverse())
    , e_(delta.size()+1)
    , transpositionIndices_(delta.size()-1)
  {
  }

  void LeastSquareObjective::qr(MatrixRef R, CondensedOrthogonalMatrix & Q, const std::vector<bool>& act, Index shift) const
  {
    int na = 0;
    for (auto b : act)
      na += b;

    return qr(R, Q, na, act, shift);
  }

  void LeastSquareObjective::qr(MatrixRef R, CondensedOrthogonalMatrix & Q, Index nact, const std::vector<bool>& act, Index shift) const
  {
    assert(act.size() == n_ + 1);
    assert(R.rows() == n_ - 1 && R.cols() == n_ - nact);
    R.setZero();

    //permutation management
    int up = 0;
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
        if (nact != n_ - 1) buildJj(R.bottomRows(n_ - 1 - startOffset), a0, n_ - 1, startType, EndType::Case3);
        else R(n_ - 2, 0) = d_(n_ - 1);
        break;
      case 1:
        buildJj(R.bottomRows(n_ - 1 - startOffset), a0, n_ - 1, startType, EndType::Case2);
        break;
      default:
        buildJj(R.middleRows(startOffset, n_ + 1 - ap - startOffset), a0, n_ - ap, startType, EndType::Case1);
        break;
      }
      return;
    }

    Index k = a0;
    Index i1 = 0;
    while (k + i1 <= n_ && !act[static_cast<size_t>(k + i1)]) ++i1;
    buildJj(R.block(startOffset, 0, i1 + addDim, i1), a0, a0 + i1 - 1, startType, EndType::Case4);

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
        case 0: buildJj(R.block(k - 1, c, ii - 1, ii), n_ + 1 - ii, n_ - 1, StartType::Case3, EndType::Case3); break;
        case 1: buildJj(R.block(k - 1, c, ii, ii), n_ - ii, n_ - 1, StartType::Case3, EndType::Case2); break;
        default: buildJj(R.block(k - 1, c, ii + 1, ii), n_ + 1 - ii - ai, n_ - ai, StartType::Case3, EndType::Case1);  break;
        }
        break;
      }
      else
      {
        buildJj(R.block(k - 1, c, ii + 1, ii + 1), k, k + ii - 1, StartType::Case3, EndType::Case4);
        c += ii;
        k += ii + ai;
      }
    }
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
    buildJj(R, dstart, dend, startType, endType);
  }
}
#pragma once

#include<iterator>

#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>
#include <Givens.h>
#include <GivensSequence.h>

namespace bms
{
  /** Performs a QR in place of the mxn Hessenberg matrix M. Supposes n>=m.
    *
    * Upon completion, M contains the R matrix and Q is a sequence of Givens rotations.
    * The function returns true if R is full rank. 0 is determined by the absolute
    * value thresh.
    */
  template<typename Derived>
  bool hessenbergQR(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, bool append = false, double thresh = 1e-15)
  {
    auto m = M.rows();
    auto n = M.cols();
    assert(n + 1 >= m);
    if (!append)
      Q.clear();
    Q.reserve(Q.size() + m - 1);
    for (Eigen::DenseIndex i = 0; i < m-1; ++i)
    {
      Q.emplace_back(M, i, i);
      Q.back().applyTo(const_cast<Eigen::MatrixBase<Derived>&>(M).rightCols(n-i));
    }
    if (m<=n)
      return abs(M(m - 1, m - 1)) > thresh;
    else
      return abs(M(n - 1, n - 1)) > thresh;
  }

  /** Performs a QR in place of the mxn tridiagonal matrix M. Supposes n>=m-1.
  *
  * Upon completion, M contains the R matrix and Q is a sequence of Givens rotations.
  * The function returns true if R is full rank. 0 is determined by the absolute
  * value thresh.
  */
  template<typename Derived>
  bool tridiagonalQR(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, bool append = false, double thresh = 1e-15)
  {
    auto m = M.rows();
    auto n = M.cols();
    assert(n >= m-1);
    if (!append)
      Q.clear();
    if (m >= 2 && n >= 2)
    {
      auto p = std::min(m, n) - 2;
      //auto p = m - 2;
      Q.reserve(Q.size() + p + 1);
      for (Eigen::DenseIndex i = 0; i < p; ++i)
      {
        Q.emplace_back(M, i, i);
        Q.back().applyTo(const_cast<Eigen::MatrixBase<Derived>&>(M).block<2, 3>(i, i), condensed_t());
      }
      for (Eigen::DenseIndex i = p; i < m - 1; ++i)
      {
        Q.emplace_back(M, i, i);
        if (i+2 == n)
          Q.back().applyTo(const_cast<Eigen::MatrixBase<Derived>&>(M).block<2, 2>(i, i), condensed_t());
        else if (i+3<=n)
          Q.back().applyTo(const_cast<Eigen::MatrixBase<Derived>&>(M).block<2, 3>(i, i), condensed_t());
        else //i+1==n
          Q.back().applyTo(const_cast<Eigen::MatrixBase<Derived>&>(M).block<2, 1>(i, i), condensed_t());
      }
      
      return abs(M(p + 1, p + 1)) > thresh;
    }
    else
    {
      assert(m != 0 && n != 0);
      if (n == 1)
      {
        if (m==1)
          return abs(M(0, 0))> thresh;
        Q.emplace_back(M, 0, 0);
        Q.back().applyTo(const_cast<Eigen::MatrixBase<Derived>&>(M).block<2, 1>(0, 0), condensed_t());
        return abs(M(0,0))> thresh;
      }
      // for m==1, we do nothing
      return M.lpNorm<Eigen::Infinity>() > thresh;
    }
  }

  /** A class for QR decomposition of matrices with particular form.*/
  class BMS_DLLAPI SpecialQR
  {
  public:
    /** kmax is the maximum value of k, as described above
      *
      */
    SpecialQR(Eigen::DenseIndex kmax);

    /** Perform a QR decomposition for M with the form
    * | -e0    e0                                              |
    * |  e0  -e0-e1   e1                                       |
    * |        e1   -e1-e2  e2                                 |
    * |               ...                                      |
    * |                     e_{k-2}  -e_{k-2}-e_{k-1}  e_{k-1} |
    * |                                  e_{k-1}      -e_{k-1} |
    *
    * or a variant with -e_{k-1}-e_k instead of -e_{k-1} for the last element.
    *
    * The matrices are (k+1)x(k+1).
    *
    * DEPRECATED
    */
    template<typename Derived>
    bool compute(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, bool variant, bool append = false, double thresh = 1e-15) const;

    /** A class that can perform a QR decomposition for matrices with the form
    * | -e0    e0                                              |
    * |  e0  -e0-e1   e1                                       |
    * |        e1   -e1-e2  e2                                 |
    * |               ...                                      |
    * |                     e_{k-2}  -e_{k-2}-e_{k-1}  e_{k-1} |
    * |                                  e_{k-1}          E    |
    *
    * where e depends on a EndType enum value:
    * - Case1: E= [-e_{k-1}; -e_k] (the matrix is extend by one row)
    * - Case2: E = -e_{k-1}-e_k
    * - Case3: E = [-e_{k-1}, -e_k] (the matrix is extend by one column)
    * - Case4: E = -e_{k-1}
    *
    * The matrices are (k+2)x(k+1) for Case1, (k+1)x(k+2) for Case3 and
    * (k+1)x(k+1) for Case2 and Case4
    */
    template<typename Derived>
    bool compute(const Eigen::MatrixBase<Derived>& M, const VectorConstRef& e, GivensSequence& Q, EndType endType, bool append = false, double thresh = 1e-15) const;

  private:
    GivensSequence G_;
    mutable Eigen::VectorXd e_;
    mutable Eigen::VectorXd c1_;
    mutable Eigen::VectorXd c2_;
    mutable Eigen::VectorXd d_;
    mutable Eigen::VectorXd l_;
    const Eigen::VectorXd dprecomp_;
    const Eigen::VectorXd lprecomp_;
  };

  template<typename Derived>
  inline bool SpecialQR::compute(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, bool variant, bool append, double thresh) const
  {
    assert(M.rows() == M.cols());
    assert(M.rows() <= e_.size() && "allocated size is not sufficient");

    auto n = M.cols();

    auto e = e_.head(n);
    auto c1 = c1_.head(n - 1);
    auto c2 = c2_.head(n - 1);
    auto d = dprecomp_.head(n - 1);
    auto l = lprecomp_.head(n - 1);

    double ek = 0;
    if (variant)
      ek = M(n - 1, n - 1) + M(n - 2, n - 1);

    e.head(n - 1) = M.diagonal<-1>();
    e(n - 1) = -ek;
    c1.array() = (d.array() - 1)*e.head(n - 1).array();
    c2.array() = d.array()*e.tail(n-1).array();
    c1.array() *= l.array();
    c2.array() *= l.array();

    const_cast<Eigen::MatrixBase<Derived>&>(M).diagonal<-1>().setZero();
    const_cast<Eigen::MatrixBase<Derived>&>(M).diagonal().head(n - 1) = c1;
    const_cast<Eigen::MatrixBase<Derived>&>(M)(n - 1, n - 1) = ek / sqrt(n);
    const_cast<Eigen::MatrixBase<Derived>&>(M).diagonal<1>() = -c1 - c2;
    const_cast<Eigen::MatrixBase<Derived>&>(M).diagonal<2>() = c2.head(n - 2);

    Q.clear();
    Q.reserve(n - 1);
    std::copy(G_.begin(), G_.begin() + (n - 1), std::back_inserter(Q));

    return variant;
  }

  template<typename Derived>
  inline bool SpecialQR::compute(const Eigen::MatrixBase<Derived>& M, const VectorConstRef& e, GivensSequence & Q, EndType endType, bool append, double thresh) const
  {
    auto n = e.size();
    assert(n > 0);
    assert(n <= c1_.size() && "allocated size is not sufficient");
    assert(endType == EndType::Case1 && M.rows() == n + 1 && M.cols() == n + 0
        || endType == EndType::Case2 && M.rows() == n + 0 && M.cols() == n + 0
        || endType == EndType::Case3 && M.rows() == n + 0 && M.cols() == n + 1
        || endType == EndType::Case4 && M.rows() == n + 1 && M.cols() == n + 1);

    auto& R = const_cast<Eigen::MatrixBase<Derived>&>(M);
    R.setZero();
    
    if (n == 1)
    {
      switch (endType)
      {
      case EndType::Case1: R(0, 0) = -e[0] * sqrt(2);  break;
      case EndType::Case2: R(0, 0) = -e[0];  break;
      case EndType::Case3: R(0, 0) = -e[0], R(0, 1) = e[0]; break;
      case EndType::Case4: R(0, 0) = -e[0] * sqrt(2); R(0, 1) = -R(0, 0); break;
      }
    }
    else
    {
      auto c1 = c1_.head(n);
      auto c2 = c2_.head(n);
      auto d = d_.head(n);
      auto l = l_.head(n);

      d = dprecomp_.head(n);
      l = lprecomp_.head(n);
      if (endType == EndType::Case2 || endType == EndType::Case3)
      {
        d(n - 1) = 0;
        l(n - 1) = 1 / sqrt(n);
      }

      c1.array() = (d.array() - 1)*e.array();
      c2.head(n - 1).array() = d.head(n - 1).array()*e.tail(n - 1).array();
      c2[n - 1] = 0;
      c1.array() *= l.array();
      c2.array() *= l.array();

      switch (endType)
      {
      case EndType::Case1:
      case EndType::Case2:
        R.diagonal() = c1;
        R.diagonal<1>() = -c1.head(n - 1) - c2.head(n - 1);
        R.diagonal<2>() = c2.head(n-2);
        break;
      case EndType::Case3:
        R.diagonal() = c1;
        R.diagonal<1>() = -c1 - c2;
        R.diagonal<2>() = c2.head(n-1);
        break;
      case EndType::Case4:
        R.topRows(n).diagonal() = c1;
        R.topRows(n).diagonal<1>() = -c1 - c2;
        R.topRows(n).diagonal<2>() = c2.head(n - 1);
        break;
      }
    }

    if (!append)
      Q.clear();
    Q.reserve(Q.size() + M.rows()-1);
    std::copy(G_.begin(), G_.begin() + M.rows() - 1, std::back_inserter(Q));

    return endType != EndType::Case4 && std::abs(M(n - 1, n - 1)) > thresh;
  }
}

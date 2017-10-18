#pragma once

#include<iterator>

#include <Eigen/Core>

#include <bms_api.h>
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
  bool hessenbergQR(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, double thresh = 1e-15)
  {
    auto m = M.rows();
    auto n = M.cols();
    assert(n >= m);
    Q.clear();
    Q.reserve(m - 1);
    for (Eigen::DenseIndex i = 0; i < m-1; ++i)
    {
      Q.emplace_back(M, i, i);
      Q.back().applyTo(const_cast<MatrixBase<Derived>&>(M).rightCols(n-i));
    }
    return abs(M(m-1,m-1)) > thresh;
  }

  /** Performs a QR in place of the mxn tridiagonal matrix M. Supposes n>=m.
  *
  * Upon completion, M contains the R matrix and Q is a sequence of Givens rotations.
  * The function returns true if R is full rank. 0 is determined by the absolute
  * value thresh.
  */
  template<typename Derived>
  bool tridiagonalQR(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, double thresh = 1e-15)
  {
    auto m = M.rows();
    auto n = M.cols();
    assert(n >= m);
    Q.clear();
    Q.reserve(m - 1);
    for (Eigen::DenseIndex i = 0; i < m-2; ++i)
    {
      Q.emplace_back(M, i, i);
      Q.back().applyTo(const_cast<MatrixBase<Derived>&>(M).block<2,3>(i,i), condensed_t());
    }
    Q.emplace_back(M, m-2, m-2);
    if (m==n)
      Q.back().applyTo(const_cast<MatrixBase<Derived>&>(M).block<2, 2>(m-2, m-2), condensed_t());
    else
      Q.back().applyTo(const_cast<MatrixBase<Derived>&>(M).block<2, 3>(m-2, m-2), condensed_t());
    return abs(M(m-1, m-1)) > thresh;
  }

  /** A class that can perform a QR decomposition for matrices with the form
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
    */
  class BMS_DLLAPI SpecialQR
  {
  public:
    /** kmax is the maximum value of k, as described above
      *
      */
    SpecialQR(int kmax);

    template<typename Derived>
    bool compute(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, bool variant, double thresh = 1e-15) const;

  private:
    GivensSequence G_;
    mutable Eigen::VectorXd e_;
    mutable Eigen::VectorXd c1_;
    mutable Eigen::VectorXd c2_;
    Eigen::VectorXd d_;
    Eigen::VectorXd l_;
  };

  template<typename Derived>
  inline bool SpecialQR::compute(const Eigen::MatrixBase<Derived>& M, GivensSequence& Q, bool variant, double thresh) const
  {
    assert(M.rows() == M.cols());
    assert(M.rows() <= e_.size() && "allocated size is not sufficient");

    auto n = M.cols();

    auto e = e_.head(n);
    auto c1 = c1_.head(n - 1);
    auto c2 = c2_.head(n - 1);
    auto d = d_.head(n - 1);
    auto l = l_.head(n - 1);

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
}
#pragma once
/* Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
 *
 * This file is part of CPS.
 *
 * CPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CPS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>

#include <cps/cps_api.h>
#include <cps/GivensSequence.h>

#include <iostream>

namespace cps
{
  /** A class representing a product Q1 Q2 ... Qk P Qh where the Qi are 
    * orthogonal matrices and P is a permutation matrix.
    * The matrices Qi are expressed as a product of Givens rotations G1 ... Gp.
    *
    * The goal of this class is twofold:
    * - providing an abstraction for easier manipulation
    * - preallocating memory to avoid memory allocation in critical code
    */
  class CPS_DLLAPI CondensedOrthogonalMatrix
  {
  public:
    /** Create an instance preallocating kmax Givens sequences with pmax Givens
      * rotation each. The represented matrix is nxn
      * If Ptranspose is true, P^T is stored instead of P.
      */
    CondensedOrthogonalMatrix(int n, int kmax, int pmax, bool Ptranspose = false);

    /**Empty version*/
    CondensedOrthogonalMatrix() {}

    /** This methods acts exactly as the operator=(), except that it supposes 
      * (and checks by assert) that other will fit in the memory already allocated
      * to this object.
      */
    CondensedOrthogonalMatrix& copyNoAlloc(const CondensedOrthogonalMatrix& other);

    /** Return the size n of the represented matrix.*/
    Eigen::DenseIndex size() const;
    /** Reset to identity*/
    void reset(bool Ptranspose = false);
    /** Resize the allocated memory. This invalidates the previous values.*/
    void resize(int n, int kmax, int pmax);

    /** Get the sequence Q1 Q2 ... Qk.*/
    std::vector<GivensSequence>& Q();
    /** Get the i+1-th element of the above sequence (0-based indexation).*/
    GivensSequence& Q(size_t i);
    /** Get the Qh matrix*/
    GivensSequence& Qh();
    /** Get the P matrix*/
    Eigen::Transpositions<Eigen::Dynamic>& P();

    /** M = Qh^T P^T Q_k^T Q_{k-1}^T .... Q_1^T M
      *
      * Don't let the const ref on M fool you: it is a (recommended) trick to
      * accept temporaries such as blocks. Internally, the const is cast away.
      *
      * We don't use Eigen::Ref here because this function will be used with
      * many small fixed - size matrices, and we want the compiler to take
      * advantage of that.
      */
    template<typename Derived>
    void applyTo(const Eigen::MatrixBase<Derived>& M) const;

    /** Computes M = M Q_1 Q_2 ... G_n P Qh */
    template<typename Derived>
    void applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const;

    /** Return the corresponding nxn orthogonal matrix
    *
    * Use only for debugging purposes.
    */
    Eigen::MatrixXd matrix();

  private:
    /** Wether P (false) or P^T (true) is stored in transpositions_.*/
    bool ptranspose_;
    /** Size of the represented matrix.*/
    Eigen::DenseIndex n_;
    /** The sequence Q1 Q2 ...  Qk*/
    std::vector<GivensSequence> sequences_;
    /** The matrix Qh*/
    GivensSequence Qh_;
    /** The matrix P*/
    Eigen::Transpositions<Eigen::Dynamic> transpositions_;
    /** Internal vector used to reset transpositions_ to identity.*/
    Eigen::VectorXi identityTransposition_;
  };

  inline CondensedOrthogonalMatrix& CondensedOrthogonalMatrix::copyNoAlloc(const CondensedOrthogonalMatrix& other)
  {
    assert(other.sequences_.size() <= sequences_.size());
    assert(other.n_ == n_);
    ptranspose_ = other.ptranspose_;
    Qh_ = other.Qh_;
    transpositions_ = other.transpositions_;
    //we want to keep the preallocated memory, so we just clear the GivensSequences we don't use
    size_t i;
    for (i = 0; i < other.sequences_.size(); ++i)
    {
      assert(other.sequences_[i].size() <= sequences_[i].capacity());
      sequences_[i] = other.sequences_[i];
    }
    size_t kmax = sequences_.size();
    for (; i < kmax; ++i)
      sequences_[i].clear();

    return *this;
  }

  inline Eigen::DenseIndex CondensedOrthogonalMatrix::size() const
  {
    return n_;
  }

  template<typename Derived>
  inline void CondensedOrthogonalMatrix::applyTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& Q : sequences_)
      Q.applyTo(M);
    if (ptranspose_)
      const_cast<Eigen::MatrixBase<Derived>&>(M) = transpositions_*M;
    else
    {
#if EIGEN_VERSION_AT_LEAST(3, 2, 90)
      /** Work around https://stackoverflow.com/questions/49637702/applying-inverse-transposition-in-eigen */
      const auto & indices = transpositions_.indices();
      auto & M_ = const_cast<Eigen::MatrixBase<Derived>&>(M);
      for(Eigen::Index i = indices.size(); i > 0; --i)
      {
        M_.row(i-1).swap(M_.row(indices(i-1)));
      }
#else
      const_cast<Eigen::MatrixBase<Derived>&>(M) = transpositions_.transpose()*M;
#endif
    }
    Qh_.applyTo(M);
  }

  template<typename Derived>
  inline void CondensedOrthogonalMatrix::applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& Q : sequences_)
      Q.applyOnTheRightTo(M);
    if (ptranspose_)
    {
#if EIGEN_VERSION_AT_LEAST(3, 2, 90)
      /** Work around https://stackoverflow.com/questions/49637702/applying-inverse-transposition-in-eigen */
      const auto & indices = transpositions_.indices();
      auto & M_ = const_cast<Eigen::MatrixBase<Derived>&>(M);
      for(Eigen::Index i = indices.size(); i > 0; --i)
      {
        M_.col(i-1).swap(M_.col(indices(i-1)));
      }
#else
      const_cast<Eigen::MatrixBase<Derived>&>(M) = M*transpositions_.transpose();
#endif
    }
    else
      const_cast<Eigen::MatrixBase<Derived>&>(M) = M*transpositions_;
    Qh_.applyOnTheRightTo(M);
  }

  inline std::vector<GivensSequence>& CondensedOrthogonalMatrix::Q()
  {
    return sequences_;
  }

  inline GivensSequence& CondensedOrthogonalMatrix::Q(size_t i)
  {
    return sequences_[i];
  }

  inline GivensSequence& CondensedOrthogonalMatrix::Qh()
  {
    return Qh_;
  }

  inline Eigen::Transpositions<Eigen::Dynamic>& CondensedOrthogonalMatrix::P()
  {
    return transpositions_;
  }
}

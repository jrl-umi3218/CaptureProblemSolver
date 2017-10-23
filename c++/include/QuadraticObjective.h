#pragma once

#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>
#include <CondensedOrthogonalMatrix.h>

namespace bms
{
  class BMS_DLLAPI LeastSquareObjective
  {
  public:
    using Index = Eigen::DenseIndex;

    LeastSquareObjective(const Eigen::VectorXd& delta);

    void qr(MatrixRef R, CondensedOrthogonalMatrix& Q, const std::vector<bool>& act) const;
    void qr(MatrixRef R, CondensedOrthogonalMatrix& Q, const std::vector<bool>& act, int nact) const;

    Eigen::MatrixXd matrix() const;
    /** Return J*Na. For debug*/
    Eigen::MatrixXd projectedMatrix(const std::vector<bool>& act) const;
    Eigen::MatrixXd projectedMatrix(const std::vector<bool>& act, Index nact) const;

  private:
    /** Given e = d_(dstart:dend), build a matrix with main body of the form
      *  S x
      *  x + *
      *    * + *
      *      ...
      *       * + y
      *         y E
      * where S and E are specified by startType and endType, the [x * * ... * y] 
      * are a sub-vector of e, and the + are such that a line * + * (or x + * or
      * y + *) has the form [e_i  -e_i-e_{i+1}  e_{i+1}] 
      * startType:
      * -1: S = [e_0; -e_0 - e_1],  x = e_1
      * -2: S = -e_0 - e_1,         x = e_1
      * -3: S = -e_0,               x = e_0 
      * stopType:
      * -1: E = [same; e(k)],       y = e_{k-1}
      * -2: E = same,               y = e_{k-1}
      * -3: E = [same, e(k)],       y = e_{k-1}
      * -4: E = same,               y = e_k
      * where same means that the value in place is not changed.
      */
    enum class StartType { Case1, Case2, Case3 };
    enum class EndType { Case1, Case2, Case3, Case4 };
    void buildJj(MatrixRef Jj, Index dstart, Index dend, StartType startType, EndType endType) const;

    Eigen::DenseIndex n_;
    Eigen::VectorXd delta_;
    Eigen::VectorXd d_;     // 1/delta

    //computation data
    mutable Eigen::VectorXd e_;
  };
}
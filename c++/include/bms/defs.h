#pragma once

#include <Eigen/Core>

namespace bms
{
  typedef Eigen::Ref<const Eigen::MatrixXd> MatrixConstRef;
  typedef Eigen::Ref<Eigen::MatrixXd>       MatrixRef;
  typedef Eigen::Ref<const Eigen::VectorXd> VectorConstRef;
  typedef Eigen::Ref<Eigen::VectorXd>       VectorRef;

  /** Enums for specifying the shape of the top-left and bottom-right corners
    * of some matrices. 
    * See LeastSquareObjective::buildJj and the class SpecialQR
    */
  enum class StartType { Case1, Case2, Case3 };
  enum class EndType { Case1, Case2, Case3, Case4 };
}
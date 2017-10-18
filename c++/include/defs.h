#pragma once

#include <Eigen/Core>

namespace bms
{
  typedef Eigen::Ref<const Eigen::MatrixXd> MatrixConstRef;
  typedef Eigen::Ref<Eigen::MatrixXd>       MatrixRef;
  typedef Eigen::Ref<const Eigen::VectorXd> VectorConstRef;
  typedef Eigen::Ref<Eigen::VectorXd>       VectorRef;
}
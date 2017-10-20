#pragma once

#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>
#include <GivensSequence.h>

namespace bms
{
  class BMS_DLLAPI LeastSquareObjective
  {
  public:
    LeastSquareObjective(const Eigen::VectorXd& delta) {}

    //void qr(MatrixRef R, ) const

  private:
    Eigen::VectorXd delta_;
    Eigen::VectorXd d_;     // 1/sqrt(delta)
  };
}
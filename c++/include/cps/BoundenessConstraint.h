#pragma once

#include <Eigen/Core>

#include <cps/cps_api.h>
#include <cps/defs.h>

namespace cps
{
  /** The function sum(delta_j/(sqrt(x(j+1))+sqrt(x(j))) - alpha sqrt(x(n)) - b*/
  class CPS_DLLAPI BoundenessConstraint
  {
  public:
    BoundenessConstraint(const Eigen::VectorXd& delta, double alpha, double b);

    void compute(double& val, const VectorConstRef& x) const;
    void compute(double& val, VectorRef grad, const VectorConstRef& x) const;

    void setAlpha(double alpha);
    void setb(double b);

  private:
    Eigen::DenseIndex n_;
    double alpha_;
    double b_;
    Eigen::VectorXd delta_;

    //intermediate data
    mutable Eigen::VectorXd y_;
  };
}

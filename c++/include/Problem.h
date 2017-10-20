#pragma once

#include <string>

#include <Eigen/Core>

#include <bms_api.h>
#include <BoundenessConstraint.h>
#include <LinearConstraints.h>
#include <QuadraticObjective.h>


namespace bms
{
  struct BMS_DLLAPI RawProblem
  {
    void read(const std::string& filepath);

    double g;
    double lmin;
    double lmax;
    Eigen::VectorXd delta;
    double wi_min;
    double wi_max;
    double zi;
    double dzi;
    double zf;

    //solution (optional)
    Eigen::VectorXd Phi_;
  };


  class BMS_DLLAPI Problem
  {
  public:
    Problem(const RawProblem& pb);

    const LeastSquareObjective& objective() const;
    const BoundenessConstraint& nonLinearConstraint() const;
    const LinearConstraints& linearConstraints() const;

  private:
    LeastSquareObjective lso_;
    BoundenessConstraint bc_;
    LinearConstraints lc_;

    RawProblem raw_;
  };
}
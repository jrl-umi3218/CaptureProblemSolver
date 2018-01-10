#pragma once

#include <string>

#include <Eigen/Core>

#include <bms/bms_api.h>
#include <bms/BoundenessConstraint.h>
#include <bms/LinearConstraints.h>
#include <bms/QuadraticObjective.h>


namespace bms
{
  struct BMS_DLLAPI RawProblem
  {
    void read(const std::string& filepath);

    double g = 9.80665; // ISO 80000-3
    double lambda_min;
    double lambda_max;
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
    LeastSquareObjective& objective();
    const BoundenessConstraint& nonLinearConstraint() const;
    BoundenessConstraint& nonLinearConstraint();
    const LinearConstraints& linearConstraints() const;
    LinearConstraints& linearConstraints();
    Eigen::VectorXd::Index size() const;

    void set_zf(double zf);
    void set_zi(double zi);
    void set_dzi(double dzi);
    void set_lambda_min(double lambda_min);
    void set_lambda_max(double lambda_max);
    void set_lambdas(double lambda_min, double lambda_max);
    void set_wi_min(double wi_min);
    void set_wi_max(double wi_max);
    void set_wi(double wi_min, double wi_max);

  private:
    void computeAndSetBounds0();
    void computeAndSetZonotopeBounds();
    void computeAndSetBoundsN();
    void computeAndSetAlpha();
    void computeAndSetb();

    LeastSquareObjective lso_;
    LinearConstraints lc_;
    BoundenessConstraint bc_;

    RawProblem raw_;
  };
}

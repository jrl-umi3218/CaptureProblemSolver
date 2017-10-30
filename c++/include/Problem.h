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
    LeastSquareObjective& objective();
    const BoundenessConstraint& nonLinearConstraint() const;
    BoundenessConstraint& nonLinearConstraint();
    const LinearConstraints& linearConstraints() const;
    LinearConstraints& linearConstraints();

    void set_zf(double zf);
    void set_zi(double zi);
    void set_dzi(double dzi);
    void set_lambda_min(double lmin);
    void set_lambda_max(double lmax);
    void set_lambdas(double lmin, double lmax);
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
    BoundenessConstraint bc_;
    LinearConstraints lc_;

    RawProblem raw_;
  };
}
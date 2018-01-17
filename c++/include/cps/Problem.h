#pragma once

#include <string>

#include <Eigen/Core>

#include <cps/cps_api.h>
#include <cps/BoundenessConstraint.h>
#include <cps/LinearConstraints.h>
#include <cps/QuadraticObjective.h>


namespace cps
{
  struct CPS_DLLAPI RawProblem
  {
    void read(const std::string& filepath);

    double g = 9.80665; // ISO 80000-3
    double lambda_min;
    double lambda_max;
    Eigen::VectorXd delta;
    double init_omega_min;
    double init_omega_max;
    double init_zbar;
    double init_zbar_deriv;
    double target_height;

    //solution (optional)
    Eigen::VectorXd Phi_;
  };


  class CPS_DLLAPI Problem
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

    void set_init_omega(double init_omega_min, double init_omega_max);
    void set_init_omega_max(double init_omega_max);
    void set_init_omega_min(double init_omega_min);
    void set_init_zbar(double init_zbar);
    void set_init_zbar_deriv(double init_zbar_deriv);
    void set_lambda_max(double lambda_max);
    void set_lambda_min(double lambda_min);
    void set_lambdas(double lambda_min, double lambda_max);
    void set_target_height(double target_height);

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

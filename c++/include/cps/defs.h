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

#include <Eigen/Core>

namespace cps
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
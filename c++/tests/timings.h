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


#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <Eigen/QR>

#include <cps/LeastSquare.h>
#include <cps/Problem.h>
#include <cps/ProblemMatrices.h>
#include <cps/QRAlgorithms.h>
#include <cps/QuadraticObjective.h>
#include <cps/SQP.h>

#include "SQPTestCommon.h"


/** n: size of matrices, N: number of tests*/
void QRPerformances(int n, const int N);
void LSPerformance(int n, const int N);
void QRJAPerformance(int n, const int N);
void SQPPerformance(const std::string& filepath, int n, const int N);

/** Test SQP timings
 *
 * \param filepath Paths to problem files
 *
 * \param n Problem sizes (-1: original, otherwise resampled to the given value)
 *
 * \param N Number of test samples
 */
void SQPTimings(const std::vector<std::string>& filepath, const std::vector<int>& n, const int N);
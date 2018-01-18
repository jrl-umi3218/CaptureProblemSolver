#pragma once

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
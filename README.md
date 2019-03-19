CPS: a Capture Problem Solver
=============================

A dedicated solver in C++ for the capture problem presented initially in 
[S. Caron, B. Mallein, "Balance control using both ZMP and COM height variations: A convex boundedness approach", ICRA 2018](https://hal.archives-ouvertes.fr/hal-01590509/document).

The solver itself is presented in
[S. Caron, A. Escande, L. Lanari, B. Mallein, "Capturability-based Analysis, Optimization and Control of 3D Bipedal Walking", under review](https://hal.archives-ouvertes.fr/hal-01689331/document).

This repository contains three folders:
 - the C++ code of the solver
 - a MATLAB code that was used for prototyping
 - a LaTeX technical document, with some of the math behind the solver (not fully complete)

Installation
------------

Compilation has been tested on Linux (gcc/clang) and Windows (Visual Studio).

### Dependencies

To compile you will need the following tools:

 * [Git](https://git-scm.com/)
 * [CMake](https://cmake.org/) >= 2.8.12
 * [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/) (use [pkg-config-lite](https://sourceforge.net/projects/pkgconfiglite/) on Windows)
 * [doxygen](http://www.doxygen.org)
 * A compiler with C++11 support
 
and the following dependencies:
 * [Boost](http://www.boost.org/) >= 1.49 (>= 1.64 for Python bindings)
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2

This repository also uses [jrl-cmakemodules](https://github.com/jrl-umi3218/jrl-cmakemodules) as a submodule.

### Building from source on Linux

Follow the standard CMake build procedure:

```sh
git clone --recursive https://github.com/jrl-umi3218/CaptureProblemSolver
cd CaptureProblemSolver
mkdir build && cd build
cmake [options] ..
make && make install
```

where the main options are:
 * `-DCMAKE_BUILD_TYPE=Release` Build in Release mode
 * `-DCMAKE_INSTALL_PREFIX=some/path/to/install` default is `/usr/local`
 * `-DPYTHON_BINDINGS=ON` Build Python bindings

Use
---

### C++ code

The C++ code is mainly intended to work as a library.

The main class is `cps::CaptureSolver` which provides the solver and is used through its method `solve` to which a `cps::Problem` instance is passed.
A simple example of use can be found in `main.cpp`.

`cps::CaptureSolver` is a thin wrapper around `cps::SQP` where the real work is done. `cps::SQP` is used the same way.

As an alternative, the compilation of the project `main` generates an executable that takes as arguments the paths of files describing capture problem.
Example of such files can be found in `c++\tests\data\`.

### Input file format
A problem can be described by a simple text file such as those found in `c++\tests\data\`.

The file parser looks for (matlab readable, semi-colum terminated) lines with the `=` character in them, and recognize the following fields (in any order):
 * Delta (![](https://latex.codecogs.com/svg.latex?\boldsymbol{\delta}) in the paper)
 * g
 * lambda_min (![](https://latex.codecogs.com/svg.latex?\lambda_{min}))
 * lambda_max (![](https://latex.codecogs.com/svg.latex?\lambda_{max}))
 * omega_i_min (![](https://latex.codecogs.com/svg.latex?\omega_{\mathrm{i},\text{min}}))
 * omega_i_max (![](https://latex.codecogs.com/svg.latex?\omega_{\mathrm{i},\text{max}}))
 * s 
 * z_bar (![](https://latex.codecogs.com/svg.latex?\bar{z}_{\mathrm{i}}))
 * zd_bar (![](https://latex.codecogs.com/svg.latex?\dot{\bar{z}}_{\mathrm{i}}))
 * z_f (![](https://latex.codecogs.com/svg.latex?\bar{z}_{\mathrm{f}}))
 * (optionally) Phi (![](https://latex.codecogs.com/svg.latex?\left[\varphi_0=0,%20\varphi_1,%20\ldots,%20\varphi_n\right])), the solution computed by any other mean.
 
All other lines are ignored.

Thanks
------

- To Vincent Samy for his help with the
  [Boost.Python](http://www.boost.org/doc/libs/1_64_0/libs/python/doc/html/) bindings

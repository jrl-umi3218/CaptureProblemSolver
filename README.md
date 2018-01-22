CPS: a Capture Problem Solver
=============================

A dedicated solver in C++ for the capture problem presented initially in 
[S. Caron, B. Mallein, "Balance control using both ZMP and COM height variations: A convex boundedness approach", ICRA 2018](https://hal.archives-ouvertes.fr/hal-01590509/document).

The solver itself is presented in
[S. Caron, A. Escande, B. Mallein, L. Lanari, "Capturability-based Analysis, Optimization and Control of 3D Bipedal Walking", under review](https://hal.archives-ouvertes.fr/hal-01689331/document).

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
 * [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
 * [doxygen](http://www.doxygen.org)
 * [Boost](http://www.boost.org/) >= 1.49 (>= 1.64 for Python bindings)
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
 * A compiler with C++11 support

Thanks
------

- To Vincent Samy for his help with the
  [Boost.Python](http://www.boost.org/doc/libs/1_64_0/libs/python/doc/html/) bindings

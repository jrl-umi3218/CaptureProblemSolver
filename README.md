Capture Problem Solver
===

This project propose a dedicated solver in C++ for the capture problem presented initially in 
[S. Caron, B. Mallein "Balance control using both ZMP and COM height variations: A convex boundedness approach", ICRA 2018](https://hal.archives-ouvertes.fr/hal-01590509v3/document)

The solver itself is presented in
S. Caron, A. Escande, B. Mallein, L. Lanari, "Capturability-based Analysis, Optimization and Control of 3D Bipedal Walking", under review

The repository contains three folders
 - the c++ code of the solver
 - a matlab code that was used for prototyping (and currently lacks of comments)
 - a latex technical documentation, with some of the math behind the solver

Installation
---
## Dependencies

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 2.8
 * [pkg-config]()
 * [doxygen]()
 * [Boost]() >= 1.49
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
 * [A compiler with C++11 support]()


Thanks
---
- To Vincent Samy for his help with the Boost.Python bindings

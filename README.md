# Nosh

[![Build Status](https://travis-ci.org/nschloe/Nosh.svg?branch=master)](https://travis-ci.org/nschloe/nosh)
[![Coverity Scan Build](https://scan.coverity.com/projects/1659/badge.svg)](https://scan.coverity.com/projects/1659)


This is Nosh,
a free and open-source implementation of numerical solutions methods
for nonlinear Schr\"odinger equations of the form $0=(\K+V+g|\psi|^2)\psi$.

Numerical parameter continuation is the main tool for solving this nonlinear partial differential
equation: Given a number of user-specified system parameters
and an initial guess for one parameter setting, Nosh will find a continuous curve
of solutions for a changing parameter value.
The discretization is a mixed volume-tetrahedral formulation and can handle arbitrarily-shaped domains.
The Jacobian system in the Newton process are solved using a preconditioned MINRES method,
providing exceptional computational efficiency that allows for computations with a large
number of unknowns such as appearing in three-dimensional systems.
Being based on the Trilinos-toolkit, Nosh runs efficiently
in highly parallel high-performance environments.

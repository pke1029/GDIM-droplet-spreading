# GDIM-droplet-spreading
Matlab code for solving the Geometric Diffuse Interface Method (GDIM) thin film equation in one spatial dimension. REF: [placeholder]

## Introduction
The classical thin film equation in 1D is given by

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_th%3D-%5Cpartial_x%5Bh%5E3%5Cpartial_%7Bx%7D%28%5Cpartial_%7Bxx%7Dh%20-%20%5CPi%28h%29%29%5D)

In this work, we proposed the regularisation 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ![equation](https://latex.codecogs.com/gif.latex?%5Cpartial_th%3D-%5Cpartial_x%5Bh%5Cbar%7Bh%7D%5E2%5Cpartial_%7Bx%7D%28%5Cpartial_%7Bxx%7D%5Cbar%7Bh%7D%20-%20%5CPi%28h%2C%5Cbar%7Bh%7D%29%29%5D)

with a new potential

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ![equation](https://latex.codecogs.com/gif.latex?%5CPi%28h%2C%5Cbar%7Bh%7D%29%3D-2%5Cchi%5Cfrac%7B%5C%7Ch%5C%7C_1%5E2%7D%7B%5Clangle%20h%2C%5Cbar%7Bh%7D%5Crangle%5E2%7D%5Cbar%7Bh%7D)

For derivations, please see ref [placeholder]. This repository presents two numerical methods for solving the regularised thin film equation, a fully implicit finite difference scheme and a mesh-free particle method. We also implemented a fast summation algorithm for the particle method that reduces the computational complexity of the particle method from `O(N^2)` to `O(N)`, where `N` is the number of particles.

## Usage
The directory `complete_wetting/` contains cases where ![equation](https://latex.codecogs.com/gif.latex?%5Cchi%3D0) (i.e. no potential energy) and `partial_wetting/` contains cases with ![equation](https://latex.codecogs.com/gif.latex?%5Cchi%3E0). We recommend starting with `complete_wetting/` to get a feeling for the workings of the numerical methods. Both folder contain the following solvers:
- `solve_particle.m` solves the GDMI TFE with the particle method using direct summation.
- `solve_sparse.m` is the same as `solve_particle.m` but drop particles with zero weight since they do not contribute to the solution.
- `solve_fast.m` solves the GDMI TFE with the particle method with the fast summation algorithm. 
- `solve_newton.m` solves the GDMI TFE using a fully implicit finite difference method with Newton's optimisation.

The folders contain other utility functions, namely `clamp.m`, `my_diag.m`, and `my_centered_array.m` to enhance redability of the code. Dependencies are noted within the comment section of every solvers.

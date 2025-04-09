# Numerical Linear Algebra Library & Examples

This repository contains a custom Fortran 90 module `LinAl.f90` implementing key numerical linear algebra algorithms. It also includes several example programs demonstrating how to apply these algorithms to solve real-world and mathematical problems.

## 📦 Contents

### 🔧 Core Module
- `LinAl.f90`: Implements foundational numerical algorithms:
    - Householder QR factorization
    - LU decomposition with partial pivoting
    - Gauss-Jacobi, Gauss-Seidel, and Conjugate Gradient iterative solvers

### 📁 Examples

| Filename                                     | Description |
|---------------------------------------------|-------------|
| `svd_image_compression.f90`                 | Compresses a grayscale image using full SVD and reconstructs low-rank approximations. |
| `iterative_solvers_comparison.f90`          | Solves a linear system using Gauss-Jacobi, Gauss-Seidel, and Conjugate Gradient. |
| `inverse_iteration_solver.f90`              | Computes the eigenvector of a matrix using inverse iteration for a given eigenvalue. |
| `qr_algorithm_eigenvalues.f90`              | Estimates eigenvalues using the QR algorithm with and without shifts. |
| `manual_householder_tridiagonalization.f90` | Applies two steps of Householder tridiagonalization to a small symmetric matrix. |
| `vandermonde_qr_least_squares.f90`          | Performs polynomial least-squares fitting using a Vandermonde matrix and Householder QR. |

---

## 🛠️ Build Instructions

Make sure you have a Fortran compiler (e.g., `gfortran`) and LAPACK installed.

To compile a sample file:

``bash
gfortran -O2 -o vandermonde vandermonde_qr_least_squares.f90 LinAl.f90
./vandermonde

# LinAl — Fortran Linear Algebra Module

This module implements a collection of fundamental linear algebra operations in Fortran, including:

- Hessenberg transformation
- Householder QR decomposition
- Cholesky decomposition and solver
- LU decomposition and solver
- Gaussian elimination and back substitution
- Matrix reading from files

## 📦 Structure

- `LinAl.f90`: The full Fortran module

## 💡 Usage

Import the module and call the subroutines:

```fortran
use LinAl
call Householder(m, n, Q, R)

Program iterative_solvers_comparison
    use Linal, only: GaussJacobi, GaussSeidel, ConjGrad
    !the program prompts the user for a the diagonal value for the matrix A
    !in Ax = b and uses one of Gauss-jacobi, Gauss-Seidel and Conjugate Gradient
    !to find the solution x.
    implicit none
    integer :: m, choice, i
    real :: D
    integer :: di
    real, allocatable :: A(:,:), b(:), x(:)
    real :: accuracy
    accuracy = 10.0**(-5)

    m = 100 ! size of the matrix
    allocate(A(m,m), b(m), x(m))
    !prompt for the value D
    print *, "Enter the value for diagonal elements (D):"
    read *, D
    !populate the values of A
    A = 1.0
    do i = 1, m
        A(i,i) = i
        b(i) = i
    end do
    x = 0.0
    !prompts for choice of alogrithm
    print *, "Choose algorithm: 1 for Gauss-Jacobi, 2 for Gauss-Seidel, 3 for Conjugate Gradient"
    read *, choice
    !conditions on choice to call algorithm from LinAl module.
    di = INT(D)
    if (choice == 1) then
        call GaussJacobi(A, b, x, m, di, accuracy)
    else if (choice == 2) then
        call GaussSeidel(A, b, x, m, di, accuracy)
    else if (choice == 3) then
        call ConjGrad(A, b, x, m, accuracy)
    else
        print *, "Invalid choice"
    end if
    !prints the solution
    print *, "The solution vector x is:"
    print *, x
    deallocate(A, b, x)
End Program iterative_solvers_comparison

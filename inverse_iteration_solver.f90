Program inverse_iteration_solver
    use LinAl, only: LUDecomposition, LUBackSubstitution
    implicit none
    integer, parameter :: n = 4
    real, dimension(n,n) :: A, B, Id
    real :: r(n), y(n), x(n), mu(n)
    integer ::  i, j, k
    integer, dimension(n) :: s
    logical :: f

    !eigenvalues
    mu = [-8.0286, 7.9329, 5.6689, -1.5732]
    
    !Choose the eigenvalue you want to find the eigenvector for
    !by j = index of mu
    j=3

    !Matrix A
    A = reshape([2,1,3,4,1,-3,1,5,3,1,6,-2,4,5,-2,-1],[n,n])
    Id = 0.0
    do i = 1, n
        Id(i,i) = 1.0
    end do
    !inverse iteration
    B = A - (mu(j) * Id)

    print*, "For eigenvalue mu =", mu(j)
    f = .false.
    !Use LU decompositon to find the inverse of B by solving for the identity matrix
    call LUDecomposition(B,n,f,s)
    call LUBackSubstitution(B, Id, n, n, s)

    !initial guess of x
    x = [1.0,1.0,1.0,1.0]

    !initialize r for the loop
    r = 1.0
    do while (sqrt(sum(r**2)) > 0.0000001)
        y = matmul(Id,x) / sqrt(sum((matmul(Id,x))**2))
        r = y - x
        x = y
    end do

    print*, "The eigenvector is:"
    print*, (x(i), i = 1, n)

    
End Program inverse_iteration_solver

Program qr_algorithm_eigenvalues

    use LinAl, only: Householder
    implicit none
    integer, parameter :: m = 3
    real, dimension(m,m) :: A, R, Q, Id, As
    real :: mu
    integer :: i, j

    A = reshape([3,1,0,1,2,1,0,1,1],[m,m])

    print*, 'Without Shift:'
    print*, 'Error convergence:'
    R = A
    !do QR algorithm without shift
    !use householder to decompose A into R and Q 
    call Householder(m,m,Q,R)
    do while (sqrt(sum((A-matmul(R,Q))**2))>0.1)
        A = matmul(R,Q)
        R = A
        call Householder(m,m,Q,R)
        print*, sqrt(sum((A-matmul(R,Q))**2))
    end do

    !returns diagonal matrix of eigenvalues of A
    print*, "Diagonal matrix:"
    do i = 1, m
        print*, (A(i,j), j = 1, m)
    end do

    Id = 0.0
    do i = 1, m
        Id(i,i) = 1.0
    end do

    A = reshape([3,1,0,1,2,1,0,1,1],[m,m])
    print*, "With Shift:"
    print*, "Error convergence:"
    
    do while (sqrt(sum((A-matmul(R,Q))**2))>0.5)
        mu = A(m,m)
        R = A - (mu * Id)
        call Householder(m,m,Q,R)
        A = matmul(R,Q) + (mu * Id)
        print*, sqrt(sum((A-matmul(R,Q))**2))
    end do

    print*, "Diagonal matrix:"
    do i = 1, m
        print*, (A(i,j), j = 1, m)
    end do

End Program qr_algorithm_eigenvalues

Program vandermonde_qr_least_squares
    use LinAl, only: Householder
    implicit none
    real, dimension(:,:), allocatable :: A,R,Q,QR,QTQ,Id
    real, dimension(:), allocatable :: y,d,x, xi, fitted_y
    integer :: m,n,i,j
    real :: fnorm

    !assign parameters
    m = 21
    n = 6
    !read atkinson.dat file
    allocate (xi(m), y(m))
    open (unit=10, file='atkinson.dat',status='old',action='read')
    do i = 1, 21
        read (10, *) xi(i), y(i)
    end do
    close (10)
    !create Vandermonde matrix
    allocate (A(m, n))
    do i = 1, m
        do j = 1, n
            A(i,j) = xi(i)**(j-1)
        end do
    end do
    !perform Householder QR factorization
    allocate (Q(m, m), R(m, n)) 
    R = A !store A in R to be decomposed to store copy of A
    !call Householder with the inputs: matrix R with A elements,
    !empty matrix Q and the size parameters
    call Householder(m,n,Q,R)
    !it returns R and Q matrices post decomposition
    
    allocate (QR(m, n))
    QR = matmul(Q,R)
    print*, 'A - QR ='
    do i=1,m
        print*, (A(i,j) - QR(i,j), j=1,n)
    end do
    fnorm = sqrt(sum((A-QR)**2))
    print*, 'Frobenius norm of A - QR =', fnorm
    allocate (QTQ(m,m), Id(m,m))
    QTQ = matmul(transpose(Q),Q)
    Id = 0.0; forall(i=1:m) Id(i,i)=1.0
    print*, 'Q^TQ - I ='
    do i = 1, m
        print*, (QTQ(i,j) - Id(i,j), j=1,m)
    end do
    fnorm = sqrt(sum((QTQ-Id)**2))
    print*, 'Frobenius norm of Q^TQ - I =', fnorm
    !solve the least squares equations projected onto the span of A
    allocate (d(m))
    d = matmul(transpose(Q),y)
    allocate(x(n))
    x(n) = d(n)/R(n,n)
    !back substitution Rx=Qty
    do i=n-1, 1, -1
        x(i)=d(i)
        do j=i+1,n
            x(i)=x(i)-R(i,j)*x(j)
        end do
        x(i)=x(i)/R(i,i)
    end do
    print*,'solution vector x:'
    do i=1,n
        print*,x(i)
    end do
    fnorm = sqrt(sum((matmul(A,x)-y)**2))
    print*, 'Frobenius norm of Ax-y =', fnorm

    allocate (fitted_y(m))
    fitted_y = matmul(A,x)
    open(2, file='fitted_curve.dat', status='replace', action='write')
    do i = 1, m
        write(2, *) x(i), fitted_y(i)
    end do
    close(2)
End Program vandermonde_qr_least_squares

module LinAl
  implicit none

contains
  !********************************************************************************
  !ConjGrad takes in matrix A, vector b, the size of them and the desired accuracy
  !it outputs the solution x
  !********************************************************************************
  subroutine ConjGrad(A, b, x, n, accuracy)
    implicit none
    integer, intent(in) :: n
    integer :: i, j, l
    real, dimension(n), intent(out) :: x
    real, dimension(n) :: b, r, p, y
    real, dimension(n,n), intent(in):: A
    real, intent(in) :: accuracy
    real :: alpha, beta, E, Enew
    logical :: symmetric

    !check that A is symmetric
    symmetric = .true.
    do i = 1, n
      do j = i + 1, n
        if (A(i,j) /= A(j,i)) then
          symmetric = .false.
          print *, "Matrix is not symmetric"
          exit
        end if
      end do
    end do

    !initialize x, r, p
    x = 1.0
    r = b - matmul(A,x)
    p = r
    E = sum(r**2)
    l = 0
    do while (sqrt(E) > accuracy) !condition on the desired accuracy
      y = matmul(A,p)
      alpha = E / dot_product(p, y)
      x = x + alpha * p
      r = r - alpha * y
      Enew = sum(r**2)
      beta = Enew / E
      p = r + beta * p
      E = Enew
      l = l + 1 !to track the number of iterations
    end do
    print *, "Number of iterations:", l
  end subroutine ConjGrad
  !**************************************************************************
  !the GaussJacobi subroutine takes the matrix A, the vector b, their size
  !and the desire accuracy
  !it outputs the solutions
  !**************************************************************************
  subroutine GaussJacobi(A, b, x, m, di, accuracy)
    implicit none
    character(len=30) :: filename
    integer, intent(in) :: m
    real, dimension(m,m), intent(in) :: A
    real, dimension(m), intent(in) :: b
    real, dimension(m,m) :: R
    real, dimension(m,m) :: Dinv
    real, dimension(m), intent(out) :: x
    real, intent(in) :: accuracy
    integer :: i, j, l
    integer, intent(in) :: di
    

    Dinv = 0.0
    R = 0.0
    !populate the inverse of diagonal matrix D
    do i = 1, m
      do j = 1, m
        if (i == j) then
          Dinv(i,j) = 1.0 / A(i,j)
        else
          R(i,j) = A(i,j)
        end if
      end do
    end do
    l = 0
    x = 1.0 !make a guess for x
    write(filename, '(A, I0, A)') 'jacobierror_', di, '.txt'
    open(unit=(10+di), file=filename, status='replace', action='write')
    do while ((sqrt(sum((b - matmul(A,x))**2)))>accuracy)
      x = matmul(Dinv,(b - matmul(R,x)))
      write(10+di,*) sqrt(sum((b - matmul(A,x))**2))
      l = l + 1
    end do
    close(10+di)
    print *, "Number of iterations:", l 
  end subroutine GaussJacobi
  !**************************************************************************
  !the GaussSeidel subroutine takes the matrix A, the vector b, their size
  !and the desire accuracy
  !it outputs the solutions
  !**************************************************************************
  subroutine GaussSeidel(A, b, x, m, di, accuracy)
    implicit none
    character(len=30) :: filename
    integer, intent(in) :: m
    real, dimension(m,m), intent(in) :: A
    real, dimension(m), intent(in) :: b
    real, dimension(m), intent(out) :: x
    real, intent(in) :: accuracy
    real :: L(m,m), U(m,m), D(m,m), DL(m,m), DLinv(m,m)
    integer :: i, j, iter
    logical :: f
    integer :: s(m)
    integer, intent(in) :: di

    !split A into L and U
    L = 0.0
    U = 0.0
    D = 0.0
    do i = 1, m
      do j = 1, m
        if (i > j) then
          L(i,j) = A(i,j)
        else if (i < j) then
          U(i,j) = A(i,j)
        else
          D(i,j) = A(i,j)
        end if
      end do
    end do

    DLinv = 0.0
    do i = 1, m
      DLinv(i,i) = 1.0
    end do

    DL = D + l
    !use my LU decomposition and back substitution modules for the inverse
    call LUDecomposition(DL,m,f,s)
    call LUBackSubstitution(DL,DLinv,m,m,s)
    x = 1.0
    iter = 0
    write(filename, '(A, I0, A)') 'seidelerror_', di, '.txt'
    open(unit=(10+di), file=filename, status='replace', action='write')
    do while ((sqrt(sum((b - matmul(A,x))**2)))>accuracy)
      write(10+di,*) sqrt(sum((b - matmul(A,x))**2))
      x = matmul(DLinv,(b - matmul(U,x)))
      iter = iter + 1
    end do
    close(10+di)
    print *, "Number of iterations:", iter
  end subroutine GaussSeidel

  !********************************************************
  ! Householder based Hessenberg reduction
  !input: matrix A, number of rows and columns
  !output: matrix A in Hessenberg form
  !********************************************************

  subroutine Hessenberg(m, n, A)
    implicit none
    integer, intent(in) :: m, n
    real, dimension(m, n), intent(inout) :: A
    real, dimension(m) :: vj
    real, dimension(m, m) :: H, Id
    real :: s
    integer :: i, j

    !loop over columns
    do j = 1, n
      !intialize the vector vj with the lower elements of given column of A
      vj = 0.0
      vj(j+1:m) = A(j+1:m, j)
      !compute signed norm
      s = sqrt(sum(vj(j+1:m) ** 2))
      if (A(j+1, j) > 0) then
        vj(j+1) = vj(j+1) + s
      else
        vj(j+1) = vj(j+1) - s
      end if
      !normalize vj
      s = sqrt(sum(vj(j+1:m) ** 2))
      vj(j+1:m) = vj(j+1:m) / s 

      H = 0.0
      Id = 0.0
      do i = 1, m
        Id(i, i) = 1.0
      end do
      !compute Householder matrix Hj
      H = Id - 2.0 * matmul(reshape(vj, (/m,1/)), reshape(vj, (/1,m/)))
      !update A
      A = matmul(H, A)
      A = matmul(A, transpose(H))
    end do
  end subroutine Hessenberg

  !********************************************************
  !Householder based QR decomposition
  !input: matrix A stored in R, number of rows and columns
  !output: matrix R, matrix Q
  !********************************************************

  subroutine Householder(m, n, Q, R)
    implicit none
    integer, intent(in) :: m, n
    real, dimension(m, m), intent(out) :: Q 
    real, dimension(m, n), intent(inout) :: R
    real, dimension(m) :: vj
    real, dimension(m, m) :: H, Id
    real :: s
    integer :: i, j

    !intialize Q as an identity matrix
    Q = 0.0
    do i = 1, m
      Q(i, i) = 1.0
    end do

    !loop over columns
    do j = 1, n
      !intialize the vector vj with the lower elements of given column of A
      vj = 0.0
      vj(j:m) = R(j:m, j)
      !compute signed norm
      s = sqrt(sum(vj(j:m) ** 2))
      if (R(j, j) > 0) then
        vj(j) = vj(j) + s
      else
        vj(j) = vj(j) - s
      end if
      !normalize vj
      s = sqrt(sum(vj(j:m) ** 2))
      vj(j:m) = vj(j:m) / s 

      H = 0.0
      Id = 0.0
      do i = 1, m
        Id(i, i) = 1.0
      end do
      !compute Householder matrix Hj
      H = Id - 2.0 * matmul(reshape(vj, (/m,1/)), reshape(vj, (/1,m/)))
      !update R
      R = matmul(H, R)
      !update Q
      Q = matmul(Q, H)
    end do
  end subroutine Householder


  !********************************************************
  ! Cholesky decomposition
  !input: matrix A, logical flag f indicating whether singular or not and whether the matrix is SPD or not
  !output: matrix A storing L in the lower triangle, flag f
  !********************************************************
 
  subroutine Cholesky(A,f,m)
    implicit none
    integer :: i,j,k
    integer, intent(in)::m
    real, dimension(m,m), intent(INOUT) :: A
    logical, intent(OUT) :: f

    f=.False.

    !loop over columns
    do j = 1, m
      !calculate new diagonal elements
      do k = 1, j-1
        A(j,j) = A(j,j)-A(j,k)*A(j,k)
      end do
      !drop flag for any problems
      if (A(j,j)<=1.0e-10) then
        f=.True.
        return
      end if
      A(j,j) = sqrt(A(j,j))
      !calculate elements below diagonal
      do i = j+1, m
        do k = 1 , j-1
          A(i,j) = A(i,j)-A(i,k)*A(j,k)
        end do
        A(i,j) = A(i,j)/A(j,j)
      end do
    end do
  end subroutine Cholesky

  !********************************************************
  ! Cholesky back substitution
  !input: matrix A, vector b, number of rows of A
  !output: matrix B
  !********************************************************
 
  subroutine CholeskyBackSubstitution(A,b,m)
    implicit none
    integer :: i,j,k
    integer, intent(in)::m
    real :: sum
    real, dimension(m,m), intent(IN) :: A
    real, dimension(m) :: y
    real, dimension(m), intent(INOUT) :: b

    !forward substitution, solving Ly = b
    do i = 1, m
      sum = b(i)
      do j=1,i-1
        sum = sum - A(i,j)*y(j)
      end do
      y(i) = sum/A(i,i)
    end do
    !backward substitution, L^*x = y 
    do i = m, 1, -1
      !stop if entry is zero, singular matrix
      if (A(i,i) < 10e-11) then
        print*,"Singular matrix"
        return
      end if
      do k = i + 1, m
        y(i) = y(i) - A(k,i)*b(k)
      end do
      b(i) = y(i)/A(i,i) !b stores the solution
    end do
  end subroutine CholeskyBackSubstitution

  !********************************************************
  ! Reads a matrix from a file
  !input: name of the file, number of rows and columns
  !output: matrix
  !********************************************************
  subroutine readMat(filename,row,col,mat)

    implicit none
    character(len=*) :: filename
    integer :: i,j,row,col
    real, dimension(row,col), intent(INOUT) :: mat
    
    

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4


    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) row,col

    ! Read matrix
    do i=1,row
       read(10,*) ( mat(i,j), j=1,col )
    enddo

    close(10)
    
  end subroutine readMat

  !*********************************************************
  ! Performs Gaussian elimination
  !input: matrices A and B, number of rows and columns for both and a flag
  !output: matrices A and B, flag
  !*********************************************************
  subroutine gaussianElimination(Amat, Bmat, Arow, Acol, Brow, Bcol, f)
    implicit none
    LOGICAL, intent(OUT) :: f
    integer :: i, j, k, Arow, Acol, Brow, Bcol
    real, dimension(Arow, Acol), intent(INOUT) :: Amat
    real, dimension(Brow, Bcol), intent(INOUT) :: Bmat
    real :: pivot
    real, dimension(Acol) :: temp_row
    real, dimension(Bcol) :: temp_Brow

    f = .False.

    do j = 1, Arow - 1
        ! pivot
        pivot = abs(Amat(j, j))
        k = j

        do i = j + 1, Arow
            if (abs(Amat(i, j)) > pivot) then
                pivot = abs(Amat(i, j))
                k = i ! Update pivot row
            end if
        end do

        ! If pivot row is different, swap rows k and j
        if (k /= j) then
            temp_row = Amat(j, :)
            Amat(j, :) = Amat(k, :)
            Amat(k, :) = temp_row

            temp_Brow = Bmat(j, :)
            Bmat(j, :) = Bmat(k, :)
            Bmat(k, :) = temp_Brow
        end if

        ! Check for singular matrix
        if (Amat(j, j) < 10**(-11)) then
            f = .True.
            return
        end if

        ! Perform row elimination for rows below the pivot row
        do i = j + 1, Arow
            pivot = Amat(i, j) / Amat(j, j)  ! Compute multiplier
            Amat(i, :) = Amat(i, :) - pivot * Amat(j, :)
            Bmat(i, :) = Bmat(i, :) - pivot * Bmat(j, :)
        end do
    end do

  end subroutine gaussianElimination

  !*********************************************************
  ! Performs back substitution
  !input: matrices A and B, number of rows and columns for both
  !output: matrix B as solution X
  !*********************************************************
  subroutine backSubstitution(Amat,Bmat,Arow,Acol,Brow,Bcol)
    implicit none
    integer :: i,j,Arow,Acol,Brow,Bcol
    real, dimension(Arow,Acol), intent(INOUT) :: Amat
    real, dimension(Brow,Bcol), intent(INOUT) :: Bmat
    !Back Substitution
    do i=Arow,1,-1
      Bmat(i,:)=Bmat(i,:)/Amat(i,i)
      do j=i-1,1,-1
        Bmat(j,:)=Bmat(j,:)-Bmat(i,:)*Amat(j,i)
      end do
    end do
  end subroutine backSubstitution

  !*********************************************************
  ! Performs LU decomposition
  !input: matrix A, number of rows and columns, flag
  !output: matrix A, flag and vector s with pivoting information
  !*********************************************************
  
  subroutine LUDecomposition(A,m,f,s)
    implicit none
    integer :: i,j,k,m,q,l
    integer, dimension(m), intent(INOUT) :: s
    real, dimension(m,m), intent(INOUT) :: A
    logical, intent(OUT) :: f
    real :: pivot
    real, dimension(m) :: temp_row
    f=.False.
    do j = 1, m
      s(j) = j
    end do
    ! pivot
    do j = 1, m - 1
      pivot = abs(A(j,j))
      k = j
      do i = j + 1, m
        if (abs(A(i,j))>pivot) then
          pivot = abs(A(i,j))
          k = i
        end if
      end do
      if (k /= j) then
        temp_row = A(j,:)
        A(j,:) = A(k,:)
        A(k,:) = temp_row
        q = s(j)
        s(j) = s(k)
        s(k) = q
      end if
      if (abs(A(j,j)) < 1.0e-11) then
        f=.True.
        return
      end if
      do i = j + 1, m
        A(i,j) = A(i,j)/A(j,j)
        do l = j + 1, m
          A(i,l) = A(i,l) - A(i,j)*A(j,l)
        end do
      end do
    end do
  end subroutine LUDecomposition

  !*********************************************************
  ! Performs LU back substitution
  !input: matrices A and B(containing multiple bi), number of rows and columns for both, vector s with pivoting information
  !output: matrix B as solution X
  !*********************************************************
  subroutine LUBackSubstitution(A, B, m, n, s)
    implicit none
    integer :: i, j, k, m, n
    integer, dimension(m) :: s
    real, dimension(m,m) :: A
    real, dimension(m,n) :: B
    real, dimension(m,n) :: y
    real, dimension(n) :: sum  

    ! Initialize y with Pb
    do j = 1, m
        y(j,:) = B(s(j),:)
    end do

    do j = 1, m - 1
        do i = j + 1, m
            y(i,:) = y(i,:) - A(i,j) * y(j,:)
        end do
    end do

    ! Backward substitution, Ux = y
    do i = m, 1, -1
        ! Stop if entry is zero, singular matrix
        if (abs(A(i,i)) < 1.0e-15) then
            stop "Singular matrix detected"
        end if

        sum = 0.0
        do k = i + 1, m
            sum = sum + A(i,k) * B(k,:) 
        end do

        ! Solve for each column of B (each b vector)
        B(i,:) = (y(i,:) - sum) / A(i,i)
    end do
  end subroutine LUBackSubstitution
  
end module LinAl

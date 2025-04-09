Program manual_householder_tridiagonalization

    implicit none
    real, dimension(4, 4) :: A
    real, dimension(4) :: vj
    real, dimension(4, 4) :: H, Id
    real :: s
    integer :: i, j

    !Assign values of A
    A = reshape([5.0, 4.0, 1.0, 1.0, 4.0, 5.0, 1.0, 1.0, 1.0, 1.0, 4.0, 2.0, 1.0, 1.0, 2.0, 4.0],[4,4])
    !loop over columns
    do j = 1, 2
      !intialize the vector vj with the lower elements of given column of A
      vj = 0.0
      vj(j+1:4) = A(j+1:4, j)
      !compute signed norm
      s = sqrt(sum(vj(j+1:4) ** 2))
      if (A(j+1, j) > 0) then
        vj(j+1) = vj(j+1) + s
      else
        vj(j+1) = vj(j+1) - s
      end if
      !normalize vj
      s = sqrt(sum(vj(j+1:4) ** 2))
      vj(j+1:4) = vj(j+1:4) / s 

      H = 0.0
      Id = 0.0
      do i = 1, 4
        Id(i, i) = 1.0
      end do
      !compute Householder matrix Hj
      H = Id - 2.0 * matmul(reshape(vj, (/4,1/)), reshape(vj, (/1,4/)))
      !update A
      A = matmul(H, A)
      A = matmul(A, transpose(H))
    end do

    print*, 'Matrix A:'
    do i = 1, 4
        print*, (A(i,j), j = 1, 4)
    end do

End Program manual_householder_tridiagonalization

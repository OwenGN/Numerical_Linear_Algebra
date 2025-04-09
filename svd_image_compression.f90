Program svd_image_compression
    use LinAl
    implicit none
    character(len=20) :: filename
    integer, parameter :: m = 1279, n = 1920
    real :: A(m, n), U(m, m), VT(n, n), S(m), Sk(m), As(m,n), E(8)
    real, dimension(:), allocatable :: work
    integer :: i, j, lwork, info, k, karray(8), l
    lwork = -1

    !read in the data
    open(1, file='dog_bw_data.dat', action='read')
    do i = 1, m
        read(1,*) (A(i, j), j=1, n)
    end do
    close(1)
    As = A

    ! Perform SVD
    ! First call to dgesvd to get the optimal work array size LAPACK tip 1
    allocate(work(1))
    call dgesvd('A', 'A', m, n, As, m, S, U, m, VT, n, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dgesvd('A', 'A', m, n, As, m, S, U, m, VT, n, work, lwork, info)
    deallocate(work)

    !reconstruct low-rank approximated images
    karray = ([10, 20, 40, 80, 160, 320, 640, 1279])
    print*, 'Singular values:'
    do i = 1, 10
        print *, i,'singular value =', S(i)
    end do
    do i = 1, 8
        print *, karray(i),'singular value =', S(karray(i))
    end do
    E = 0.0
    do k = 1, 8
        do i = 1, m
            Sk(i) = 0.0
        end do
        do i = 1, karray(k)
            Sk(i) = S(i) !only assign the first k singular values
        end do

        write(filename, '(A, I0, A)') 'Image_appn_', karray(k), '.dat'
        open(unit=(10+k), file=filename, status='replace', action='write')
        do i = 1, m
            do j = 1, n
                As(i, j) = 0.0
                do l = 1, karray(k)
                    As(i, j) = As(i, j) + U(i, l) * Sk(l) * VT(l, j) !only use the first k columns of U and first k rows of VT in U(Sigma)V^T reconstruction
                end do
            end do
            write((10+k),*) (As(i, j), j=1, n)
        end do
        E(k) = sqrt(sum((A-As)**2))/ m*n !Error for each k
        close(10+k)
    end do
    print *, 'Errors: '
    do i = 1, 8
        print *, 'k = ', karray(i), ' Error = ', E(i)
    end do

End Program svd_image_compression

program example
    use params
    !use add_lib
    use KS_econ
    use KS_lib
    implicit none
    double precision :: LAM(M)
    double precision :: G(M)
    double precision :: Z_L(N), Z_U(N)
    integer*8 :: IPROBLEM
    integer*8 :: IPCREATE
    integer :: IERR, NEW_X
    integer :: IPSOLVE, IPADDSTROPTION
    integer :: IPADDNUMOPTION, IPADDINTOPTION
    integer :: IPOPENOUTPUTFILE
    double precision :: F
    integer :: ii,jj, TASK, NZ
    integer :: ACON(N*M), AVAR(N*M)

    double precision :: NORM_ERR


    external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS

! =============================================================================
!                   Declarations and set up of the problem
! =============================================================================

    call linspace(K_GRID,k_min,k_max,k_size)
    K_POL = 2.0 + 0.8*K_GRID
    L = 0.75d0
    call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)
    do jj=1,znum
        L_z(jj) = (z0(znum) - z0(jj))/(z0(znum) - z0(1))
    end do
    print*, L_z

    print*,pr_mat_z
    print*, DOT_PRODUCT(L_Z,pr_mat_z(1,:))





! =============================================================================
!                    Computation of household solution
! =============================================================================

NORM_ERR = 4.0

do while (NORM_ERR>5e-1)
    do ii=1,k_size
        K1_agg = 1.0 + 0.9*K_GRID(ii)
        call data_vectors(DAT,IDAT,K_GRID,K_POL,alpha,delta, sigma, beta, Z,L,K1_agg,pr_mat_z(1,:),k_size,znum,ii)
        X = K_GRID(ii)

        IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,&
                & IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)

        IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
        IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 1)
        IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation', 'limited-memory')
        IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
        call IPFREE(IPROBLEM)


        if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
            K_POL_STORE(ii) = X(1)
        else
            K_POL_STORE(ii) = 0.d0
        end if

    end do
    NORM_ERR = norm2(K_POL-K_POL_STORE)
    print*, NORM_ERR
    K_POL = K_POL_STORE
end do

print*, K_POL_STORE





end program


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

    L = 0.75d0
    call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)
    do jj=1,znum
        L_z(jj) = (z0(znum) - z0(jj))/(z0(znum) - z0(1))
        K_POL(:,jj) =  0.85*K_GRID
    end do

    !print*, DOT_PRODUCT(pr_mat_z(1,:),L_Z)

    !print*, DOT_PRODUCT(pr_mat_z(1,:),UPRIME)

    !call data_vectors(DAT,IDAT,K_GRID,K_POL,alpha,delta, sigma, beta, Z,L,K1_agg,pr_mat_z(2,:),L_Z,k_size,znum,ii)
    !print*, DAT
    !call EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
    !print*, G

    !ii=1
    !jj=1
    !K1_agg = 0.5 + 0.9*K_GRID(10)
    !call data_vectors(DAT,IDAT,K_GRID,K_POL,alpha,delta, sigma, beta, Z,L,K1_agg,pr_mat_z(1,:),L_Z,k_size,znum,ii,jj)


    !print*, DAT
    !print*, IDAT

    !X = K_GRID(ii)/2
    !IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,&
    !            & IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)

    !IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
    !IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 3)
    !IERR = IPADDSTROPTION(IPROBLEM, 'linear_solver', 'mumps')
    !IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation', 'limited-memory')

    !IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
    !call IPFREE(IPROBLEM)
    !print*, "solution at ", K_GRID(ii), "when L is " , L_Z(jj), "is ", X


! =============================================================================
!                    Computation of household solution
! =============================================================================

NORM_ERR = 3.0
do while (NORM_ERR>1e-1)
    do jj=1,znum
    do ii=1,k_size
        K1_agg = K_GRID(18)
        call data_vectors(DAT,IDAT,K_GRID,K_POL,alpha,delta, sigma, beta, Z,L,K1_agg,pr_mat_z(jj,:),L_Z,k_size,znum,ii,jj)
        X = 0.85*K_GRID(ii)

        IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,&
                & IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)

        IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
        IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 1)
        IERR = IPADDSTROPTION(IPROBLEM, 'linear_solver', 'mumps')
        IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation', 'limited-memory')
        IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
        call IPFREE(IPROBLEM)


        if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
            K_POL_STORE(ii,jj) = X(1)
        else
            K_POL_STORE(ii,jj) = K_POL(ii,jj)
        end if

    end do
    print*, "z is" , jj
    print*, K_POL_STORE(:,jj)
    end do
    NORM_ERR = norm2(K_POL(:,1)-K_POL_STORE(:,1)) + norm2(K_POL(:,2)-K_POL_STORE(:,2))
    print*, "norm err", NORM_ERR
    K_POL = K_POL_STORE
end do


do ii=1,k_size
    print*, "policy at", K_GRID(ii), "is ", K_POL_STORE(ii,:)
end do





end program


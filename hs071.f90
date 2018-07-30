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
    integer :: ii,jj,ttt, TASK, NZ
    integer :: ACON(N*M), AVAR(N*M)

    double precision :: NORM_ERR, parameters(6), ZZ(k_size,1)
    integer :: interp, L_state(n_agents), LL, ITER, high_size, low_size

    double precision :: T(n_cheb,k_size), theta(n_cheb,znum)
    double precision :: K_TRY(n_agents), TT(n_cheb), K_SIM(n_sim), switch(n_sim)
    double precision :: rand(n_agents), agg_rand, empl

    double precision, allocatable :: Kagg_low(:,:), Kagg_high(:,:), Kagg_lowlag(:,:), Kagg_highlag(:,:)

    !double precision :: x_prod(reg_dim,reg_dim), x_inv(reg_dim,reg_dim), predict(reg_dim,1)
    !double precision :: x_prodd(reg_dim,reg_dim), x_invv(reg_dim,reg_dim), predictt(reg_dim,1)

    double precision :: beta_high(reg_dim,1), beta_low(reg_dim,1)



    external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS

! =============================================================================
!                   Declarations and set up of the problem
! =============================================================================
    L = 0.5d0
    Zstate(1:Lnum) = 0.5
    Lstat_index(1:Lnum,1) = 1


    parameters = (/alpha,delta, sigma, beta, k_max,k_min/)


    interp = 2  ! if = 1, linear interpolation; if = 2, chebyshev interpolation

    if (interp .eq. 1) then
        ! Cartesian Grid
        call linspace(K_GRID,k_min,k_max,k_size)
    else if (interp .eq. 2) then
        ! Chebyshev Nodes
        do ii=1,k_size
            ZZ(ii,1) = -cos(pi*(2*ii-1)/(2*k_size))
            K_GRID(ii) = (ZZ(ii,1) + 1)*((k_max-k_min)/2) + k_min
        end do
        call TMATRIX(T,ZZ,n_cheb,k_size)
    end if

    call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)
    do jj=1,znum
        L_z(jj) = (z0(znum) - z0(jj))/(z0(znum) - z0(1))
    end do

    call tauchen(Lnum,rhoagg,sigmaagg,nstdevagg,pr_mat_L,z0)





! =============================================================================
!                    Computation of household solution
! =============================================================================

K_agg = 1.1*K_GRID((k_size+3)/2)

do ii = 1,n_agents
    L_STATE(ii) = 0
    K_TRY(ii) = K_GRID((k_size+1)/2)
end do
L_STATE(n_agents/2+1:n_agents) = 1

high_size = 0
low_size = 0

do ttt = 1,n_sim
    call RANDOM_NUMBER(agg_rand)
    if (agg_rand .GE. pr_mat_z(Lstat_index(1,ttt),1)) then
        Lstat_index(2,ttt) = 2
        Zstate(2) = 0.54d0 !high value
        if (ttt .ge. 2) then
            high_size = high_size + 1
        end if
        switch(ttt) = abs(Lstat_index(2,ttt)-Lstat_index(1,ttt))
        K1_agg = exp(0.67418757011682917 + 0.77135610233800866*log(K_agg) -3.2308659286596875E-002*switch(ttt))

    else
        Lstat_index(2,ttt) = 1
        Zstate(2) = 0.46d0 ! low value
        if (ttt .ge. 2) then
            low_size = low_size + 1
        end if
        switch(ttt) = abs(Lstat_index(2,ttt)-Lstat_index(1,ttt))
        K1_agg = exp(0.69214043631004074 + 0.76367097151457131 *log(K_agg) -9.8659257720790094E-005*switch(ttt))

    end if



    print*, switch(ttt)

    NORM_ERR = 3.0
    ITER = 0
    do jj=1,znum
    K_POL(:,jj) =  0.85*K_GRID
    end do

    do while (NORM_ERR>1e-1 .AND. ITER<40)
        ITER = ITER + 1
        do jj = 1,znum
            do ii = 1,n_cheb
                theta(ii,jj) = dot_product(K_POL(:,jj),T(ii,:))/dot_product(T(ii,:),T(ii,:))
            end do
        end do

        do jj=1,znum
        do ii=1,k_size

            call data_vectors(DAT,IDAT,K_GRID,K_POL,ZZ,parameters,znum,k_size,n_cheb, &
                    & L,K_agg,K1_agg,pr_mat_z(jj,:),theta,L_Z,ii,jj,interp,Zstate)
            X = K_POL(ii,jj)

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
                if (interp.eq.1) then
                    K_POL_STORE(ii,jj) = K_POL(ii,jj)
                else if (interp .eq. 2) then
                    K_POL_STORE(ii,jj) = 0.9*K_GRID(ii)
                end if
            end if

        end do
        !print*, "z is" , jj
        !print*, K_POL_STORE(:,jj)
        end do
        NORM_ERR = norm2(K_POL(:,1)-K_POL_STORE(:,1)) + norm2(K_POL(:,2)-K_POL_STORE(:,2))
        print*, "norm err for iteration", ITER , "is", NORM_ERR
        K_POL = 0.5*(K_POL_STORE+K_POL)
    end do

    print*, "************************************************************"
    print*, "iteration number", ttt
    print*, "solutions for", k_size , "values when L is", L_Z
    print*, "************************************************************"
    do ii=1,k_size
        print*, "policy at", K_GRID(ii), "is ", K_POL_STORE(ii,:)
    end do


    do ii = 1,n_agents
        call RANDOM_NUMBER(rand)
        !print *, rand(ii), pr_mat_z(L_STATE(ii)+1,1)
        if (rand(ii) .GE. pr_mat_z(L_STATE(ii)+1,1)) then
            L_STATE(ii) = 1
        else
            L_STATE(ii) = 0
        end if
        call TMATRIX(TT,2*( K_TRY(ii)- k_min)/(k_max - k_min) - 1,n_cheb,1)
        K_TRY(ii) = dot_product(theta(:, L_STATE(ii)),TT)
    end do
    empl = sum(L_STATE)
    print*, "average employment state at time ", ttt,  "is", empl/n_agents
    print*, "aggregate state is", Lstat_index(1,ttt), "today and",  Lstat_index(2,ttt), "tomorrow"
    print*, "when K is ", K_agg, "forecast is ", K1_agg ,"and aggregate policy is " ,sum(K_TRY)/n_agents

    K_SIM(ttt) = sum(K_TRY)/n_agents
    K_agg = K_SIM(ttt)
    !K1_agg = 0.25  + 0.95*K_agg
    Lstat_index(1,ttt+1) = Lstat_index(2,ttt)
    Zstate(1) = Zstate(2)


end do

print*, "there have been", high_size , "high productivity periods"
print*, "and", low_size , "low productivity periods"

allocate( Kagg_low(low_size,1), Kagg_high(high_size,1) , Kagg_lowlag(low_size,reg_dim), Kagg_highlag(high_size,reg_dim) )

low_size = 0
high_size = 0
do ttt = 2,n_sim
    print*, "************************************************************"
    print*, "current state", Lstat_index(1,ttt), "future state", Lstat_index(2,ttt), "capital", K_SIM(ttt)
    print*, "************************************************************"
    if (Lstat_index(2,ttt) .eq. 2) then
        high_size = high_size + 1
        Kagg_high(high_size,1) = log(K_SIM(ttt))
        Kagg_highlag(high_size,1) = 1.0d0
        Kagg_highlag(high_size,2) = log(K_SIM(ttt-1))
        Kagg_highlag(high_size,3) = switch(ttt)
    else
        low_size = low_size + 1
        Kagg_low(low_size,1) = log(K_SIM(ttt))
        Kagg_lowlag(low_size,1) = 1.0d0
        Kagg_lowlag(low_size,2) = log(K_SIM(ttt-1))
        Kagg_lowlag(low_size,3) = switch(ttt)
    end if
end do

print*, "there have been", high_size , "high productivity periods"
print*, "and", low_size , "low productivity periods"


!x_prod = matmul(transpose(Kagg_highlag),Kagg_highlag)
!call inverse(x_prod,x_inv,reg_dim)
!predict =  matmul(transpose(Kagg_highlag),Kagg_high)
!beta_high = matmul(x_inv,predict)
!print*, "x_prod", x_prod
!print*, "x_inv", x_inv
!print*, "predict", predict
!print*, "betas_high", beta_high

call OLS(beta_high,Kagg_highlag,Kagg_high,high_size,reg_dim)
print*, "betas_high", beta_high


print*, "so predicted value of" , K_GRID((k_size)/2), "in good times"
print*, "ends up being", exp(beta_high(1,1) + beta_high(2,1)*log(K_GRID((k_size)/2)))

!x_prodd = matmul(transpose(Kagg_lowlag),Kagg_lowlag)
!call inverse(x_prodd,x_invv,reg_dim)
!predictt =  matmul(transpose(Kagg_lowlag),Kagg_low)
!beta_low = matmul(x_invv,predictt)
!print*, "x_prod", x_prodd
!print*, "x_inv", x_invv
!print*, "predict", predictt
!print*, "betas_low", beta_low

call OLS(beta_low,Kagg_lowlag,Kagg_low,low_size,reg_dim)
print*, "betas_high", beta_low


print*, "and predicted value of" , K_GRID((k_size)/2), "in bad times"
print*, "ends up being", exp(beta_low(1,1) + beta_low(2,1)*log(K_GRID((k_size)/2)))


! =============================================================================
!                              First Simulation
! =============================================================================







end program






        module KS_lib
            use KS_econ
            implicit none
	        integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
	        parameter  (N = 1, M = 1, NELE_JAC = 1, NELE_HESS = 1)
    	    parameter  (IDX_STY = 1 )
	        double precision :: X(N),X_L(N),X_U(N)
	        double precision :: G_L(M),G_U(M)
            double precision :: DAT((znum+2)*k_size+2*znum+ 9 +n_cheb*znum)
            integer :: IDAT(6)


    	    data X   / 10d0 /
	        data X_L /  0d0 /
	        data X_U / 100d0 /

            !  Set bounds for the constraints

            data G_L / 0d0 /
	        data G_U / 0d0 /

        end


! =============================================================================
!                    Computation of data vectors
! =============================================================================
        subroutine data_vectors(DAT,IDAT,K_GRID,K_POL,ZZ,params,znum, &
        & k_size,n_cheb,L,K_agg,K1_agg,pr_vec,theta,L_VEC,ii,jj,itp,Zstate)
            implicit none
            integer, intent(in) :: ii,jj,znum,k_size,n_cheb,itp
            integer :: IDAT(6), kkk
            double precision, intent(in) :: params(*),K_GRID(k_size), K_POL(k_size,znum),ZZ(k_size)
            double precision, intent(in) :: theta(n_cheb,znum),K_agg,Zstate(2)
            double precision ::  K_POLS(k_size*znum), thetas(n_cheb*znum), pr_vec(znum), L_VEC(znum), k_max, k_min
            double precision :: DAT((znum+1)*k_size+2*znum+9), alpha, delta, sigma, beta, Z, L, K1_agg, K_SS
            double precision :: R, W, R1, W1
            double precision ::  T(n_cheb,k_size)

            alpha = params(1)
            delta = params(2)
            sigma = params(3)
            beta = params(4)
            k_max = params(5)
            k_min = params(6)

            IDAT(1) = ii ! position over K grid
            IDAT(2) = k_size
            IDAT(3) = znum
            IDAT(4) = jj ! position over z grid
            IDAT(5) = itp ! interpolation type
            IDAT(6) = n_cheb

            K_SS =  K_agg  ! IDAT(1)

            K_POLS = reshape(K_POL,(/k_size*znum/))
            thetas = reshape(theta,(/n_cheb*znum/))

            call int_rate(R,alpha,Zstate(1), K_SS,L)
            call wage(W,alpha,Zstate(1), K_SS,L)
            call int_rate(R1,alpha,Zstate(2),K1_agg,L)
            call wage(W1,alpha,Zstate(2),K1_agg,L)

            !if (interp_type.eq.2) then
            !    ! Chebyshev Basis Functions
            !    call TMATRIX(T,ZZ,n_cheb,k_size)
                ! Coeff theta
            !end if

            DAT(1:k_size) = K_GRID
            DAT(     k_size+1    : (znum+1)*k_size) = K_POLS
            DAT((znum+1)*k_size+1: (znum+2)*k_size) = ZZ
            DAT((znum+2)*k_size+1) = R
            DAT((znum+2)*k_size+2) = W
            DAT((znum+2)*k_size+3) = R1
            DAT((znum+2)*k_size+4) = W1
            DAT((znum+2)*k_size+5) = delta
            DAT((znum+2)*k_size+6) = sigma
            DAT((znum+2)*k_size+7) = beta
            DAT((znum+2)*k_size+8) = k_min
            DAT((znum+2)*k_size+9) = k_max
            DAT((znum+2)*k_size+9 +1 : (znum+2)*k_size +9  +znum) = pr_vec
            DAT((znum+2)*k_size+9 +1 + znum :  (znum+2)*k_size+9 + 2*znum) = L_VEC
            DAT((znum+2)*k_size+9 +2*znum + 1 : (znum+2)*k_size+9 +2*znum + znum*n_cheb ) = thetas


        end



! =============================================================================

!                    Computation of objective function

! =============================================================================
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      F = 0.0
      IERR = 0
      return
      end
! =============================================================================

!                Computation of gradient of objective function

! =============================================================================
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
      implicit none
      integer i, N, NEW_X
      double precision GRAD(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

       do i=1,N
       GRAD(i) = 0.0
       end do
      IERR = 0
      return
      end
! =============================================================================

!                     Computation of equality constraints

! =============================================================================
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      implicit none
      integer :: N, NEW_X, M
      double precision :: G(M), X(N)
      double precision :: DAT(*)
      integer :: IDAT(*), ZNUM
      integer :: IERR, POS, ZPOS, NUM, i, interp, n_cheb
      double precision :: SIGMA, KAPPA, K_TILDE, K_INT(1)
      double precision, allocatable :: K_GRID(:),K_POL(:,:),K_POLS(:), ZZ(:)
      double precision, allocatable :: PR_VEC(:), L_VEC(:), thetas(:), theta(:,:), TT(:)
      double precision :: R,W,R1,W1, delta, beta, u_prime, UPRIME, k_min, k_max

        !  inputs:
        POS = IDAT(1)
        NUM = IDAT(2)
        ZNUM = IDAT(3)
        ZPOS = IDAT(4)
        interp = IDAT(5)
        n_cheb = IDAT(6)

        allocate(PR_VEC(ZNUM),L_VEC(ZNUM))
        allocate(K_GRID(NUM),K_POLS(NUM*ZNUM),K_POL(NUM,ZNUM),ZZ(NUM))
        allocate(thetas(n_cheb*ZNUM), theta(n_cheb,ZNUM), TT(n_cheb))

        K_GRID = DAT(1:NUM)
        K_POLS = DAT(NUM+1:(ZNUM+1)*NUM)
        ZZ = DAT((ZNUM+1)*NUM+1: (ZNUM+2)*NUM)
        R = DAT((ZNUM+2)*NUM + 1)
        W = DAT((ZNUM+2)*NUM + 2)
        R1 = DAT((ZNUM+2)*NUM + 3)
        W1 = DAT((ZNUM+2)*NUM + 4)
        delta = DAT((ZNUM+2)*NUM + 5)
        sigma = DAT((ZNUM+2)*NUM + 6)
        beta = DAT((ZNUM+2)*NUM + 7)
        k_min = DAT((ZNUM+2)*NUM + 8)
        k_max = DAT((ZNUM+2)*NUM + 9)
        PR_VEC = DAT((ZNUM+2)*NUM+9 + 1 : (ZNUM+2)*NUM+ 9 + ZNUM)
        L_VEC = DAT( (ZNUM+2)*NUM+9 + ZNUM + 1 :  (ZNUM+2)*NUM+ 9 + 2*ZNUM)
        thetas = DAT((znum+2)*NUM+9 +2*ZNUM + 1 : (ZNUM+2)*NUM+9 +2*ZNUM + n_cheb*ZNUM )

        K_POL = reshape(K_POLS,(/NUM,ZNUM/))
        theta = reshape(thetas,(/n_cheb,ZNUM/))


        UPRIME = 0.0
        do i=1,ZNUM
            if (interp.eq.1) then
                call pwl_value_1d ( NUM, K_GRID, K_POL(:,i), 1, X, K_INT )
            else if (interp .eq. 2) then
                !call pwl_value_1d ( NUM, K_GRID, K_POL(:,i), 1, X, K_INT )
                call TMATRIX(TT,2*( X - k_min)/(k_max - k_min) - 1,n_cheb,1)
                K_INT = dot_product(theta(:,i),TT)
            end if
            K_TILDE = K_INT(1)
            UPRIME = UPRIME + PR_VEC(i)*u_prime( L_VEC(i)*W1 + (1-delta+R1)*X - K_TILDE  , sigma)
        end do

        G =  u_prime((1-delta+R)*K_GRID(POS) + W*L_VEC(ZPOS) - X, sigma) -  &
               & beta*(1-delta+R1)*UPRIME

      !u_prime( W1*L_VEC(ZPOS) + (1-delta+R1)*X - K_TILDE  , sigma)
      IERR = 0
      return
      end

! =============================================================================

!                Computation of Jacobian of equality constraints

! =============================================================================
subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A, &
        IDAT, DAT, IERR)
    integer TASK, N, NEW_X, M, NZ
    double precision X(N), A(NZ), xAdj(N), Gx1(M), Gx2(M)
    integer ACON(NZ), AVAR(NZ), I
    double precision DAT(*)
    integer IDAT(*)
    integer IERR
    !
    !     structure of Jacobian:
    !
    integer AVAR1(M*N), ACON1(M*N)
    AVAR1(1) = 1
    ACON1(1) = 1
    !save  AVAR1, ACON1
    !
    if(TASK.eq.0) then
        do I = 1, M*N
            AVAR(I) = AVAR1(I)
            ACON(I) = ACON1(I)
        enddo
    else
        call EV_G(N, X, NEW_X, M, Gx1, IDAT, DAT, IERR)
        h = 1e-4

        do ii = 1,N
            do jj = 1,M
                xAdj = X
                xAdj(ii) = xAdj(ii) + h
                call EV_G(N, xAdj, NEW_X, M, Gx2, IDAT, DAT, IERR)
                A(ii + (jj-1)*N) = (Gx2(jj) -  Gx1(jj))/h
            end do
        end do
    endif
    IERR = 0
    return

end

! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM, &
        NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
    implicit none
    integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i
    double precision X(N), OBJFACT, LAM(M), HESS(NNZH)
    integer IRNH(NNZH), ICNH(NNZH)
    double precision DAT(*)
    integer IDAT(*)
    integer IERR
    !
    !     structure of Hessian:
    !
    integer IRNH1(10), ICNH1(10)
    data  IRNH1 /1, 2, 2, 3, 3, 3, 4, 4, 4, 4/
    data  ICNH1 /1, 1, 2, 1, 2, 3, 1, 2, 3, 4/
    save  IRNH1, ICNH1

    if(TASK.eq.0) then
        do i = 1, 10
            IRNH(i) = IRNH1(i)
            ICNH(i) = ICNH1(i)
        enddo
    else
        do i = 1, 10
            HESS(i) = 0d0
        enddo
        !
        !     objective function
        !
        !         HESS(1) = OBJFACT * 2d0*X(4)
        !         HESS(2) = OBJFACT * X(4)
        !         HESS(4) = OBJFACT * X(4)
        !         HESS(7) = OBJFACT * (2d0*X(1) + X(2) + X(3))
        !         HESS(8) = OBJFACT * X(1)
        !         HESS(9) = OBJFACT * X(1)
        !
        !     first constraint
        !
        !         HESS(2) = HESS(2) + LAM(1) * X(3)*X(4)
        !         HESS(4) = HESS(4) + LAM(1) * X(2)*X(4)
        !         HESS(5) = HESS(5) + LAM(1) * X(1)*X(4)
        !         HESS(7) = HESS(7) + LAM(1) * X(2)*X(3)
        !         HESS(8) = HESS(8) + LAM(1) * X(1)*X(3)
        !         HESS(9) = HESS(9) + LAM(1) * X(1)*X(2)
        !
        !     second constraint
        !
        !         HESS(1) = HESS(1) + LAM(2) * 2d0
        !         HESS(3) = HESS(3) + LAM(2) * 2d0
        !         HESS(6) = HESS(6) + LAM(2) * 2d0
        !         HESS(10)= HESS(10)+ LAM(2) * 2d0
    endif
    IERR = 0
    return
end
!
! =============================================================================
!
!                   Callback method called once per iteration
!
! =============================================================================
!
subroutine ITER_CB(ALG_MODE, ITER_COUNT, OBJVAL, INF_PR, INF_DU, &
        MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, &
        DAT, ISTOP)
    implicit none
    integer ALG_MODE, ITER_COUNT, LS_TRIAL
    double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
    double precision ALPHA_DU, ALPHA_PR
    double precision DAT(*)
    integer IDAT(*)
    integer ISTOP
    !
    !     You can put some output here
    !
    write(*, *) 'Testing callback function in iteration ', ITER_COUNT
    !
    !     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
    !     simple example.
    !
    if (INF_PR.le.1D-04) ISTOP = 1

    return
end






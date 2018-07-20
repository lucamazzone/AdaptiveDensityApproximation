module add_lib
    use KS_econ
    implicit none
	integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
	parameter  (N = 2, M = 2, NELE_JAC = 4, NELE_HESS = 10)
	parameter  (IDX_STY = 1 )
	double precision :: X(N),X_L(N),X_U(N)
	double precision :: G_L(M),G_U(M)
    double precision :: DAT(2)
    integer :: IDAT(1)


	data X   / 1d0, 5d0 /
	data X_L / 0d0, 0d0 /
	data X_U / 100d0, 100d0 /

!     Set bounds for the constraints

	data G_L / 0d0, 0d0 /
	data G_U / 0d0, 0d0 /
end


! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!
!
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
!
!                Computation of gradient of objective function
!
! =============================================================================
!
subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
    implicit none
    integer i, N, NEW_X
    double precision GRAD(N), X(N)
    double precision DAT(*)
    integer IDAT(*)
    integer IERR
    !      GRAD(1) = X(4)*(2d0*X(1)+X(2)+X(3))
    !      GRAD(2) = X(1)*X(4)
    !      GRAD(3) = X(1)*X(4) + 1d0
    !      GRAD(4) = X(1)*(X(1)+X(2)+X(3))

    do i = 1, N
        GRAD(i) = 0.0
    end do
    IERR = 0
    return
end
!
! =============================================================================
!
!                     Computation of equality constraints
!
! =============================================================================
!
subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
    implicit none
    integer N, NEW_X, M
    double precision G(M), X(N)
    double precision DAT(*)
    integer IDAT(*)
    integer IERR
    G(1) = 1.5 * 0.3 * x(1)**(-0.7) * x(2)**(0.7) - x(1) + DAT(1)
    G(2) = 1.5 * 0.3 * x(1)**(0.3) * x(2)**(-0.3) - x(2)
    IERR = 0
    return
end
!
! =============================================================================
!
!                Computation of Jacobian of equality constraints
!
! =============================================================================
!
subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A, &
        IDAT, DAT, IERR)
    integer TASK, N, NEW_X, M, NZ
    double precision X(N), A(NZ)
    integer ACON(NZ), AVAR(NZ), I
    double precision DAT(*)
    integer IDAT(*)
    integer IERR
    !
    !     structure of Jacobian:
    !
    integer AVAR1(4), ACON1(4)
    data  AVAR1 /1, 1, 2, 2/
    data  ACON1 /1, 2, 1, 2/
    save  AVAR1, ACON1
    !
    if(TASK.eq.0) then
        do I = 1, 4
            AVAR(I) = AVAR1(I)
            ACON(I) = ACON1(I)
        enddo
    else
        A(1) = 1.5 * (-0.21) * x(1)**(-1.7) * x(2)**(0.7) - 1
        A(2) = 1.5 * (0.21) * x(1)**(-0.7) * x(2)**(-0.3)
        A(3) = 1.5 * (0.21) * x(1)**(-0.7) * x(2)**(-0.3)
        A(4) = 1.5 * (-0.21) * x(1)**(0.3) * x(2)**(-1.3) - 1
        !        A(5) = 2d0*X(1)
        !        A(6) = 2d0*X(2)
        !        A(7) = 2d0*X(3)
        !        A(8) = 2d0*X(4)
    endif
    IERR = 0
    return
end
!
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


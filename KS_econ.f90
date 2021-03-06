module KS_econ
implicit none
    integer, parameter :: k_size = 41, ni = 1, n_agents = 2000, n_sim = 151
    integer, parameter :: znum = 2, Lnum = 2
    integer, parameter :: n_cheb = k_size-1 !k_size-1
    integer, parameter :: reg_dim = 3
    integer :: Lstat_index(Lnum,n_sim+1)
    double precision :: K_FUN(k_size) , K_GRID(k_size), K_POL(k_size,znum)
    double precision :: K_POL_STORE(k_size,znum)
    double precision, parameter :: k_min = 0.d0, k_max = 22.d0
    double precision, parameter :: beta = 0.99 ! discount rate
    double precision, parameter :: alpha = 0.36 ! capital share in C.D. production function
    double precision, parameter :: delta = 0.03 ! capital depreciation
    double precision, parameter :: sigma = 1.0 ! IES / risk aversion parameter
    double precision, parameter ::  pi = 3.1415926535  ! greek pi
    double precision :: Zstate(Lnum) ! aggregate TFP parameter
    double precision :: xi(1), yi(1), K_SS(1), L, R, W,  K_agg, K1_agg(1), R1, W1

    double precision, parameter :: rhoz = 0.7 ! serial corr of idiosync shocks
    double precision, parameter :: nstdevz = 1.0  ! stuff for tauchen
    double precision, parameter :: sigmaz = 0.075 ! std of idiosync shocks

    double precision, parameter :: rhoagg = 0.9 ! serial corr of agg shocks
    double precision, parameter :: nstdevagg = 1.0  ! stuff for tauchen
    double precision, parameter :: sigmaagg = 0.075 ! std of agg shocks


    double precision :: pr_mat_z(znum,znum), z0(znum), L_z(znum), pr_mat_L(Lnum,Lnum)



    contains




    !k_min = 0.d0
    !k_max = 10.d0


    !call linspace(K_GRID,k_min,k_max,k_size)
    !K_FUN = K_GRID**2

    !xi = 2.4
    !call pwl_value_1d ( k_size, K_GRID, K_FUN, ni, xi, yi )
    !print*, yi(1), xi**2


end module

subroutine init(K_SS,beta,delta,alpha)
        double precision, intent(in) :: beta, delta, alpha
        double precision :: K_SS
    K_SS = ((1/beta - (1-delta))/alpha)**(-1/(1-alpha))
end

! =============================================================================
!                           CHEBYSHEV INTERPOLATION
! =============================================================================

! ************************** BASIS FUNCTIONS **************************
subroutine TMATRIX(Tmat,zz,mm,n)
    implicit none
    integer, intent(in) :: mm,n
    double precision , intent(in) :: zz(n,1)
    double precision :: Tmat(mm,n)
    integer :: i,j

    do i=1,n
        Tmat(1,i) = 1.
        Tmat(2,i) = zz(i,1)
        do j = 3,mm
            Tmat(j,i) =  2.*zz(i,1)*Tmat(j-1,i) - Tmat(j-2,i)
        end do
    end do
end
! =============================================================================
!                              OLS REGRESSION
! =============================================================================
subroutine OLS(beta,X,y,n,dim)
    implicit none
    integer, intent(in) :: n, dim
    double precision, intent(in) :: X(n,dim), y(n,1)
    double precision :: beta(dim,1)
    double precision :: x_prod(dim,dim), x_inv(dim,dim), predict(dim,1)

    x_prod = matmul(transpose(X),X)
    call inverse(x_prod,x_inv,dim)
    predict =  matmul(transpose(X),y)
    beta = matmul(x_inv,predict)

end subroutine OLS
! =============================================================================
!                                   LINSPACE
! =============================================================================

subroutine linspace(z,x,y,n)
    implicit none
    !n = the dimension of the output vector
    !z = the n x 1 output vector, with equally spaced points between x and y
    !x = the minimum of the linear grid
    !y = the maximum of the linear grid

    !input/output declarations
    integer :: n
    double precision :: z(n),x,y
    !local declarations
    integer :: i
    double precision :: d

    d = (y-x)/dble(n-1)
    z(1) = x
    do i = 2,n-1
        z(i) = z(i-1) + d
    end do
    z(n) = y

    return
 end subroutine linspace

! =============================================================================
!                           LINEAR INTERPOLATION
! =============================================================================
subroutine pwl_value_1d ( nd, xd, yd, ni, xi, yi )
!    piecewise linear function which interpolates data (XD(I),YD(I))
!    Input, real ( kind = 8 ) XD(ND), the data points.
!    Input, real ( kind = 8 ) YD(ND), the data values.
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
    implicit none

    integer ( kind = 4 ) nd
    integer ( kind = 4 ) ni
    integer ( kind = 4 ) i
    integer ( kind = 4 ) k
    real ( kind = 8 ) t
    real ( kind = 8 ) xd(nd)
    real ( kind = 8 ) yd(nd)
    real ( kind = 8 ) xi(ni)
    real ( kind = 8 ) yi(ni)

    yi(1:ni) = 0.0D+00

    if ( nd == 1 ) then
        yi(1:ni) = yd(1)
        return
    end if

    do i = 1, ni
        if ( xi(i) <= xd(1) ) then
            t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
            yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)
        else if ( xd(nd) <= xi(i) ) then
            t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
            yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)
        else
            do k = 2, nd
                if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then
                  t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
                  yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
                exit
                end if
            end do
        end if
    end do

    return
end
! =============================================================================
!                           UTILITY FUNCTION
! =============================================================================
double precision function u(c,gamma)
implicit none

!input/output declarations
double precision :: c, gamma

!other declarations

    if (gamma .eq. 1.0) then
        u = log(c)
        else
        u = (c**(1-gamma))/(1-gamma)
    end if

end function u
! =============================================================================
!                           d U / d C
! =============================================================================
double precision function u_prime(c,gamma)
implicit none

!input/output declarations
double precision :: c, gamma


!other declarations
    if (gamma .eq. 1.0) then
        u_prime = 1/c
        else
        u_prime = c**(-gamma)
    end if

end function u_prime
! =============================================================================
!                           Interest Rate
! =============================================================================
subroutine int_rate(R, alpha,Z,K,L)
implicit none

!input/output declarations
double precision :: alpha, Z, K, L, R

!other declarations
    R = alpha*Z*(K/L)**(alpha-1)

end subroutine
! =============================================================================
!                               Wage
! =============================================================================
subroutine wage(W,alpha,A,K,L)
implicit none

!input/output declarations
double precision :: alpha, A, K, L,W

!other declarations
    W = (1-alpha)*A*(K/L)**(alpha)
end subroutine
! =============================================================================
!                      Discretization of Shocks
! =============================================================================
double precision function erfcc(x)
implicit none

!input/output declarations
double precision :: x

!other declarations
double precision :: t,z

z = abs(x)
t = 1.0 / (1.0 + 0.5*z)

erfcc = t * exp(-z * z-1.26551223+t*(1.00002368+t*(0.37409196+&
    t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+&
    t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))

if (x.lt.0.0) erfcc = 2.0 - erfcc
end function erfcc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function normcdf(x,mu,sigma)
implicit none

!input/output declarations
!x: the input value for the normal cdf
!mu: the mean of the normal distribution
!sigma: the standard deviation of the normal distribution
double precision :: x,mu,sigma,erfcc

!other declarations
double precision :: z

!standardized value ~ N(0,1)
z = (x-mu)/sigma

normcdf =  0.5 * erfcc( ( -1.0 * z ) / sqrt(2.0) )

end function normcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)
implicit none

!input/output declarations
integer, intent(in) :: znum
double precision, intent(in) :: rhoz,sigmaz,nstdevz
double precision, intent(out) :: pr_mat_z(znum,znum),z0(znum)

!other declarations
integer :: zct,zprimect
double precision :: zmin,zmax,gridinc,stdev,normcdf

!determine end points of the grid (log space)
stdev =  ((sigmaz**2.0)/(1-rhoz**2.0))**0.5
zmin = - nstdevz * stdev
zmax = nstdevz * stdev

!insert points into the grid (log space)
call linspace(z0,zmin,zmax,znum)
gridinc = z0(2)-z0(1)

!loop over z states
do zct=1,znum

    !insert transition matrix middle rows
    do zprimect=2,(znum-1)
        pr_mat_z(zct,zprimect) = &
            normcdf(z0(zprimect)+gridinc/2.0,rhoz*z0(zct),sigmaz) - &
            normcdf(z0(zprimect)-gridinc/2.0,rhoz*z0(zct),sigmaz)
    end do !zct

    !first interval and last interval take the rest of the weight
    pr_mat_z(zct,1) = normcdf(z0(1)+gridinc/2.0,rhoz*z0(zct),sigmaz)
    pr_mat_z(zct,znum) = 1.0 - normcdf(z0(znum)-gridinc/2.0,rhoz*z0(zct),sigmaz)

end do !zct

!round the transition matrix
do zct=1,znum
    pr_mat_z(zct,:) = pr_mat_z(zct,:)/sum(pr_mat_z(zct,:))
end do !zct

!convert grid back to z-space
!z0 = exp(z0)
end subroutine tauchen

! =============================================================================
!                      Computation of Inverse
! =============================================================================

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
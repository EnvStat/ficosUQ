subroutine rk4vec ( t0, m, u0, dt, f, u, dm)
!*****************************************************************************80
!
!! RK4VEC takes one Runge-Kutta step for a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T0, the current time.
!
!    Input, integer ( kind = 4 ) M, the dimension of the system.
!
!    Input, real ( kind = 8 ) U0(M), the solution estimate at the current time.
!
!    Input, real ( kind = 8 ) DT, the time step.
!
!    Input, external F, a subroutine of the form 
!      subroutine f ( t, m, u, uprime ) 
!    which evaluates the derivative UPRIME(1:M) given the time T and
!    solution vector U(1:M).
!
!    Output, real ( kind = 8 ) U(M), the fourth-order Runge-Kutta solution 
!    estimate at time T0+DT.
!
  use prec

  implicit none

!---Arguments
  real(kind=dp), intent(in)                 :: t0 ! Current time t0
  integer(kind=4), intent(in)               :: m  ! Number of equations
  real(kind=dp), dimension(m), intent(in)   :: u0 ! Solution estimate at current time
  real(kind=dp), intent(in)                 :: dt ! Time step dt
  external                                  :: f  ! Equations
  real(kind=dp), dimension(m), intent(out)  :: u  ! Solution at time t0+dt
  logical, intent(in)                       :: dm ! debug mode passthrough flag

!---Local variables
  real(kind=dp), dimension(m) :: f0, f1, f2, f3   ! Derivatives
  real(kind=dp), dimension(m) :: u1, u2, u3       ! Intermediate solutions
  real(kind=dp)               :: t1, t2, t3       ! Times


! Initialize u-vectors 
  u1 = u0
  u2 = u0
  u3 = u0

! Get four sample values of the derivative.
  call f ( t0, m, u0, f0, dm)

  t1 = t0 + dt / 2.0D+00
  u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
  call f ( t1, m, u1, f1, dm)

  t2 = t0 + dt / 2.0D+00
  u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
  call f ( t2, m, u2, f2, dm)

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  call f ( t3, m, u3, f3, dm)

! Combine them to estimate the solution U at time T1.
  u(1:m) = u0(1:m) + dt * ( f0(1:m) + 2.0D+00 * f1(1:m) + 2.0D+00 * f2(1:m) &
    + f3(1:m) ) / 6.0D+00

  return
end

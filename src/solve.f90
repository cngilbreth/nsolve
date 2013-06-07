program test
  implicit none
  integer,  parameter :: rk = kind(1d0)
  real(rk), parameter :: pi = 3.141592653589793_rk

  real(rk), parameter :: h = 1._rk/256
  real(rk), parameter :: xub = 8._rk
  real(rk), parameter :: xlb = 0._rk
  integer,  parameter :: l = 1
  integer,  parameter :: n = nint((xub-xlb)/h) + 1
  real(rk), parameter :: t = 1.e-10_rk

  real(rk) :: v(0:n), u(0:n), x, elb, eub, e, g
  integer :: i

  ! Compute v(x)
  ! Must include effective potential from orbital angular momentum here
  do i=0,n
     x = xlb + i*h
     v(i) = 0.5_rk * x**2 + 0.5_rk * l*(l+1)/(x+t)**2
  end do
  g = 0.5_rk

  ! Initial values (approximate)
  u(0) = 0._rk
  u(1) = h**(l+1)
  u(n) = 0._rk
  u(n-1) = exp(-(xub)**2/2)

  ! Eigenvalue window
  elb = 2.25_rk
  eub = 2.625_rk

  call solve1(h,g,v,elb,eub,u,e)
  open(unit=10,file='u.dat')
  do i=1,n
     write (10,*) i, u(i)
  end do
  close(10)

contains


  subroutine solve1(h,G,v,elb,eub,u,e)
    ! Find an eigenvalue e and eigenfunction u of the radial Schrodinger
    ! equation with effective potential v:
    !   -G u''(x) + (v(x)-e)u = 0,  x ≥ 0, u(0) = u0, u(inf) = uinf
    ! Input:
    !   h:   Grid step size
    !   G:   Multiplier of d²u/dx² term
    !   v:   Effective potential, evaluated at grid points (must include
    !        ℏ² l(l+1)/(2 m r^2) term)
    ! Input/output:
    !   elb, eub:
    !        On input: Initial upper and lower bounds for the eigenvalue. The
    !        [elb,eub] Must bracket exactly 1 eigenvalue, and u must have the
    !        same number of nodes on this interval.
    !        On output: A small interval [elb,eub] containing the eigenvalue
    !   u:   On input, u(0), u(1), u(n-1), and u(n) are set to appropriate
    !        values for starting the integrations.
    !        On output, the wavefunction evaluated at the grid points.
    ! Output:
    !   e:   Eigenvalue estimate, e = (elb + eub)/2
    implicit none
    real(rk), intent(in) :: h, G, v(:)
    real(rk), intent(inout) :: elb, eub, u(:)
    real(rk), intent(out) :: e

    real(rk), parameter :: accuracy = 1d-12

    real(rk) :: fmid, emid, flb, fub

    call fnumerov1(h,g,v,elb,u,flb)
    open(unit=10,file='ulb.dat')
    do i=1,n
       write (10,*) i, u(i)
    end do
    close(10)
    call fnumerov1(h,g,v,eub,u,fub)
    open(unit=11,file='uub.dat')
    do i=1,n
       write (11,*) i, u(i)
    end do
    close(11)
    if (.not. flb*fub < 0._rk) stop "solve1: Zero not properly bracketed"

    do while (abs(eub - elb) > accuracy)
       ! Try the midpoint in [xlb,xub]
       emid = (elb + eub) * 0.5_rk
       call fnumerov1(h,g,v,emid,u,fmid)
       write (0,*) elb, eub, fmid
       ! Decrease interval
       if (fmid*flb >= 0._rk) then
          elb = emid
       else
          eub = emid
       end if
    end do
    e = (elb + eub) * 0.5_rk
    ! Postconditions:
    !   1. emid = (xlb + xub)/2
    !   2. abs(eub - elb) < accuracy
    !   3. sign(f(elb)) .ne. sign(f(eub))
  end subroutine solve1



  subroutine fnumerov1(h,G,v,e,u,f)
    ! Numerically intergrate the homogeneous differential eigenvalue equation,
    !   -G u''(x) + (v(x)-e)u = 0, u(0) = u0, u(inf) = uinf
    ! from two sides, and compute the fitness function f associated with
    ! matching the two solutions at a single point.
    ! Input:
    !   h:       Grid step size (x(i+1)-x(i))
    !   G:       Arbitrary constant
    !   v:       Points v(x(i)), i=1,..,n
    !   e:       Trial eigenvalue
    ! Input/output:
    !   u:       On input, u(0), u(1), u(n-1), and u(n) are set to appropriate
    !            values for starting the integrations.
    !            On output, the two numerically integrated functions, rescaled
    !            so they match at the separation point.
    ! Notes:
    !   1. u is a solution to the equation with eigenvalue e when the fitness
    !      function f is zero (i.e., changes sign).
    !   2. This routine uses a 7-point formula for the numerical derivative.
    implicit none
    real(rk), intent(in) :: h, G, v(:), e
    real(rk), intent(inout) :: u(:)
    real(rk), intent(out) :: f

    real(rk) :: q(size(v)), S(size(v)), up1, up2, usave(4)
    integer :: i,n,isep

    n = size(v)
    isep = 0
    do i=1,n
       q(i) = (e - v(i))/G
       s(i) = 0._rk
       ! Select leftmost turning point
       if (q(i)*q(1) < 0 .and. isep == 0) isep = i
    end do
    if (isep == 0) stop "No turning point found"

    ! Left solution
    call numerov(h,q(1:isep+3),S(1:isep+3),+1,u(1:isep+3))
    u(1:isep+3) = u(1:isep+3)/u(isep)
    ! Numerical derivative
    up1 = (-u(isep-3) + 9*u(isep-2) - 45*u(isep-1) + 45*u(isep+1) - 9*u(isep+2) + u(isep+3))/(60*h)
    usave(1:4) = u(isep-3:isep)

    ! Right solution
    call numerov(h,q(isep-3:n),S(isep-3:n),-1,u(isep-3:n))
    u(isep-3:n) = u(isep-3:n)/u(isep)
    ! Numerical derivative
    up2 = (-u(isep-3) + 9*u(isep-2) - 45*u(isep-1) + 45*u(isep+1) - 9*u(isep+2) + u(isep+3))/(60*h)

    ! Restore overwritten values of u from first solution
    u(isep-3:isep) = usave(1:4)

    ! Compute matching condition
    ! f = u'_<(x0)/u_<(x0) - u'_>(x0)/u_>(x0)
    f = up1 - up2
  end subroutine fnumerov1



  subroutine numerov(h,q,S,dir,u)
    ! Integrate a second-order ODE initial value problem using the Numerov
    ! method.
    !   u''(x) + q(x)u = S(x), u(lb) = u0, u(lb+h) = u1.
    ! or
    !   u''(x) + q(x)u = S(x), u(ub) = u0, u(ub-h) = u1.
    ! Input:
    !   h:    dx
    !   q:    function q(x) evaluated at grid points
    !   S:    function S(x) evaluated at grid points
    !   dir:  Direction to integrate: backward (-1) or forward (+1)
    ! Output:
    !   u:    Solution of the differential equation
    implicit none
    real(rk), intent(in) :: h, q(:), S(:)
    integer,  intent(in) :: dir
    real(rk), intent(inout) :: u(:)

    real(rk) :: c
    integer :: n

    n = size(q)
    if (size(q) .ne. size(u) .or. size(q) .ne. size(S)) &
         stop "Error in numerov()"
    c = h**2/12._rk

    ! Eq. (3.8) in Koonin and Meredith
    if (dir > 0) then
       ! GIVEN: u(1), u(2)
       do i=2,n-1
          u(i+1) = (2._rk-10*c*q(i))*u(i) - (1+c*q(i-1))*u(i-1) &
               + c*(S(i+1) + 10*S(i) + S(i-1))
          u(i+1) = u(i+1) / (1+c*q(i+1))
       end do
    else
       ! GIVEN: u(n), u(n-1)
       do i=n-1,2,-1
          u(i-1) = (2._rk-10*c*q(i))*u(i) - (1+c*q(i+1))*u(i+1) &
               + c*(S(i+1) + 10*S(i) + S(i-1))
          u(i-1) = u(i-1) / (1+c*q(i-1))
       end do
    end if
  end subroutine numerov

end program test

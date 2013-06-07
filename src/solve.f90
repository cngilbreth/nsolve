program test
  implicit none
  integer, parameter :: rk = kind(1d0)

  real(rk) :: pi = 3.141592653589793_rk

  integer, parameter  :: n = 2048
  real(rk) :: v(0:n), u(0:n), xlb, xub, h, x, elb, eub, e, f, g
  integer :: i

  ! Bounds
  xlb = 0._rk; xub = 8._rk
  h = (xub-xlb)/n
  ! Compute v(x)
  do i=0,n
     x = xlb + i*h
     v(i) = 0.5d0 * x**2
  end do
  g = 0.5d0
  ! Initial values
  u(0) = 0.d0
  u(1) = h
  u(n) = 0.d0
  u(n-1) = exp(-(xub)**2/2)
  elb = 3.25d0
  eub = 3.625d0

  call solve(h,g,v,elb,eub,u,e)
  open(unit=10,file='u.dat')
  do i=1,n
     write (10,*) i, u(i)
  end do
  close(10)

contains


  subroutine solve(h,G,v,elb,eub,u,e)
    implicit none
    real(rk), intent(in) :: h, G, v(:)
    real(rk), intent(inout) :: elb, eub, u(:)
    real(rk), intent(out) :: e

    real(rk), parameter :: accuracy = 1d-12

    real(rk) :: fmid, emid, flb, fub

    call fnumerov0(h,g,v,elb,-1,u,flb)
    open(unit=10,file='ulb.dat')
    do i=1,n
       write (10,*) i, u(i)
    end do
    close(10)
    call fnumerov0(h,g,v,eub,-1,u,fub)
    open(unit=11,file='uub.dat')
    do i=1,n
       write (11,*) i, u(i)
    end do
    close(11)
    if (.not. flb*fub < 0._rk) stop "solve: Zero not properly bracketed"

    do while (abs(eub - elb) > accuracy)
       ! Try the midpoint in [xlb,xub]
       emid = (elb + eub) * 0.5_rk
       call fnumerov0(h,g,v,emid,-1,u,fmid)
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
  end subroutine solve



  subroutine fnumerov0(h,G,v,e,dir,u,f)
    ! Numerically intergrate the homogeneous differential eigenvalue equation,
    !   -G u''(x) + (v(x)-e)u = 0, u(0) = u0, u(inf) = uinf
    ! from one side only, and compute a fitness function f associated with the
    ! boundary condition.
    ! Input:
    !   h:       Grid step size (x(i+1)-x(i))
    !   G:       Arbitrary constant
    !   v:       Points v(x(i)), i=1,..,n
    !   e:       Trial eigenvalue
    ! Input/output:
    !   u:       On input, u(1), u(n) are set to appropriate boundary values.
    !            If dir = +1, u(2) is set to an initial value to begin the
    !            integration.
    !            If dir = -1, u(n-1) is set to an initial value to begin
    !            the integration.
    !            On output, u(1:n-2) or (3:n) is the numerically integrated
    !            function. The input points are left unchanged.
    ! Notes:
    !   u is a solution to the equation with eigenvalue e when the fitness
    !   function f is zero (i.e., changes sign).
    implicit none
    real(rk), intent(in) :: h, G, v(:), e
    integer,  intent(in) :: dir
    real(rk), intent(inout) :: u(:)
    real(rk), intent(out) :: f

    real(rk) :: q(size(v)), S(size(v)), u0, un


    q = (e - v)/G
    S = 0
    if (dir < 0) then
       u0 = u(1)
       call numerov(h,q,S,dir,u)
       !u(1:n-2) = u(1:n-2) / sqrt(dot_product(u(1:n-2),u(1:n-2)))
       f = u(1)
       u(1) = u0
    else
       un = u(n)
       call numerov(h,q,S,dir,u)
       !u(3:n-1) = u(3:n-1) / sqrt(dot_product(u(3:n-1),u(3:n-1)))
       f = u(n)
       u(n) = un
    end if
  end subroutine fnumerov0



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
    !   u is a solution to the equation with eigenvalue e when the fitness
    !   function f is zero (i.e., changes sign).
    implicit none
    real(rk), intent(in) :: h, G, v(:), e
    real(rk), intent(inout) :: u(:)
    real(rk), intent(out) :: f

    real(rk) :: q(size(v)), S(size(v)), uisep, up1, up
    integer :: i,n,isep

    n = size(v)
    isep = 0
    do i=1,n
       q(i) = (e - v(i))/G
       s(i) = 0.d0
       ! Select leftmost turning point
       if (q(i)*q(1) < 0 .and. isep == 0) isep = i
    end do
    if (isep == 0) stop "No turning point found"
    call numerov(h,q(1:isep),S(1:isep),+1,u(1:isep))

    u(1:isep) = u(1:isep)/u(isep)
    call numerov(h,q(isep:n),S(isep:n),-1,u(isep:n))
    u(isep:n) = u(isep:n)/u(isep)

    ! Compute matching condition
    ! f = u'_<(x0)/u_<(x0) - u'_>(x0)/u_>(x0)
    up  = (u(isep) - u(isep-1))/h
    up1 = (u(isep+1) - u(isep))/h
    f = up - up1
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

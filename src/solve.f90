! Module for solving the 1-dimension Schrodinger equation, e.g., the radial
! equation for a central potential.
module solve
  implicit none

  integer,  parameter, private :: rk = kind(1d0)
  real(rk), parameter, private :: pi = 3.141592653589793_rk

  ! Accuracy of window for bracketing nodes of the wavefunction
  real(rk), parameter, private :: node_accuracy = 1.e-4_rk
  ! Accuracy with which to determine the eigenvalues
  real(rk), parameter, private :: eval_accuracy = 1.e-12_rk
  ! Maximum window [node_emin,node_emax] for bracketing wavefunction with
  ! a fixed number of nodes
  real(rk), parameter, private :: node_emin = -100._rk
  real(rk), parameter, private :: node_emax = 100._rk

  public :: bracket_nodes2, solve2

contains


  subroutine bracket_nodes2(h,v,k,u,elb,eub,ierr)
    ! Find the largest energy interval [elb,eub] between which the wavefunction
    ! has a given number of nodes.
    ! Input:
    !   h:  Grid spacing
    !   v:  Effective potential, evaluated at grid points (must include
    !       ℏ² l(l+1)/(2 m r^2) term)
    !   k:  Number of nodes
    !   u:  On input, u(0), u(1), u(n-1), and u(n) are set to appropriate
    !       values for starting the integrations.
    !       On output, the rest of u is destroyed.
    ! Output:
    !   elb,eub: If ierr == 0, energy window within which the wavefunction
    !            has exactly k nodes.
    !   ierr:  0 on success, nonzero on failure.
    ! Notes:
    !   1. Only works for bound states currently.
    !   2. Can have difficulty with highly singular potentials.
    !   3. It assumes the number of nodes is a monotonic function of the energy.
    implicit none
    real(rk), intent(in) :: h, v(:)
    integer,  intent(in) :: k
    real(rk), intent(inout) :: u(:)
    real(rk), intent(out) :: elb,eub
    integer,  intent(out) :: ierr

    real(rk) :: elb_ub, eub_lb, emid, emin, emax
    integer  :: count, i

    ! Absolute min/max for energy window
    ! These values are determined as the classically allowed region,
    ! plus some wiggle room.
    emin = max(node_emin, minval(v) + abs(minval(v))*0.1)
    if (abs(emin) < 1E-12_rk) emin = 0.0001_rk
    emax = min(node_emax, maxval(v) - abs(maxval(v))*0.1)

    call count_nodes2(h,v,emin,u,count)
    if (count .ge. max(k,1)) then
       write (0,'(1x,a,i0,a)') &
            "Error in bracket_nodes: Can't bracket all ", k, "-node &
            &wavefunctions within [emin,emax]"
       write (0,'(1x,a)') "  Try adjusting node_emin"
       write (0,'(1x,a,es17.10,a,es17.10)') &
            "Currently: emin=",emin, ", emax=", emax
       write (0,'(1x,a)') "  See file 'umin.dat' for wavefunction at emin"
       open(unit=10,file='umin.dat')
       do i=1,size(u)
          write (10,*) i, u(i)
       end do
       close(10)
       ierr = -1
       return
    end if

    call count_nodes2(h,v,emax,u,count)
    if (count <= k) then
       write (0,'(1x,a,i0,a)') &
            "Can't bracket all ", k, "-node wavefunctions within [emin,emax]"
       write (0,'(1x,a)') "Try adjusting node_emax"
       write (0,'(1x,a,es17.10,a,es17.10)') &
            "Currently: emin=",emin, ", emax=", emax
       write (0,'(1x,a)') "See file 'umax.dat' for wavefunction at emax"
       open(unit=10,file='umax.dat')
       do i=1,size(u)
          write (10,*) i, u(i)
       end do
       close(10)
       ierr = -1
       return
    end if
    ierr = 0

    ! Find accurate lower bound
    ! For k=0 (no nodes), we keep 'emin' as the lower bound.
    elb = emin
    if (k > 0) then
       elb_ub = emax
       do while (abs(elb - elb_ub) > node_accuracy)
          ! count(elb) < k
          ! count(elb_ub) >= k
          emid = (elb + elb_ub)/2
          call count_nodes2(h,v,emid,u,count)
          if (count < k) then
             elb = emid
          else
             elb_ub = emid
          end if
       end do
       elb = elb_ub
    end if

    ! Find accurate upper bound
    eub_lb = emin
    eub = emax
    do while (abs(eub - eub_lb) > node_accuracy)
       ! count(ub_lb) <= k
       ! count(eub) > k
       emid = (eub + eub_lb)/2
       call count_nodes2(h,v,emid,u,count)
       if (count <= k) then
          eub_lb = emid
       else
          eub = emid
       end if
    end do
    eub = eub_lb
  end subroutine bracket_nodes2


  subroutine count_nodes2(h,v,e,u,count)
    ! Count the number of nodes in the two-sided integrated solution to the
    ! differential equation
    !   -(1/2) u''(x) + (v(x)-e)u = 0,  x ≥ 0, u(0) = u0, u(inf) = uinf
    implicit none
    real(rk), intent(in) :: h
    real(rk), intent(in) :: v(:), e
    real(rk), intent(inout) :: u(:)
    integer,  intent(out) :: count

    integer  :: i, n
    real(rk) :: f

    call fnumerov2(h,v,e,u,f)
    n = size(v)
    count = 0
    ! We assume "nodes" are not right at the edges. If they are, need to
    ! decrease grid spacing.
    do i=2,n-2
       if (u(i+1)*u(i) < 0 .or. u(i) == 0._rk) then
          count = count + 1
       end if
    end do
  end subroutine count_nodes2



  subroutine solve2(h,v,elb,eub,u,e,ierr)
    ! Find an eigenvalue e and eigenfunction u of the radial Schrodinger
    ! equation,
    !    -(1/2) u''(x) + (v(x)-e)u = 0,  x ≥ 0, u(0) = u0, u(inf) = uinf
    !  with effective potential v, via the Numerov method.
    ! Input:
    !   h:   Grid step size
    !   v:   Effective potential, evaluated at grid points (must include
    !        ℏ² l(l+1)/(2 m r^2) term)
    ! Input/output:
    !   elb, eub:
    !        On input: Initial upper and lower bounds for the eigenvalue. The
    !        interval [elb,eub] Must bracket exactly 1 eigenvalue, and u must
    !        have the same number of nodes on this interval. See bracket_nodes.
    !        On output: A small interval [elb,eub] containing the eigenvalue
    !   u:   On input, u(0), u(1), u(n-1), and u(n) are set to appropriate
    !        boundary values for starting the integrations.
    !        On output, the wavefunction evaluated at the grid points.
    ! Output:
    !   e:   Eigenvalue estimate
    implicit none
    real(rk), intent(in) :: h, v(:)
    real(rk), intent(inout) :: elb, eub, u(:)
    real(rk), intent(out) :: e
    integer,  intent(out) :: ierr

    real(rk) :: fmid, emid, flb, fub
    integer :: i, n

    n = size(v)

    call fnumerov2(h,v,elb,u,flb)

    call fnumerov2(h,v,eub,u,fub)

    if (.not. flb*fub < 0._rk) then
       write (0,*) "Error in solve1: eigenvalue not bracketed in [elb,eub]"
       write (0,*) "  See 'ulb.dat', 'uub.dat' and 'v.dat' for more info"
       open(unit=10,file='v.dat')
       do i=1,n
          write (10,*) i, v(i)
       end do
       close(10)
       open(unit=10,file='ulb.dat')
       write (10,'(a,es20.10)') "# elb: ", elb
       do i=1,n
          write (10,*) i, u(i)
       end do
       close(10)
       ierr = 1
       return
       open(unit=10,file='uub.dat')
       write (10,'(a,es20.10)') "# eub: ", eub
       do i=1,n
          write (10,*) i, u(i)
       end do
       close(10)
    end if

    write (*,'(a,a25,tr2,a25,tr4,a13)') '#', "elb", "eub", "f(midpoint)"
    do while (abs(eub - elb) > eval_accuracy)
       ! Try the midpoint in [xlb,xub]
       emid = (elb + eub) * 0.5_rk
       call fnumerov2(h,v,emid,u,fmid)
       write (*,'(a,es25.16,tr2,es25.16,tr4,es13.5)') '#', elb, eub, fmid
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
    open(unit=10,file='u.dat')
    write (10,'(a)') '#  index            u(index)'
    do i=1,n
       write (10,'(i8,es20.10)') i, u(i)
    end do
    close(10)
    write (*,'(a)') "# wavefunction written to u.dat"
  end subroutine solve2



  subroutine fnumerov2(h,v,e,u,f)
    ! Numerically intergrate the homogeneous differential eigenvalue equation,
    !   - u''(x)/2 + (v(x)-e)u = 0, u(0) = u0, u(inf) = uinf
    ! from two sides, and compute the fitness function f associated with
    ! matching the two solutions at a single point.
    ! Input:
    !   h:       Grid step size (x(i+1)-x(i))
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
    !   2. This routine uses a 7-point formula for the numerical derivative,
    !      which is accurate to the same order as the Numerov integration.
    implicit none
    real(rk), intent(in) :: h, v(:), e
    real(rk), intent(inout) :: u(:)
    real(rk), intent(out) :: f

    real(rk) :: q(size(v)), S(size(v)), up1, up2, usave(4)
    integer :: i,n,isep

    n = size(v)
    isep = 0
    do i=1,n
       q(i) = (e - v(i))*2
       s(i) = 0._rk
       ! Select leftmost turning point
       if (q(i)*q(1) < 0 .and. isep == 0) isep = i
    end do
    if (isep == 0 .or. isep <= 3 .or. isep >= n-3) then
       write (0,*) "fnumerov2(): No turning point found."
       stop
    end if

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
  end subroutine fnumerov2



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
    integer :: n, i

    n = size(q)
    if (size(q) .ne. size(u) .or. size(q) .ne. size(S)) &
         stop "Error in numerov()"
    c = h**2/12._rk

    ! Eq. (3.8) in Koonin and Meredith
    ! TODO: Use compensated summation for this
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

end module solve

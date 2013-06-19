! ho.f90: Find eigenstates of 3D harmonic oscillator
program ho
  use solve
  implicit none

  integer,  parameter :: rk = kind(1d0)
  ! Grid stepsize (in units of oscillator length d = \sqrt{ℏ/mω})
  real(rk), parameter :: h  = 1._rk/512
  ! Number of nodes desired in the wavefunction
  integer,  parameter :: nodes = 2
  ! Angular momentum
  integer,  parameter :: l = 0
  ! A small parameter for computing l(l+1)/(2 m r^2 + t)
  ! (avoids divergence)
  real(rk), parameter :: t = h**2/16
  ! Interval [xlb,xub] over which to solve the Schrodinger eqn
  ! (in units of oscillator length)
  real(rk), parameter :: xlb = 0._rk
  real(rk), parameter :: xub = 10._rk

  real(rk) :: e
  integer :: ierr

  call solve_ho(h,l,nodes,e,ierr)
  if (ierr .eq. 0) then
     write (*,*) "energy: ", e
  else
     write (*,*) "Error"
  end if

contains


  subroutine solve_ho(h,l,k,e,ierr)
    ! Find an energy eigenstate of a particle in a 3d harmonic oscillator
    ! potential by solving the radial Schrodinger equation.
    ! Input:
    !   h:  Grid spacing
    !   l:  Angular momentum
    !   k:  Number of nodes of the radial wavefunction
    ! Output:
    !   e:     Energy eigenvalue
    !   ierr:  If ierr == 0, routine was successful. If err .ne. 0,
    !          unsuccessful.
    implicit none
    real(rk), intent(in)  :: h
    integer,  intent(in)  :: l, k
    real(rk), intent(out) :: e
    integer,  intent(out) :: ierr

    real(rk) :: x, elb, eub
    integer  :: i, n
    real(rk), allocatable :: v(:), u(:)

    n = idnint((xub-xlb)/h) + 1
    allocate(v(0:n),u(0:n))

    ! Potential v(x)
    ! Must include effective potential from orbital angular momentum here
    do i=0,n
       x = xlb + i*h
       v(i) = 0.5_rk * x**2 + 0.5_rk * l*(l+1)/(x**2 + t)
    end do

    ! Initial values (approximate)
    u(0) = 0._rk
    u(1) = h**(l+1)
    u(n) = 0._rk
    u(n-1) = exp(-(xub)**2/2)

    ! Determine the energy window [elb,eub] within which the two-sided
    ! wavefunction (not eigenfunction) determined by the Numerov method has
    ! exactly k nodes
    call bracket_nodes2(h,v,nodes,u,elb,eub,ierr)
    if (ierr .ne. 0) then
       write (0,*) "Couldn't bracket # of nodes for wavefunction"
       return
    end if

    ! Solve the equation
    call solve2(h,v,elb,eub,u,e,ierr)
  end subroutine solve_ho

end program ho

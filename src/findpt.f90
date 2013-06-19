! findpt.f90: Program for finding the interaction strength v0 for the
! poschl-teller potential in the unitary limit.
! v0.6, June 19, 2013
!
! Copyright (c) 2013 Christopher N. Gilbreth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

program findpt
  use solve
  implicit none
  integer,  parameter :: rk = kind(1d0)
  ! Grid stepsize
  real(rk), parameter :: h  = 1._rk/65536
  ! Number of nodes in wavefunction desired
  integer,  parameter :: nodes = 0
  ! Exact energy to reproduce
  real(rk), parameter :: exact_energy = 0.5_rk
  ! Effective range of interaction
  real(rk), parameter :: r0 = 1.0d0
  ! Interval [v0lb,v0ub] in which to search for v0
  real(rk), parameter :: v0lb = -10._rk
  real(rk), parameter :: v0ub = -5.0_rk
  ! Interval [xlb,xub] over which to solve the equation
  real(rk), parameter :: xlb = 0._rk
  real(rk), parameter :: xub = 12._rk
  ! A small parameter for computing l(l+1)/(2 m r^2 + t)
  real(rk), parameter :: t = h**2/4


  ! With h = 1/65536, exact_energy = 1/2:
  !  r0                  V0
  ! 1.0        -6.565820993
  ! 0.9        -6.114178264
  ! 0.8        -5.701018261
  ! 0.7        -5.327309888
  ! 0.6        -4.994431426
  ! 0.5        -4.704281279
  ! 0.4        -4.459377729
  ! 0.3        -4.262890615
  ! 0.2        -4.118495286
  ! 0.1        -4.029901582
  ! 0.08       -4.019159458
  ! 0.06       -4.010787114
  ! 0.04       -4.004797446
  ! 0.02       -4.001199840
  ! 0.01       -4.000299990
  ! 0.005      -4.000074999

  ! Note the window [v0lb,v0ub] must sometimes be adjusted to a narrow window
  ! around v0 for this to work well.

  real(rk) :: v0
  integer :: ierr

  call findv0_pt(r0,h,nodes,exact_energy,v0,ierr)
  if (ierr .ne. 0) then
     write (0,*) "Error finding v0."
  else
     write (*,*) v0
  end if

contains

  subroutine findv0_pt(r0,h,k,e,v0,ierr)
    implicit none
    real(rk), intent(in) :: r0,h
    integer,  intent(in) :: k
    real(rk), intent(in) :: e
    real(rk), intent(out) :: v0
    integer,  intent(out) :: ierr

    real(rk), parameter :: accuracy = 1d-10

    real(rk), allocatable :: u(:)
    real(rk) :: emid, v0mid, elb, eub, v0ub1, v0lb1
    integer :: i

    call solve_pt(v0lb,r0,h,k,u,elb,ierr)
    if (ierr .ne. 0) then
       write (0,*) "Error solving equation for v0lb"
       return
    end if
    call solve_pt(v0ub,r0,h,k,u,eub,ierr)
    if (ierr .ne. 0) then
       write (0,*) "Error solving equation for v0ub"
       return
    end if

    if ((elb-e)*(eub-e) >= 0) then
       write (*,*) "findv0_pt: Desired energy not bracketed. &
            &Need to expand [v0lb,v0ub]"
       ierr = -1
       return
    end if

    v0ub1 = v0ub
    v0lb1 = v0lb
    do while (abs(v0ub1 - v0lb1) > accuracy)
       v0mid = (v0ub1 + v0lb1)/2
       call solve_pt(v0mid,r0,h,k,u,emid,ierr)
       if (ierr .ne. 0) return
       if ((emid-e)*(elb-e) >= 0) then
          v0lb1 = v0mid
       else
          v0ub1 = v0mid
       end if
    end do
    v0 = (v0lb1 + v0ub1)/2

    open(unit=10,file='u.dat')
    do i=1,size(u)
       write (10,'(i0,t10,es17.10)') i, u(i)
    end do
    close(10)
  end subroutine findv0_pt


  subroutine solve_pt(v0,r0,h,k,u,e,ierr)
    ! Find an energy eigenstate of a two particles in a harmonic trap with a
    ! poschl-teller interaction potential.
    ! Input:
    !   v0: Strength of interaction
    !   r0: Range of interaction
    !   h:  Grid spacing
    !   k:  Number of nodes of the wavefunction
    ! Output:
    !   e:     Lowest relative-motion energy eigenvalue
    !   ierr:  If ierr == 0, routine was successful. If err .ne. 0,
    !          unsuccessful.
    implicit none
    real(rk), intent(in)  :: v0, r0, h
    integer,  intent(in)  :: k
    real(rk), allocatable, intent(out) :: u(:)
    real(rk), intent(out) :: e
    integer,  intent(out) :: ierr

    real(rk) :: x, elb, eub
    integer  :: i, l, n
    real(rk), allocatable :: v(:)

    write (0,'(a,es20.10,a)') &
         "Solving PT equation with v0 = ", v0, " ... "

    n = idnint((xub-xlb)/h) + 1
    allocate(v(0:n),u(0:n))

    ! Potential v(x)
    ! Must include effective potential from orbital angular momentum here
    l = 0
    do i=0,n
       x = xlb + i*h
       v(i) = 0.5_rk * x**2 + 0.5_rk * l*(l+1)/(x+t)**2 &
            + 0.5_rk * v0 / (r0**2 * cosh(sqrt(2._rk)*x/r0)**2)
    end do

    ! Initial values (approximate)
    u(0) = 0._rk
    u(1) = h**(l+1)
    u(n) = 0._rk
    u(n-1) = exp(-(xub)**2/2)

    ! Determine energy window [elb,eub] within which the two-sided wavefunction
    ! determined by the Numerov method has exactly k nodes
    call bracket_nodes2(h,v,k,u,elb,eub,ierr)
    if (ierr .ne. 0) then
       write (0,*) "Couldn't bracket # of nodes for wavefunction"
       return
    end if
    !elb = 0.3
    !eub = 0.7

    call solve2(h,v,elb,eub,u,e,ierr)
  end subroutine solve_pt


end program findpt

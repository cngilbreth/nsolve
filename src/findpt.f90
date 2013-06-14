! Program for finding the interaction strength v0 for the poschl-teller
! potential in the unitary limit.
program findpt
  use solve
  implicit none
  integer,  parameter :: rk = kind(1d0)
  real(rk), parameter :: h  = 1._rk/65536
  integer,  parameter :: nodes = 0
  real(rk), parameter :: exact_energy = 0.5_rk
  real(rk), parameter :: r0 = 0.3d0
  real(rk), parameter :: v0lb = -4.4_rk
  real(rk), parameter :: v0ub = -4.02_rk

  ! With h = 1/65536:
  !  r0                  V0
  ! 0.3        -4.262890615
  ! 0.2        -4.118495286
  ! 0.1        -4.029901582
  ! 0.08       -4.019159458
  ! 0.06       -4.010787114
  ! 0.04       -4.004797446
  ! 0.02       -4.001199840
  ! 0.01       -4.000299990
  ! 0.005      -4.000074999

  real(rk) :: e, v0
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
    real(rk) :: emid, v0mid, elb, eub, v0ub1, v0lb1


    call solve_pt(v0lb,r0,h,k,elb,ierr)
    if (ierr .ne. 0) then
       write (0,*) "Error solving equation for v0lb"
       return
    end if
    call solve_pt(v0ub,r0,h,k,eub,ierr)
    if (ierr .ne. 0) then
       write (0,*) "Error solving equation for v0ub"
       return
    end if

    if ((elb-e)*(eub-e) >= 0) then
       write (*,*) "findv0_pt: Desired energy not bracketed. Need to expand [v0lb,v0ub]"
       ierr = -1
       return
    end if

    v0ub1 = v0ub
    v0lb1 = v0lb
    do while (abs(v0ub1 - v0lb1) > accuracy)
       v0mid = (v0ub1 + v0lb1)/2
       call solve_pt(v0mid,r0,h,k,emid,ierr)
       if (ierr .ne. 0) return
       if ((emid-e)*(elb-e) >= 0) then
          v0lb1 = v0mid
       else
          v0ub1 = v0mid
       end if
    end do
    v0 = (v0lb1 + v0ub1)/2
  end subroutine findv0_pt


end program findpt

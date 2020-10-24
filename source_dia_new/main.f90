!==============================================================================
program surfgen
  use makesurfdata
  implicit none

!  call random_seed()
  
  call readinput()
  call readdisps
  call makesurf

  stop
end program surfgen
!==============================================================================
SUBROUTINE readinput()
  use progdata
  use hddata
  use pip
  use makesurfdata
  implicit none
  integer :: fid
  integer :: fits(10)
  integer :: i, ifits
  namelist /fitting/ npoints, epmax, nc, useGPU, rhog, fits, restart, scalemod, &
      morse

  !progdata
  natoms=4
  printlvl=5

  !hddata
  nstates=3
  call readcoordinput()
  ncoord = 30 + noopc ! 30 PIP + Number of oop coords in file 
  print *, "Total number of coordinates: ", ncoord
  print *, "Number of user defined OOP coordinates:", noopc
  
  !readdisps
  enfptn  =  'energy.all'
  gmfptn  =  'geom.all'
  hdfptn  =  'hd.all'
  eshift = 76.6476575743

  ! pip modes
  scalemod=-1
  morse=1.0
  
  ! Set up energy scaling
  energyT(1:5) = (/        10.0D3, 20.0D3, 30.0D3,  100.0D3,  1000.D3 /) / ev2cm
  highEScale(1:5) = (/ 1.0D0, 0.5D0,  0.2D0,  0.1D0,   0.1D-1  /)

  ! restart from WB-restart.txt ?
  restart = .false.
  
  fits=0
  call FLUnit(fid)
  open(fid,file='fitting.para',delim='APOSTROPHE')
  read(unit=fid,nml=fitting)
  close(fid)

  call FLUnit(fid)
  open(fid,file='fitting.nml',delim='APOSTROPHE')
  write(unit=fid,nml=fitting)
  close(fid)

  if (.not. allocated(trainings)) allocate(trainings(nc))

  ! Run all trainings by default
  trainings = 1
  ! Run only trainings selected by user
  ifits = 0
  do i = 1, 10
     ifits = ifits + fits(i)
  end do
  if (ifits .ne. 0) then
     print *, "Selecting trainings from full set."
     trainings =  0
     do i = 1, 10
        if (fits(i) .ne. 0) then
           trainings(fits(i)) = 1
        end if
     end do
  end if

  !pip
  npip=ncoord-noopc
  pipbasfl='MOL_1_3_4.BAS'
  call initpip

  !ANN
  ANN%structfl='struct'
  call ANN%init()
  print *, "ANN%RX = ", ANN%RX
  print *, "npip   = ", npip
  if (ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if (ANN%RY .ne. lastlayer) then
      print "(A,i3)", "ANN%RY    = ", ANN%RY
      print "(A,i3)", "lastlayer = ", lastlayer
      stop 'ANN%RY .ne. lastlayer'
  end if
  

  return
end SUBROUTINE readinput
!================================================================================

!==============================================================================
program surfgen
  use makesurfdata
  use pip
  implicit none

  call random_seed()

  call readinput
  call readdisps
  call makesurf

  stop
end program surfgen
!==============================================================================
SUBROUTINE readinput
  use progdata
  use hddata
  use pip
  use makesurfdata
  implicit none
  integer :: fid

  namelist /fitting/ npoints, epmax, useGPU, w_energy, w_grad, w_fij,   &
    lambda, ptwthreshold, energyT, highEScale, deg_cap,&
    gradcutoff, cpcutoff, intGradT, intGradS, nrmediff, ediffcutoff,    &
    scalemod, morse

  !readEner: uniform shift on ab initio energies
  eshift = 76.6476575743 !Hartree
  deg_cap= 2.d-3 !au
  gorder=200 / EV2CM !eV
  gradcutoff= 100000.!eV
  cpcutoff = -1.
  nrmediff = 2000 / EV2CM
  ediffcutoff = 2000 / EV2CM
  
  !progdata
  natoms=4
  printlvl=5

  !hddata
  nstates=3
  call readcoordinput()
  ncoord=30 + noopc
  print *, "Total number of coordinates: ", ncoord
  print *, "Number of user defined OOP coordinates:", noopc

  !readdisps
  enfptn   =  'energy.all'
  gmfptn   =  'geom.all'
  grdfptn  =  'cartgrd.drt1.state$.all'
  cpfptn   =  'cartgrd_total.drt1.state$.drt1.state$.all'

  ! pip modes
  scalemod=-1
  morse=1.0
  
  !makeLocalIntCoord
  intGradT = 1.d-2
  intGradS = 1.d-2

  !number of internal degrees
  nvibs=3*natoms-6

  !OrthGH_ab
  deggrdbinding = .true.
  gcutoff = 1.d-7

  call FLUnit(fid)
  open(fid,file='fitting.para',delim='APOSTROPHE')
  read(unit=fid,nml=fitting)
  close(fid)

  call FLUnit(fid)
  open(fid,file='fitting.nml',delim='APOSTROPHE')
  write(unit=fid,nml=fitting)
  close(fid)

  energyT=energyT/ev2cm
  energyT_en=energyT
  highEScale_en=highEScale

  deg_cap = deg_cap * AU2EV
  gradcutoff = gradcutoff / EV2CM
  cpcutoff =  cpcutoff / EV2CM
  if (cpcutoff .le. 0) cpcutoff = gradcutoff

  !pip
  npip=ncoord - noopc
  pipbasfl='MOL_1_3_4.BAS'

  call initpip
  !ANN
  ANN%structfl='struct'
  call ANN%init()
  if(ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if(ANN%RY .ne. lastlayer) stop 'ANN%RY .ne. 19'


  return
end SUBROUTINE readinput
!================================================================================

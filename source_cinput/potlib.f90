!************************************************************
! potdata: module containing potential evaluation variables
module potdata
  implicit none
  
  ! nonlinear parameters
  integer :: ncons
  double precision, dimension(:), allocatable :: coef, coefc
  
contains
  
end module potdata

!*
! evaluate_nnsurf: evaluate the nn surface potential.
! Input:
!  cgeom = cartesian geometry
! Output:
!  e     = adiabatic energy
!  gd    = adiabatic energy gradients
!  hd    = diabatic hamiltonian
!  dhd   = gradient of diabatic hamiltonian
!*
subroutine evaluate_nnsurf(cgeom, e, gd, hd, dhd)
  use progdata, only: natoms
  use hddata, only: nstates, ncoord, EvaluateHD0, ANN
  implicit none

  ! .. Input/Output variables ..
  double precision, dimension(3*natoms), intent(in) :: cgeom
  double precision, dimension(nstates),  intent(out):: e
  double precision, dimension(nstates,nstates), intent(out) :: hd
  double precision, dimension(3*natoms,nstates,nstates), intent(out) :: gd
  double precision, dimension(3*natoms,nstates,nstates), intent(out) :: dhd

  ! .. Local variables ..
  double precision, dimension(nstates,nstates)        :: evec, hmat2, dhtmp
  double precision, dimension(nstates)                :: eerr
  double precision, dimension(ncoord,nstates,nstates) :: dhmat
  double precision, dimension(ncoord)                 :: igeom
  double precision, dimension(ncoord,3*natoms)        :: bmat
  integer                                             :: m, lwork, liwork, info
  integer, dimension(2*nstates)                       :: isuppz
  double precision, dimension(nstates*(nstates+26))   :: work
  integer, dimension(10*nstates)                      :: iwork
  double precision                                    :: bohr2ang
  
  lwork = nstates*(nstates + 26)
  liwork= nstates*10
  bohr2ang = 0.529177249

  call buildWBmat(natoms, ncoord, cgeom, igeom, bmat)

  call EvaluateHd0(igeom, hd, dhd)
  
  
  e   = 0d0
  gd  = 0d0
  hd  = 0d0
  dhd = 0d0

  
  
  return
end subroutine evaluate_nnsurf

!*
! evaluate_nnsurf77: FORTRAN77 wrapper for evaluate_nnsurf.
! Input:
!  mnatm = maximum number of atoms
!  mnst  = maximum number of states
!  cgeom = cartesian geometry
! Output:
!  e     = adiabatic energy
!  gd    = adiabatic energy gradients
!  hd    = diabatic hamiltonian
!  dhd   = gradient of diabatic hamiltonian
!*
subroutine evaluate_nnsurf77(mnatm,mnst,cgeom,e,gd,hd,dhd)
  use progdata, only: natoms
  use hddata, only: nstates
  implicit none
  integer, intent(in) :: mnatm, mnst
  double precision, dimension(3*mnatm), intent(in) :: cgeom
  double precision, dimension(mnst),      intent(out) :: e
  double precision, dimension(mnst,mnst), intent(out) :: hd
  double precision, dimension(3*mnatm,mnst,mnst),intent(out) :: gd, dhd
  ! local packed variables
  double precision, dimension(3*natoms)        :: cgpck
  double precision, dimension(nstates,nstates) :: hdpck
  double precision, dimension(3*natoms,nstates,nstates) :: gdpck, dhdpck
  
  if ((natoms .gt. mnatm) .or. (nstates .gt. mnst)) then
    print *, "evaluate_nnsurf77: natoms > max atoms || nstates > max states."
    stop
  end if
  cgpck = cgeom(1:3*natoms)
  call evaluate_nnsurf(cgpck, e, gdpck, hdpck, dhdpck)
  hd(1:nstates,1:nstates) = hdpck
  gd( 1:3*natoms,1:nstates,1:nstates) = gdpck
  dhd(1:3*natoms,1:nstates,1:nstates) = dhdpck
  return
end subroutine evaluate_nnsurf77

!*
! getinfo: get atom number and state information
!*
subroutine getinfo(nat, nst)
  use hddata,   only: nstates
  use progdata, only: natoms
  integer, intent(out) :: nat, nst
  nat = natoms
  nst = nstates
end subroutine getinfo

!*
! initpotential: initialize potential energy surface
!*
subroutine initpotential()
  call prepot()
end subroutine initpotential
  
!*
! prepot: initialize the nn potential energy surface.
!*
subroutine prepot()
  use hddata
  use potdata
  use progdata
  use pip
  use nn_class
  implicit none
  integer :: ios
  integer :: jobtype

  ! progdata
  natoms=4
  ! hddata
  nstates=3
  call readcoordinput()
  ncoord=30 +  noopc
  print *, "Total number of coordinates: ", ncoord
  print *, "Number of user defined OOP coordinates:", noopc

  ! PIP
  npip=ncoord - noopc
  pipbasfl='MOL_1_3_4.BAS'
  scalemod=-1
  call initpip()

  ! ANN
  ANN%structfl='struct'
  call ANN%init()
  if (ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if (ANN%RY .ne. lastlayer) stop 'ANN%RX .ne. npip'
  !allocate(coef(ncons),coefc(ncons))
  
  inquire(file='NN.txt',exist=ios)
  if (.not. ios) stop "No NN.txt file found."
  print *, "Reading NN.txt..."
  call ANN%savenn('NN.txt',0)
 
  return
end subroutine prepot


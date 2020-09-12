!==============================================================================
program surfgen
  use makesurfdata
  use hddata
  use pip
  implicit none
  real*8, allocatable :: igeom1(:), igeom2(:), x(:,:), x0(:,:)
  real*8, allocatable :: t1(:),t2(:),t3(:),t4(:),t5(:),t6(:)
  real*8, allocatable :: t(:,:),Bmat(:,:),Bmat0(:,:)
  real*8 :: dx
  integer :: i,j,k

  call random_seed()
  call readinput
  call initpip

  allocate(igeom1(ncoord),igeom2(ncoord),x(3,natoms),x0(3,natoms))
  allocate(t1(ncoord),t2(ncoord),t3(ncoord),t4(ncoord),t5(ncoord),t6(ncoord))
  allocate(t(ncoord,6),Bmat(ncoord,3*natoms),Bmat0(ncoord,3*natoms))

  call random_number(x0)

  !E
  x=x0
  call buildWBmat(natoms,ncoord,x,t1,Bmat)
  print*, t1
  print*, '--------------------------------------------------------'

  !(12)
  x=x0
  x(:,2)=x0(:,3)
  x(:,3)=x0(:,2)
  call buildWBmat(natoms,ncoord,x,t2,Bmat)
  print*, t2
  print*, '--------------------------------------------------------'

  !(23)
  x=x0
  x(:,3)=x0(:,4)
  x(:,4)=x0(:,3)
  call buildWBmat(natoms,ncoord,x,t3,Bmat)
  print*, t3
  print*, '--------------------------------------------------------'

  !(13)
  x=x0
  x(:,2)=x0(:,4)
  x(:,4)=x0(:,2)
  call buildWBmat(natoms,ncoord,x,t4,Bmat)
  print*, t4
  print*, '--------------------------------------------------------'

  !(123)
  x=x0
  x(:,2)=x0(:,3)
  x(:,3)=x0(:,4)
  x(:,4)=x0(:,2)
  call buildWBmat(natoms,ncoord,x,t5,Bmat)
  print*, t5
  print*, '--------------------------------------------------------'

  !(132)
  x=x0
  x(:,2)=x0(:,4)
  x(:,3)=x0(:,2)
  x(:,4)=x0(:,3)
  call buildWBmat(natoms,ncoord,x,t6,Bmat)
  print*, t6
  print*, '--------------------------------------------------------'

  t(:,1)=abs(t2-t1)
  t(:,2)=abs(t3-t2)
  t(:,3)=abs(t4-t3)
  t(:,4)=abs(t5-t4)
  t(:,5)=abs(t6-t5)
  t(:,6)=abs(t1-t6)

  print*, maxval(abs(t),2)
  print*, '--------------------------------------------------------'

  dx=1.d-6
  k=0
  do i=1,natoms
    do j=1,3
      k=k+1
      x=x0
      x(j,i)=x0(j,i)-dx
      call buildWBmat(natoms,ncoord,x,igeom1,Bmat0)
      x=x0
      x(j,i)=x0(j,i)+dx
      call buildWBmat(natoms,ncoord,x,igeom2,Bmat0)
      Bmat(:,k)=(igeom2-igeom1)/(2.d0*dx)
    end do
  end do

  call buildWBmat(natoms,ncoord,x0,igeom1,Bmat0)

  print*, maxval(abs(Bmat-Bmat0),2)

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
  namelist /fitting/ npoints, epmax, useGPU, w_energy, w_grad, w_fij, lambda, &
                     ptwthreshold, energyT, highEScale

  !readEner: uniform shift on ab initio energies
  eshift = 56.475986184916d0 !Hartree
  deg_cap=5.d-2 !eV
  gorder=deg_cap !eV

  !progdata
  natoms=4
  printlvl=5

  !hddata
  !ncoord = ANN%RX + 1, the first RX coordinates are PIPs, the last is OOP
  nstates=2
  ncoord=31
  !pip
  npip=ncoord-1
  pipbasfl='MOL_1_3_4.BAS'
  scalemod=3
  morse=2.d0
  !ANN
  ANN%structfl='struct'
  call ANN%init()
  if(ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if(ANN%RY .ne. 3) stop 'ANN%RY .ne. 3'

  !readdisps
  enfptn   =  'energy.all'
  gmfptn   =  'geom.all'
  grdfptn  =  'cartgrd.drt1.state$.all'
  cpfptn   =  'cartgrd_total.drt1.state$.drt1.state$.all'

  !makeLocalIntCoord
  intGradT = 1.d-10
  intGradS = intGradT

  !number of internal degrees
  nvibs=3*natoms-6

  !OrthGH_ab
  deggrdbinding = .true.
  gcutoff = 1.d-8

  call FLUnit(fid)
  open(fid,file='fitting.para',delim='APOSTROPHE')
  read(unit=fid,nml=fitting)
  close(fid)

  call FLUnit(fid)
  open(fid,file='fitting.nml',delim='APOSTROPHE')
  write(unit=fid,nml=fitting)
  close(fid)

  return
end SUBROUTINE readinput
!================================================================================

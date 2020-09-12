!==========================================================================
program main
  use nnpes
  implicit none
  real*8, parameter :: au2cm=219474d0
  character(3) :: sym
  real*8 :: anums
  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: h
  real*8, dimension(:,:,:), allocatable :: geom
  real*8, dimension(:,:,:), allocatable :: cg, dcg
  real*8, allocatable :: energy(:,:)
  real*8 :: r(6)
  integer :: i,j,k,total,ios,idx

  !initialize NN PES
  call pesinit

  total=200
  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  allocate(geom(3,natoms,total))
  allocate(energy(nstates,total))

  open(unit=100,file='geom.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=101,file='energy.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)

  do i=1,total
    do j=1,natoms
      read(100,*,iostat=ios) sym,anums,geom(:,j,i)
    end do
    read(101,*,iostat=ios) energy(:,i)
    if(ios.ne.0) then
      total=i-1
      exit
    end if
  end do

  print*,total, 'valid points.'
  !energy shift
  energy=energy+56.475986184916d0

  open(200,file='energy.dat')
  do i=1,total
    !NN PES evaluation
    call NNEvaluate(geom(:,:,i),e,cg,h,dcg)
    call x2r(4,geom(:,:,i),r)
    if(r(1).gt.8.d0) cycle
    write(200,"(9e20.10)") r(1),energy(:,i)*au2cm,e*au2cm,h*au2cm
  end do

  stop
end
!==========================================================================
subroutine x2r(natoms,x,r)
  !calculate the bond distances
  implicit none
  integer, intent(in) :: natoms
  real*8, intent(in) :: x(3,natoms)
  real*8, intent(out) :: r(natoms*(natoms-1)/2)
  real*8 :: dx(3)
  integer :: i,j,k

  k=0
  do i=1,natoms-1
    do j=i+1,natoms
      k=k+1
      dx=x(:,i)-x(:,j)
      r(k)=sqrt(dot_product(dx,dx))
    end do
  end do

  return
end
!==========================================================================

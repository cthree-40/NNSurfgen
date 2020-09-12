!===================================================================================
!Module for NN PES
module nnpes
  use nn_class
  implicit none
  integer :: natoms
  type(nn) :: ANN !Feed-Forward Neural Network
  integer :: ncoord   !total number of internal coords
  integer :: nstates  !number of electronic states
end module nnpes
!===================================================================================
!Module for PIP
module pip
  use nnpes, only: natoms
  implicit none
  integer :: npip
  integer :: nbond
  integer :: scalemod
  integer :: maxn
  real*8 :: morse
  character(len=99) :: pipbasfl

  type pipstruct
    integer :: n
    integer, allocatable :: bondlist(:,:)
  end type pipstruct

  type(pipstruct), allocatable :: pipdat(:)

  contains
  !---------------------------------------------------------------------------
  subroutine initpip
    implicit none
    integer :: fid,ios
    character(1) :: x
    integer :: id,order,num,i,j,idprev
    integer, allocatable :: list(:)

    nbond=natoms*(natoms-1)/2
    allocate(list(nbond))
    if(allocated(pipdat)) deallocate(pipdat)
    allocate(pipdat(npip))

    call FLUnit(fid)
    open(unit=fid,file=trim(adjustl(pipbasfl)),access='sequential',&
         form='formatted',STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
    if(ios .ne. 0) then
      print *,"Failed to open file [", trim(adjustl(pipbasfl)),"]"
      stop
    end if

    idprev=0
    maxn=-1
    do while(ios.eq.0)
      read(fid,*,iostat=ios) id,order,num,i,x,list
      if(order.eq.0) cycle
      if(id .gt. npip) exit
      if(id.ne.idprev) then
        !new pip term
        idprev=id
        pipdat(id)%n=num
        maxn=max(num,maxn)
        allocate(pipdat(id)%bondlist(num,nbond))
        pipdat(id)%bondlist(i,:)=list
      else
        !the same pip term
        pipdat(id)%bondlist(i,:)=list
      end if
    end do
    if(idprev .ne. npip) stop 'idprev .ne. npip!'
    close(fid)

    return
  end subroutine initpip
  !---------------------------------------------------------------------------
  subroutine evaluatepip(cgeom,pip,pipbmat)
    implicit none
    real*8, intent(in) :: cgeom(3*natoms)
    real*8, intent(out) :: pip(npip),pipbmat(npip,3*natoms)
    real*8, allocatable :: r(:),Bmat0(:,:),Bmat(:,:),mono(:,:),dmono(:,:,:),dt(:)
    real*8 :: dx(3)
    integer :: i,j,k,l,p

    allocate(r(nbond),mono(maxn,nbond),dmono(maxn,nbond,3*natoms),dt(3*natoms))
    allocate(Bmat0(nbond,3*natoms),Bmat(nbond,3*natoms))

    Bmat0=0.d0
    k=0
    do i=1,natoms-1
      do j=i+1,natoms
        k=k+1
        dx=cgeom(3*i-2:3*i)-cgeom(3*j-2:3*j)
        r(k)=sqrt(dot_product(dx,dx))
        Bmat0(k,3*i-2:3*i)=dx/r(k)
        Bmat0(k,3*j-2:3*j)=-Bmat0(k,3*i-2:3*i)
      end do
    end do

    select case(scalemod)
    case(-1)
      !inverse
      r=1.d0/r
      do i=1,nbond
        Bmat(i,:)=-Bmat0(i,:)*r(i)**2
      end do
    case(0)
      !no scale
      Bmat=Bmat0
    case(1)
      !morse
      r=exp(-morse*r)
      do i=1,nbond
        Bmat(i,:)=-Bmat0(i,:)*r(i)*morse
      end do
    case(2)
      !(exp(-ar)/r
      mono(1,1:nbond)=1.d0/r
      r=exp(-morse*r)/r
      do i=1,nbond
        Bmat(i,:)=-Bmat0(i,:)*r(i)*(mono(1,i)+morse)
      end do
    case default
      stop 'Wrong scalemod in PIP!'
    end select

    do i=1,npip
      pip(i)=0.d0
      do j=1,pipdat(i)%n
        do k=1,nbond
          p=pipdat(i)%bondlist(j,k)
          if(p.eq.0) then
            mono(j,k)=1.d0
            dmono(j,k,:)=0.d0
          else if(p.gt.0) then
            mono(j,k)=r(k)**p
            dmono(j,k,:)=Bmat(k,:)*dble(p)*r(k)**(p-1)
          else
            stop 'p.lt.0!'
          end if
        end do
        pip(i)=pip(i)+product(mono(j,:))
      end do
      pipbmat(i,:)=0.d0
      do j=1,pipdat(i)%n
        do k=1,nbond
          dt=dmono(j,k,:)
          do l=1,nbond
            if(l.eq.k) cycle
            dt=dt*mono(j,l)
          end do
          pipbmat(i,:)=pipbmat(i,:)+dt
        end do
      end do
    end do

    return
  end subroutine evaluatepip
  !---------------------------------------------------------------------------
end
!================================================================================
subroutine buildWBmat(natom,ncoords,cgeom,igeom,Bmat)
  use pip
  implicit none
  integer, intent(in) :: natom, ncoords
  real*8, intent(in) :: cgeom(3*natom)
  real*8, intent(out) :: igeom(ncoords)
  real*8, intent(out) :: Bmat(ncoords,3*natoms)
  real*8 :: dcoord(12)

  !PIP
  call evaluatepip(cgeom,igeom(1:npip),Bmat(1:npip,:))

  !f(ncoords)
  !N-H-H-H
  call OOP(cgeom(1:3),cgeom(4:6),cgeom(7:9),cgeom(10:12),igeom(ncoords),dcoord,1)
  Bmat(ncoords,:)=dcoord

  return
end
!===============================================================================
function triple_product(a,b,c)
  implicit none
  real*8 :: triple_product
  real*8, intent(in) :: a(3),b(3),c(3)
  triple_product = a(1)*b(2)*c(3) - a(1)*b(3)*c(2) - a(2)*b(1)*c(3) + &
                   a(2)*b(3)*c(1) + a(3)*b(1)*c(2) - a(3)*b(2)*c(1)
  return
end
!=============================================================================
subroutine OOP(a,b,c,d,coord,dcoord,id)
  implicit none
  real*8, intent(in) :: a(3),b(3),c(3),d(3)
  integer, intent(in) :: id
  real*8, intent(out) :: coord, dcoord(12)
  real*8, external :: triple_product
  real*8 :: ab(3), ac(3), ad(3), rab, rac, rad
  real*8 :: abx(12),acx(12),adx(12),dtp(12)

  ab=a-b
  ac=a-c
  ad=a-d
  rab=sqrt(dot_product(ab,ab))
  rac=sqrt(dot_product(ac,ac))
  rad=sqrt(dot_product(ad,ad))
  coord=triple_product(ab,ac,ad)/(rab*rac*rad)

  if(id.ne.0) then
    abx=0.d0
    abx(1:3)=ab/rab
    abx(4:6)=-abx(1:3)
    acx=0.d0
    acx(1:3)=ac/rac
    acx(7:9)=-acx(1:3)
    adx=0.d0
    adx(1:3)=ad/rad
    adx(10:12)=-adx(1:3)
    !tp = ab(1)*ac(2)*ad(3) - ab(1)*ac(3)*ad(2) - ab(2)*ac(1)*ad(3) + &
    !     ab(2)*ac(3)*ad(1) + ab(3)*ac(1)*ad(2) - ab(3)*ac(2)*ad(1)
    dtp(1)=ac(2)*ad(3)-ac(3)*ad(2)-ab(2)*ad(3)+ab(2)*ac(3)+ab(3)*ad(2)-ab(3)*ac(2)
    dtp(2)=ab(1)*ad(3)-ab(1)*ac(3)-ac(1)*ad(3)+ac(3)*ad(1)+ab(3)*ac(1)-ab(3)*ad(1)
    dtp(3)=ab(1)*ac(2)-ab(1)*ad(2)-ab(2)*ac(1)+ab(2)*ad(1)+ac(1)*ad(2)-ac(2)*ad(1)
    dtp(4)=-ac(2)*ad(3)+ac(3)*ad(2)
    dtp(5)=ac(1)*ad(3)-ac(3)*ad(1)
    dtp(6)=-ac(1)*ad(2)+ac(2)*ad(1)
    dtp(7)=ab(2)*ad(3)-ab(3)*ad(2)
    dtp(8)=-ab(1)*ad(3)+ab(3)*ad(1)
    dtp(9)=ab(1)*ad(2)-ab(2)*ad(1)
    dtp(10)=-ab(2)*ac(3)+ab(3)*ac(2)
    dtp(11)=ab(1)*ac(3)-ab(3)*ac(1)
    dtp(12)=-ab(1)*ac(2)+ab(2)*ac(1)
    dcoord=dtp/(rab*rac*rad)-coord*(abx/rab+acx/rac+adx/rad)
  end if

  return
end
!==================================================================================
!initialize NN PES
subroutine pesinit
  !struct: NN structure definition file
  !WB.txt: Weights and Biases
  !MOL_1_3_4.BAS: PIP definition file
  use nnpes
  use pip
  implicit none

  natoms=4
  nstates=2

  !ncoord = ANN%RX + 1, the first RX coordinates are PIPs, the last is OOP
  ncoord=31

  !pip
  npip=ncoord-1
  pipbasfl='MOL_1_3_4.BAS'
  scalemod=-1
  call initpip

  !ANN
  ANN%structfl='struct'
  call ANN%init()
  if(ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if(ANN%RY .ne. 4) stop 'ANN%RY .ne. 4'
  call ANN%savenn('WB.txt',0)

  return
end subroutine pesinit
!==================================================================================
!Evaluate diabatic Hamiltonian and its gradients w.r.t. internal coordinates
SUBROUTINE EvaluateHd0(igeom,hmat,dhmat)
  !dhmat: gradients with respect to internal coordinate
  use nnpes
  IMPLICIT NONE
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: hmat(nstates,nstates)
  real*8, intent(out) :: dhmat(ncoord,nstates,nstates)
  real*8 :: oop
  integer :: R

  R=ANN%RX

  !hmat
  oop=igeom(ncoord)
  call ANN%output(igeom(1:R))
  hmat(1,1) = ANN%y(1)
  hmat(1,2) = ANN%y(2)*oop + ANN%y(3)*oop**3
  hmat(2,1) = hmat(1,2)
  hmat(2,2) = ANN%y(4)

  !dhmat
  call ANN%cal_dadx

  dhmat(1:R,1,1)=ANN%dydx(:,1)
  dhmat(ncoord,1,1)=0.d0

  dhmat(1:R,1,2)=ANN%dydx(:,2)*oop + ANN%dydx(:,3)*oop**3
  dhmat(ncoord,1,2)=ANN%y(2) + 3.d0*oop**2*ANN%y(3)

  dhmat(:,2,1)=dhmat(:,1,2)

  dhmat(1:R,2,2)=ANN%dydx(:,4)
  dhmat(ncoord,2,2)=0.d0

  return
END SUBROUTINE EvaluateHd0
!==================================================================================
!NN PES evaluation subroutine
!atomic units bohr and hartree are used
!geom: cartesian coordinates
!e: adiabatic energies
!h: diabatic Hamiltonian
!cg: energy gradients and derivative couplings
!dcg: gradients of diabatic Hamiltonian
subroutine NNEvaluate(geom,e,cg,h,dcg)
  use nnpes
  use pip
  implicit none
  real*8, parameter :: ev2au=1.d0/27.21138602d0
  real*8, intent(in) :: geom(3,natoms)
  real*8, intent(out) :: e(nstates)
  real*8, intent(out) :: h(nstates,nstates)
  real*8, intent(out) :: cg(3*natoms,nstates,nstates)
  real*8, intent(out) :: dcg(3*natoms,nstates,nstates)
  real*8, allocatable :: igeom(:),dh(:,:,:),Bmat(:,:),ckl(:,:),work(:)
  real*8 :: oop,de
  integer :: i,j,info,lwork

  lwork=999
  allocate(igeom(ncoord),dh(ncoord,nstates,nstates),Bmat(ncoord,3*natoms))
  allocate(ckl(nstates,nstates),work(lwork))

  call buildWBmat(natoms,ncoord,geom,igeom,Bmat)
  call EvaluateHd0(igeom,h,dh)

  h=h*ev2au
  dh=dh*ev2au

  do i=1,nstates
     do j=i,nstates
       call dgemv('T',ncoord,3*natoms,1.d0,Bmat,ncoord,dh(:,i,j),1,0.d0,dcg(:,i,j),1)
       if(j.ne.i) dcg(:,j,i)=dcg(:,i,j)
     end do
  end do

  ckl=h
  call dsyev('V','U',nstates,ckl,nstates,e,work,lwork,info)
  if(info.ne.0) stop 'Failed to solve eigenvectors!'

  do i=1,3*natoms
    cg(i,:,:)=matmul(transpose(ckl),matmul(dcg(i,:,:),ckl))
  end do

  !derivative coupling
  do i=1,nstates-1
    do j=i+1,nstates
      de=e(j)-e(i)
      if(abs(de).lt.1d-30) de=1d-30
      cg(:,i,j)=cg(:,i,j)/de
      cg(:,j,i)=-cg(:,i,j)
    end do
  end do

  return
end
!==================================================================================

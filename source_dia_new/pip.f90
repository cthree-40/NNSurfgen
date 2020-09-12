!=============================================================================
module pip
  use progdata, only: natoms
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

    call FLUnit(fid)
    open(fid,file='piplist')
    do i=1,npip
      do j=1,pipdat(i)%n
        write(fid,"(<nbond+3>i4)") i,pipdat(i)%n,j,pipdat(i)%bondlist(j,:)
      end do
    end do
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
!=============================================================================

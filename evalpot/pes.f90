!===================================================================================
!Module for NN PES
module nnpes
  use nn_class
  implicit none
  integer :: natoms
  type(nn) :: ANN !Feed-Forward Neural Network
  integer :: ncoord   !total number of internal coords
  integer :: nstates  !number of electronic states

  ! Coordinate definition type
  type, public :: cdef
     integer :: ctype
     integer :: nperm
     integer :: alist(4,6)
     integer :: acoef(6)   ! coefficient for atom permlist
     integer :: scaling    ! scaling type
     double precision :: sccoef(2)  ! scaling coefs
  end type cdef
  ! Form of Hd
  type, public :: hdform
     integer   :: nc
     integer   :: terms(4,10)
  end type hdform

  type(cdef),  dimension(:),   allocatable :: coord ! Coordinates
  type(hdform),dimension(:,:), allocatable :: hddef ! Form of Hd

  integer :: lastlayer   ! Neurons in last layer
  integer :: noopc       ! Number of OOP coordinates defined in file  

  double precision, dimension(:,:), allocatable :: evecstore

  character(3), dimension(:), allocatable :: atomlabels
  
  integer :: EvalCount
  integer :: NEval
  
  ! parsing of trajectory variables
  logical :: parsing
  logical :: calcmind
  logical :: paused
  integer :: GUNIT
  double precision :: timetraj
  integer          :: isurftraj
  ! number of reference points
  integer          :: num_rpts
  double precision, dimension(:,:), allocatable :: refpoints
  ! distance coordinates
  integer                :: ndcoord
  integer, dimension(20) :: dcoordls
  integer                :: NeighborID
  double precision, dimension(:), allocatable :: lastrgeom
  double precision, dimension(:), allocatable :: ldbounds, udbounds
  logical :: firstpt
contains
  !----------------------------------------------------------------------
  subroutine readcoordinput ()
    implicit none
    integer :: ncoords
    character(255), dimension(:), allocatable :: title
    integer :: i, j, k, l, ios
    open(file="coord.input",unit=11,status="old",action="read", iostat=ios)
    read(unit=11,fmt=1001) ! comment line
    read(unit=11,fmt=1000,iostat=ios) noopc
    ncoords = noopc
    if (ios .ne. 0) then
      print *, "Error reading coord.input file."
    end if
    print *, "Reading ", ncoords, " coordinate definitions."
    if (.not. allocated(coord)) allocate(coord(ncoords))
    allocate(title(ncoords))
    do i = 1, ncoords
      read(unit=11,fmt=1001,iostat=ios) title(i)
      read(unit=11,fmt=1002,iostat=ios) &
        coord(i)%ctype, coord(i)%scaling, coord(i)%nperm
      if (coord(i)%ctype .eq. -4) then
        print "(A,i2,A)", "Coordinate ", i, " is a product coordinate."
      end if
      do j = 1, coord(i)%nperm
        read(unit=11,fmt=1003,iostat=ios) coord(i)%alist(:,j), coord(i)%acoef(j)
      end do
      read(unit=11,fmt=1004) coord(i)%sccoef(1), coord(i)%sccoef(2)
    end do
    print *, "Coordinate definitions: "
    do i = 1, ncoords
      print "(A)", trim(adjustl(title(i)))
      print *, coord(i)%ctype, coord(i)%nperm
      print "('Scaling info: (type, coef(1), coef(2))',i3,f6.2,f6.2)", &
        coord(i)%scaling, coord(i)%sccoef(1), coord(i)%sccoef(2)
      do j = 1, coord(i)%nperm
        print "(A,4i3,A,i3)", " Atoms: ", coord(i)%alist(:,j), &
          "  Coefficient: ", coord(i)%acoef(j)
      end do
    end do
    ! Get state info (for block number)
    read(unit=11,fmt=1000,iostat=ios) nstates
    print *, " NSTATES = ", nstates
    if (.not. allocated(hddef)) allocate(hddef(nstates,nstates))
    do i = 1, nstates
      do j = i, nstates
        read(unit=11,fmt=1000,iostat=ios)
        read(unit=11,fmt=1000,iostat=ios) hddef(j,i)%nc
        if (hddef(j,i)%nc .gt. 10) stop "nc > 10"
        do k = 1, hddef(j,i)%nc
          read(unit=11,fmt=1003,iostat=ios) hddef(j,i)%terms(:,k)
          hddef(i,j) = hddef(j,i)
        end do
      end do
    end do
    print *, "Block by block defs: "
    do i = 1, nstates
      do j = i, nstates
        print "('Block(',i2,',',i2,'):')", j, i
        print "('Number of coords: ',i4)", hddef(j,i)%nc
        if (hddef(j,i)%nc .lt. 0) then
          do k = 1, abs(hddef(j,i)%nc)
            print *, " PIP "
          end do
        else
          do k = 1, hddef(j,i)%nc
            print "(4i5)", hddef(j,i)%terms(:,k)
          end do
        end if
      end do
    end do
    
    ! Compute lastlayer size
    lastlayer = 0
    do i = 1, nstates
      do j = i, nstates
        if (hddef(j,i)%nc .lt. 0) then
          ! PIP
          do k = 1, abs(hddef(j,i)%nc)
            lastlayer = lastlayer + 1
          end do
        else
          do k = 1, hddef(j,i)%nc
            do l = 1, 4
              if (hddef(j,i)%terms(l,k) .ne. 0) then
                ! Term
                lastlayer = lastlayer + 1
              end if
            end do
          end do
        end if
      end do
    end do
    close(unit=11)
    return
    
1000 format(i10)      ! first line
1001 format(A)        ! coordinate titles
1002 format(3i6)      ! ctype, scaling type, perms
1003 format(4i6,i6)   ! atom list defs, coefficient
1004 format(2f6.3)    ! scaling coeficients
  end subroutine readcoordinput
  
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
end module pip
!=============================================================================
! evaluateoopcoord: evaluate user defined oop coordinates
! Input:
!  n     = coordinate number
!  natm  = number of atoms
!  cgeom = geometry
! Output:
!  val   = value of coordinate
!  dval  = derivative of coordinate
! This is a recursive subroutine to evaluate product coordinates.
recursive subroutine evaluateoopcoord(n, natm, cg, val, dval)
  use nnpes
  implicit none
  integer, intent(in) :: n, natm
  double precision, intent(in)  :: cg(natm*3)
  double precision, intent(out) :: val, dval(natm*3)

  double precision :: v, dv(natm*3)
  double precision :: v1, dv1(natm*3), v2, dv2(natm*3)
  integer :: a1, a2, a3, a4
  integer :: c(2)
  integer :: i, j, offs

  val  = 0d0
  dval = 0d0

  select case(coord(n)%ctype)

  case(-2)

     ! Umbrella OOP coordinate
     do i = 1, coord(n)%nperm

        call oop3(natm, coord(n)%alist(:,i), cg, v, dv, &
             coord(n)%sccoef(1), coord(n)%scaling, coord(n)%sccoef(2))
        val  = val  + v*coord(n)%acoef(i)

        do j = 1, 4
            offs = coord(n)%alist(j,i)
            dval(offs*3-2:offs*3) = dval(offs*3-2:offs*3) + dv(j*3-2:j*3)*coord(n)%acoef(i)
        end do
        
     end do

     val = val / coord(n)%nperm
     dval=dval / coord(n)%nperm
     
  case(-3)

     ! Four-center dot product
     do i = 1, coord(n)%nperm

        call dot4(natm, coord(n)%alist(:,i), cg, v, dv, &
             coord(n)%sccoef, coord(n)%scaling)

        val  = val  + v*coord(n)%acoef(i)

        do j = 1, 4
            offs = coord(n)%alist(j,i)
            dval(offs*3-2:offs*3) = dval(offs*3-2:offs*3) + dv(j*3-2:j*3)*coord(n)%acoef(i)
        end do

    end do
     
     val  = val  / coord(n)%nperm
     dval = dval / coord(n)%nperm

  case(-4)
     
     ! Product coordinate
     ! Set products
     c(1) = coord(n)%alist(1,1)
     c(2) = coord(n)%alist(1,2)

     ! Evaluate first coordinate
     call evaluateoopcoord(c(1), natm, cg, v1, dv1)
     ! Evaluate second coordinate
     call evaluateoopcoord(c(2), natm, cg, v2, dv2)

     val = v1 * v2
     dval = v1 * dv2(:) + v2 * dv1(:)
     
  case default
     print *, "Coordtype = ", coord(n)%ctype
     stop "Coordinate type not defined!"
  end select
     
  return
end subroutine evaluateoopcoord
!*******************************************************************************
! 4 center dot product with scalings 
! f = (a1-a2).(a3-a4) / norm(a1-a2) /norm(a3-a4)
! scalings : =0 no scaling
!            >0 multiply the product of 4 scaled distances d13,d14,d23,d24
!            <0 multiply the harmonic mean of 4 scaled distances d13,d14,d23,d24
!*******************************************************************************
SUBROUTINE dot4(natoms,atms,cgeom,fval,bval,coef,scal)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,scal
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(3,natoms),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: fval
  DOUBLE PRECISION,DIMENSION(2),INTENT(IN)        :: coef
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)       :: bval
  
  double precision  :: f, r12,r34  ! unscaled dot prod and the two distances 
  double precision,dimension(12) ::b, w12,w34 ! gradients of f, r12 and r34
  double precision,dimension(3)  ::d12, d34  !d12=a1-a2, d34=d3-d4
  double precision  :: w(6) , dwdR(12,6)  , prod
  double precision  :: coef2(2)
  integer           :: i,j,k
  
  d12 = cgeom(:,atms(1))-cgeom(:,atms(2))
  d34 = cgeom(:,atms(3))-cgeom(:,atms(4))
  
  f = dot_product(d12,d34)
  b(1:3)   = d34
  b(4:6)   =-d34
  b(7:9)   = d12
  b(10:12) =-d12
  
  w12 = 0d0
  coef2(1)  = 1d0
  coef2(2)  = 0d0
  call calcwij(0,atms(1),atms(2),coef2,cgeom,r12,w12)
  if(abs(r12)<1D-30) r12=sign(1D-30,r12)
  f = f/r12
  b = b/r12 -f/r12*w12
  
  w34 = 0d0
  call calcwij(0,atms(3),atms(4),coef2,cgeom,r34,w34)
  if(abs(r34)<1D-30) r34=sign(1D-30,r34)
  f = f/r34
  b = b/r34
  b(7:12) = b(7:12) - f/r34*w34(1:6)
  
  k=1
  dwdR = 0d0
  if(scal/=0)then
     do i=1,3
        do j=i+1,4
           call calcwij(abs(scal),atms(i),atms(j),coef,cgeom,w(k),w12)
           if(abs(w(k))<1d-30)w(k)=sign(1D-30,w(k))
           dwdR(i*3-2:i*3,k) = w12(1:3) 
           dwdR(j*3-2:j*3,k) = w12(4:6) 
           k=k+1
        end do!j
     end do!i
  end if
  
  ! perform scaling
  !  0   : unscaled
  ! >0   : scaling using products of scaled dists 
  ! <0   : scaling using harmonic mean of scaled dists
  select case(scal)
  case(0)
     fval = f
     bval = b
  case(1:)
     prod = product(w)
     fval = f*prod
     bval = b*prod +  fval*  &
          (dwdR(:,1)/w(1)+dwdR(:,2)/w(2)+dwdR(:,3)/w(3)+dwdR(:,4)/w(4)&
          +dwdR(:,5)/w(5)+dwdR(:,6)/w(6))
  case(:-1)
     prod = 1/sum(1/w)
     fval = 4*f*prod
     bval = 4*b*prod + fval*prod*&
          (dwdR(:,1)/w(1)**2+dwdR(:,2)/w(2)**2+dwdR(:,3)/w(3)**2 &
          +dwdR(:,4)/w(4)**2+dwdR(:,5)/w(5)**2+dwdR(:,6)/w(6)**2)
  case default
  end select! case(scal)
END SUBROUTINE dot4
!*******************************************************************************
! SUBROUTINE oop3 
! OOP angle values and derivatives, with the expoential scaling method
! f= C1*g*Product[Scale[rij]]
! where a is the scaling factor, g is the triple product of unit displacement vectors
! here, one of the four atoms is made a special vertex and only the three distances
! launching from this atom is used for scaling.  This the three atom atoms are
! permutationally equivalent while this one is special.
! ******************************************************************************
! natom s [input] INTEGER
!         Total number of atoms
! atms    [input] INTEGER, dimension(4) 
!         The indices of the four reference atoms
! cgeom   [input] DOUBLE PRECISION, dimension(3*natoms)
!         Cartesian geometry of all the atoms
! qval    [output] DOUBLE PRECISION
!         Value of the scaled out-of-plane coordinate
! bval    [output] DOUBLE PRECISION, dimension(12)
!         Derivatives of the OOP coordinate with respect to the cartesian coordinate
! scale   [input] DOUBLE PRECISION
!         Scaling factor used in the exponetials 
! sctype  [input] INTEGER
!         Type of the scaling function used.
SUBROUTINE oop3(natoms,atms,cgeom,qval,bval,scale,sctype,factor)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,sctype
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(natoms*3),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: qval
  DOUBLE PRECISION,INTENT(IN)                     :: scale
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)      :: bval
  double precision,intent(IN)                     :: factor 
  
  INTEGER           :: i, j, sgns(6,4)
  DOUBLE PRECISION  :: dnorm2(6), dvecs(6,3), s(3),  pw, tdv(3,6),  &
       geom(4,3), ddvec(3,3), ss,    dw(12),    &
       coef(2), dwdR(3,12), w(3), dwdR2(3,12), w2(3)
  ! eprod = product of the scaling exponentials
  ! fval  = value of unscaled oop coordinate (the triple product)
  ! fbmat = derivatives of unscaled OOP coordinate
  DOUBLE PRECISION  :: eprod, fval, fbmat(12),a,p
  
  ! The oriented direction vectors between all 6 atom pairs
  sgns(1,:)=(/ 1, 0, 0,-1/)
  sgns(2,:)=(/ 0, 1, 0,-1/)
  sgns(3,:)=(/ 0, 0, 1,-1/)
  sgns(4,:)=(/ 1, 0,-1, 0/)
  sgns(5,:)=(/ 0, 1,-1, 0/)
  sgns(6,:)=(/ 1,-1, 0, 0/)
  
  ! calculate the three scaling factors and their derivatives
  ! scaled distances for the three bonds A1-A2, A1-A3 and A1-A4 are used
  dwdR = 0d0
  coef(1)   = 1d0
  coef(2)   = 0d0
  do i=2,4
     ! scaling mode=0   (unscaled).  unscaled distances are produced here
     call calcwij(int(0),atms(1),atms(i),coef,cgeom,w(i-1),dw)
     if(abs(w(i-1))<1D-30) w(i-1)=sign(1D-30,w(i-1))
     dwdR(i-1,1:3) = dw(1:3)
     dwdR(i-1,i*3-2:i*3) = dw(4:6)
  end do
  
  ! here are the scaled distances
  coef(1) =scale
  coef(2) =factor
  dwdR2 = 0d0
  if(sctype.ne.0)then
     do i=2,4
        call calcwij(abs(sctype),atms(1),atms(i),coef,cgeom,w2(i-1),dw)
        if(abs(w2(i-1))<1D-50) w2(i-1)=sign(1D-50,w2(i-1))
        dwdR2(i-1,1:3)= dw(1:3)
        dwdR2(i-1,i*3-2:i*3)= dw(4:6)
     end do
  end if
  
  ! extract the geometry of the 4 involved atoms
  do i=1,4
     geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do
  
  ! calculate displacement vectors between these atoms
  dvecs=matmul(sgns,geom)
  do i=1,6
     dnorm2(i)=dot_product(dvecs(i,:),dvecs(i,:))
  end do
  
  ! calculate value of unscaled OOP coordinate and its derivatives
  fval=det3(dvecs)
  do i=1,3
     do j=1,3
        ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
     end do
  end do
  tdv=transpose(dvecs)
  fbmat(1:3)=cross(tdv(:,2),tdv(:,3))
  fbmat(4:6)=cross(tdv(:,3),tdv(:,1))
  fbmat(7:9)=cross(tdv(:,1),tdv(:,2))
  fbmat(10:12)=-fbmat(1:3)-fbmat(4:6)-fbmat(7:9)
  
  pw = product(w)
  fval=fval/pw
  fbmat=fbmat/pw-fval*(dwdR(1,:)/w(1)+dwdR(2,:)/w(2)+dwdr(3,:)/w(3))
  
  ! calculate the scaled oop
  if(sctype.lt.0)then
     ! use the harmonic mean of scaled distances to scale OOP
     s(1) = w2(2)*w2(3)
     s(2) = w2(3)*w2(1)
     s(3) = w2(1)*w2(2)
     ss = sum(s)
     a = w2(1)*s(1)/ss
     qval= fval*a
     bval=fbmat*a+fval*(s(1)**2*dwdR2(1,:)+s(2)**2*dwdR2(2,:)+s(3)**2*dwdR2(3,:))/ss**2
  else if(sctype.eq.0)then
     ! reciprocal or order scaling
     do i=1,3
        if(abs(w(i))<1D-20)w(i)=sign(1D-20,w(i))
     end do
     p=coef(2)**3/product(w)
     a=p**coef(1)
     qval = a*fval
     bval = a*fbmat-coef(1)*qval*(dwdR(1,:)/w(1)+dwdR(2,:)/w(2)+dwdR(3,:)/w(3))
  else
     ! multiply all 3 scaled length scaling functions
     eprod = product(w2)
     qval = fval*eprod
     bval = fbmat*eprod+fval*(w2(2)*w2(3)*dwdR2(1,:)+w2(1)*w2(3)*dwdR2(2,:)+w2(1)*w2(2)*dwdR2(3,:))
  end if
  
CONTAINS
  
  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m)
    IMPLICIT NONE
    double precision ::  det3, m(:,:)
    det3=(m(1,2)*m(2,3)-m(1,3)*m(2,2))*m(3,1)+&
         (m(1,3)*m(2,1)-m(1,1)*m(2,3))*m(3,2)+&
         (m(1,1)*m(2,2)-m(1,2)*m(2,1))*m(3,3)
  END FUNCTION det3
  
  !this function calculate the cross product of two 3D vectors
  FUNCTION cross(v1,v2)RESULT(v3)
    implicit none
    double precision,dimension(3) :: v1,v2,v3
    v3(1)=-v1(3)*v2(2)+v1(2)*v2(3)
    v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
    v3(3)=-v1(2)*v2(1)+v1(1)*v2(2)
  END FUNCTION cross
  
END SUBROUTINE oop3
SUBROUTINE dvec(na,a1,a2,cgeom,vec,vnorm)
  IMPLICIT NONE
  INTEGER, INTENT(IN)                         :: na,a1,a2
  DOUBLE PRECISION,INTENT(IN)                 :: cgeom(3*na)
  DOUBLE PRECISION,INTENT(INOUT)              :: vnorm
  DOUBLE PRECISION,dimension(3),INTENT(INOUT) :: vec

  INTEGER                                     :: i
  DOUBLE PRECISION                            :: denom

  vnorm = dble(0)
  do i = 1,3
   vec(i) = cgeom(3*(a1-1)+i)-cgeom(3*(a2-1)+i)
   vnorm = vnorm + vec(i)**2
  enddo

  vnorm = Sqrt(vnorm)
  if(vnorm<1d-30)then
   denom = dble(1)
  else
   denom = vnorm
  endif

  vec=vec/denom

  return
end SUBROUTINE dvec

!===============================================================================
!********************************************************************************
! CalcWij  
! Purpose
!         Calculate the value and gradient of one scaled Rij coordinate
! Arguments
! scaling  [input] INTEGER
!          Type of coordinate scaling method
! a1,a2    [input] INTEGER
!          Indices of end point atoms that defines the distance
! coef     [input] DOUBLE PRECISION,dimension(2)
!          Scaling coefficients
! cgeom    [input] 
SUBROUTINE calcwij(scaling,a1,a2,coef,cgeom,w,dwdR)
  use nnpes, only: natoms
  IMPLICIT NONE
  INTEGER,INTENT(IN)                                :: scaling,a1,a2
  DOUBLE PRECISION,dimension(2),INTENT(IN)          :: coef
  DOUBLE PRECISION,dimension(3*natoms),INTENT(IN)   :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                      :: w
  DOUBLE PRECISION,dimension(12),INTENT(OUT)        :: dwdR

! fval     Value of current coordinate in its unscaled form
! fbmat    Cartesian gradients of current coordinate in its unscaled form
! bval     Derivative of a unscaled coordinate with respect to cartesian 
!            coordinates of the four reference atoms.

  DOUBLE PRECISION     ::  fval, bval(12), g

!---Plain or scaled Rij coordinates--->>>

  call stre(natoms,a1,a2,cgeom,fval,bval)
  if(scaling<0) stop  "Invalid scaling type in CalcWij"
  select case(scaling)
    !  no scalings
    case(0)
      w    = (fval-coef(2))*coef(1)
      dwdR = bval*coef(1)

    !  Morse functions   Exp(-c1*(r-c2))
    case(1)
      w    =  exp(-coef(1)*(fval-coef(2)))
      dwdR = -bval*w*coef(1)

    !  Gaussian functions Exp(-c1*(r-c2)**2)
    case(2)
      w    =  exp(-coef(1)*(fval-coef(2))**2)
      dwdR = -bval*w*coef(1)*2*(fval-coef(2))

    !  Screened Columb potential Exp(-c1*(r-c2))/r     (Yukawa, leading term)
    case(3)
      g    =  exp(-coef(1)*(fval-coef(2)))
      w    =  g/(fval+1D-40)
      dwdR = -bval*w*(coef(1)+1/(fval+1D-40))

    !  Long range term of screened Columb w=Exp[c1*(c2-x)]*(x/c2)^(c1*c2)
    !  its derivative: w*c1*(c2-x)/x*x'
    case(4)
      w    = exp(coef(1)*(coef(2)-fval))*(fval/coef(2))**(coef(1)*coef(2))
      dwdR = bval*w*coef(1)*(coef(2)-fval)/fval

    !  Lennard Jones functions (c2/r)**c1
    case(5)
      w    =  (coef(2)/fval)**coef(1)
      dwdR = -bval*w*coef(1)/fval

    !  Shifted(chasmless) Yukawa exp(-c1*(r-c2))/(r+c2)
    case(6)
      g    =  exp(-coef(1)*(fval-coef(2)))
      w    =  g/(fval+coef(2))
      dwdR = -bval*w*(coef(1)+1/(fval+coef(2)))

    ! segmentation function tanh
    case(7)
      w    = tanh((fval-coef(2))/coef(1))/2
      dwdR = bval/(cosh((fval-coef(2))/coef(1))**2*coef(1))/2

    case default
      print *,"scaling = ",scaling
      stop "Unsupported bond distance scaling method."

  end select !case(scaling)
  if(abs(w)<1d-30)w=1d-30
END SUBROUTINE calcwij

!This subroutine calculates the bond length between two reference atoms and the derivative of
! this length with respect to the cartesian coordinates of the two atoms.
SUBROUTINE stre(na,a1,a2,cgeom,qval,bval)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                           :: na,a1,a2
  DOUBLE PRECISION,INTENT(IN)                  :: cgeom(3*na)
  DOUBLE PRECISION,INTENT(INOUT)               :: qval
  DOUBLE PRECISION,dimension(12),INTENT(INOUT) :: bval
  INTEGER                                      :: i
  DOUBLE PRECISION,dimension(3)                :: vec1

  call dvec(na,a1,a2,cgeom,vec1,qval)
  do i = 1,3
   bval(i)   = vec1(i)
   bval(3+i) = -vec1(i)
   bval(6+i) = 0D0
   bval(9+i) = 0D0
  enddo

  return
end SUBROUTINE stre


!================================================================================
subroutine buildWBmat(natom,ncoords,cgeom,igeom,Bmat)
  use pip
  use nnpes
  implicit none
  integer, intent(in) :: natom, ncoords
  real*8, intent(in) :: cgeom(3*natom)
  real*8, intent(out) :: igeom(ncoords)
  real*8, intent(out) :: Bmat(ncoords,3*natoms)

  real*8 :: dcoord(12)
  real*8 :: ccoord
  
  integer :: i
  !PIP
  call evaluatepip(cgeom,igeom(1:npip),Bmat(1:npip,:))

  !f(ncoords)
  ! User defined coordinates from coord.input file.
  do i = 1, noopc
      call evaluateoopcoord(i, natom, cgeom, ccoord, dcoord)
      igeom((ncoord-noopc) + i) = ccoord
      Bmat((ncoord-noopc) + i,:) = dcoord
  end do

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
  integer :: ios

  natoms=4
  nstates=3
  call readcoordinput()
  ncoord = 30 + noopc
  print *, "Total number of coordinates: ", ncoord
  print *, "Number of user defined OOP coordinates: ", noopc

  call readpotentialinput()
  
  !pip
  npip=ncoord-noopc
  pipbasfl='MOL_1_3_4.BAS'
  !scalemod=-1
  call initpip

  !ANN
  ANN%structfl='struct'
  call ANN%init()
  if(ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if(ANN%RY .ne. lastlayer) stop 'ANN%RY .ne. lastlayer'
  inquire(file='NN.txt',exist=ios)
  if (.not. ios) stop "No NN.txt file found."
  print *, "Reading NN.txt..."
  call ANN%savenn('NN.txt',0)

  return
end subroutine pesinit
!==================================================================================
!initialize NN PES (same as pesinit)
subroutine prepot
  !struct: NN structure definition file
  !WB.txt: Weights and Biases
  !MOL_1_3_4.BAS: PIP definition file
  use nnpes
  use pip
  implicit none
  integer :: ios

  natoms=4
  nstates=3
  if (allocated(atomlabels)) deallocate(atomlabels)
  allocate(atomlabels(natoms))
  atomlabels(1) = "O"
  atomlabels(2) = "H"
  atomlabels(3) = "H"
  atomlabels(4) = "H"

  call readcoordinput()
  ncoord = 30 + noopc
  print *, "Total number of coordinates: ", ncoord
  print *, "Number of user defined OOP coordinates: ", noopc

  call readpotentialinput()
  
  !pip
  npip=ncoord-noopc
  pipbasfl='MOL_1_3_4.BAS'
  !scalemod=-1
  call initpip

  !ANN
  ANN%structfl='struct'
  call ANN%init()
  if(ANN%RX .ne. npip) stop 'ANN%RX .ne. npip'
  if(ANN%RY .ne. lastlayer) stop 'ANN%RY .ne. lastlayer'
  inquire(file='NN.txt',exist=ios)
  if (.not. ios) stop "No NN.txt file found."
  print *, "Reading NN.txt..."
  call ANN%savenn('NN.txt',0)

  ! Set variables for trajectory evaluations
  if (allocated(evecstore)) deallocate(evecstore)
  allocate(evecstore(nstates,nstates))
  paused = .false.
  parsing = .true.
  GUNIT = -1
  if (parsing) call openTrajFile(0)
  call ResetEvalCount()
  call inittraj()

  
  calcmind = .true.
  if (calcmind) then
     ndcoord = 6
     dcoordls(1) = 1
     dcoordls(2) = 2
     dcoordls(3) = 3
     dcoordls(4) = 4
     dcoordls(5) = 5
     dcoordls(6) = 6
     call loadrefgeoms()
     call initmindist()
  end if
  return
end subroutine prepot
!==========================================================================
subroutine readpotentialinput
  use nnpes
  use pip
  implicit none
  character(255) :: filename
  integer        :: fileunit
  integer        :: ios
  namelist /general/ scalemod, morse
  filename = "potential.input"
  scalemod=-1
  morse=1.0d0
  call FLUnit(fileunit)
  open(file=filename,unit=fileunit,status="unknown",action="read",iostat=ios)
  if (ios .ne. 0) then
     print *, "Could not open 'potential.input' to read &general."
     print *, " WARNING: Continuing with default values."
     return
  end if
  read(unit=fileunit,nml=general)
  close(unit=fileunit)
  print *, " scalemod: ", scalemod
  print *, " morse:    ", morse
  return
end subroutine readpotentialinput
!==========================================================================
subroutine initmindist
  use nnpes
  implicit none
  integer :: i
  if (allocated(lastrgeom)) deallocate(lastrgeom)
  if (allocated(ldbounds))  deallocate(ldbounds)
  if (allocated(udbounds))  deallocate(udbounds)
  allocate(lastrgeom(ncoord))
  allocate(ldbounds(num_rpts))
  allocate(udbounds(num_rpts))
  lastrgeom = dble(0)
  do i = 1, num_rpts
     ldbounds(i) = dble(0)
     udbounds(i) = dble(0)
  end do
  firstpt = .true.
  NeighborID = 0
end subroutine initmindist
!==========================================================================
! initialze the traj calcs                                                        
SUBROUTINE inittraj()
  use nnpes
  IMPLICIT NONE                                                          
  timetraj  = -1D0                                                       
  isurftraj = 0
END SUBROUTINE inittraj
!=====================================================================
! Reset the evaluation counter                                       
SUBROUTINE ResetEvalCount()                                        
  EvalCount = 0                                                   
END SUBROUTINE ResetEvalCount
!=====================================================================
! Parsing subroutines
SUBROUTINE pauseParsing()                                         
  use nnpes
  paused = parsing .or. paused                                    
  parsing = .false.                                               
END SUBROUTINE                                                    
SUBROUTINE resumeParsing()                                        
  use nnpes
  parsing = parsing .or. paused                                   
  paused = .false.                                                
END SUBROUTINE                                                    
SUBROUTINE enableParsing()                                        
  use nnpes
  parsing = .true.                                                
END SUBROUTINE                                                    
SUBROUTINE disableParsing()                                       
  use nnpes
  parsing = .false.                                               
END SUBROUTINE                                                    
!======================================================================
! Retrieve the evaluation counter                                    
SUBROUTINE GetEvalCount(count)                                     
  INTEGER,INTENT(OUT)  :: count                                   
  count= EvalCount                                                
END SUBROUTINE GetEvalCount
!========================================================================
! loadrefgeoms: load refgeom file
subroutine loadrefgeoms()
  use nnpes
  use pip
  implicit none
  integer, parameter :: maxpts = 15000
  double precision, dimension(3*natoms, maxpts) :: cgeom
  double precision, dimension(ncoord) :: igeom
  double precision, dimension(ncoord, 3*natoms) :: bmat
  integer :: fid, ios
  integer :: i, j
  
  call FLUnit(fid)
  open(unit=fid,file="refgeom",status="old",action="read",access="sequential",&
       form="formatted",iostat=ios)
  if (ios .ne. 0) stop "Could not open refgeom."
  ! Get number of points
  do
     read(unit=fid,fmt=*,iostat=ios)
     if (ios .ne. 0) exit
     num_rpts = num_rpts + 1
  end do
  num_rpts = num_rpts / (natoms + 1)
  if (num_rpts .gt. maxpts) stop "Points in refgeom > 15000!"
  print "(A,i7,A)", "Reading ", num_rpts, " points"
  ! Read reference points
  rewind(unit=fid,iostat=ios)
  if (ios .ne. 0) stop "Rewind failed."
  read(unit=fid,fmt=*,iostat=ios) cgeom
  if (ios .gt. 0) then
     print *, "ios = ", ios
     stop "Failed to read refgeom."
  end if
  close(unit=fid)
  
  ! Allocate refpoints and process internal coordinates
  if (allocated(refpoints)) deallocate(refpoints)
  allocate(refpoints(ncoord, num_rpts))
  do i = 1, num_rpts
     call buildWBMat(natoms, ncoord, cgeom(:,i), igeom, bmat)
     refpoints(:,i) = igeom
  end do
  
  return
end subroutine loadrefgeoms

!========================================================================
! create the files corresponding to a new trajectory
SUBROUTINE openTrajFile(itraj)
  use nnpes
  use pip
  IMPLICIT NONE            
  INTEGER,intent(IN)     ::  itraj                                                      
  integer cnt                                                                    
  character(4)   ::  str                                                               
  
  write(str,"(I4)")itraj                                               
  ! prepare output for geometry, energy and distance print out
  if(GUNIT<0)then
     GUNIT = 427
     close(GUNIT)
  end if
  open(GUNIT,file='trajdata'//trim(adjustl(str))//'.csv',status='REPLACE',action='write')           
  write(str,"(I4)")natoms
  write(GUNIT,&
       "(I5,',',I5,"//trim(adjustl(str))//"(',',A3))")natoms,nstates,atomlabels(1:natoms)
  
END SUBROUTINE openTrajFile
!==================================================================================
!Evaluate diabatic Hamiltonian and its gradients w.r.t. internal coordinates
SUBROUTINE EvaluateHd0(igeom,hmat,dhmat)
  !dhmat: gradients with respect to internal coordinate
  use nnpes
  IMPLICIT NONE
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: hmat(nstates,nstates)
  real*8, intent(out) :: dhmat(ncoord,nstates,nstates)
  real*8 :: oop(noopc) ! user-defined oop coordinates
  integer :: R
  integer :: i, j, k, o
  integer :: nrn

  hmat = 0.d0
  dhmat= 0.d0
  
  R=ANN%RX

  !hmat
  do i = 1, noopc
    oop(i)=igeom((ncoord-noopc)+i)
  end do

  call ANN%output(igeom(1:R))

  ! Loop over blocks
  nrn = 1
  do i = 1, nstates
    ! PIP on diagonal
    do j = 1, abs(hddef(i,i)%nc)
      hmat(i,i) = hmat(i,i) + ANN%y(nrn)
      nrn = nrn + 1
    end do
    do j = i+1, nstates
      ! PIP*OOP off diagonal
      do k = 1, hddef(i,j)%nc
        do o = 1, 4 ! order
          if (hddef(i,j)%terms(o,k) .eq. 0) cycle
          hmat(i,j) = hmat(i,j) + &
            ANN%y(nrn)*oop(hddef(i,j)%terms(o,k))**o
          nrn = nrn + 1
        end do
      end do
      hmat(j,i) = hmat(i,j)
    end do
  end do

  !dhmat
  call ANN%cal_dadx
  ! Loop over blocks
  nrn = 1
  do i = 1, nstates
    ! PIP on diagonal
    do j = 1, abs(hddef(i,i)%nc)
      ! R part of derivative
      dhmat(1:R,i,i) = dhmat(1:R,i,i) + ANN%dydx(:,nrn)
      ! OOP part of derivative
      dhmat(R+1:ncoord,i,i) = 0.0d0
      nrn = nrn + 1
    end do
    do j = i + 1, nstates
      ! PIP*OOP off diagonal
      do k = 1, hddef(i,j)%nc
        do o = 1, 4 ! order
          if (hddef(i,j)%terms(o,k) .eq. 0) cycle
          ! R part of derivative
          dhmat(1:R,i,j) = dhmat(1:R,i,j) + &
            ANN%dydx(:,nrn)*oop(hddef(i,j)%terms(o,k))**o
          ! OOP part of derivative
          dhmat((R + hddef(i,j)%terms(o,k)),i,j) = &
            dhmat((R + hddef(i,j)%terms(o,k)),i,j) + &
            o*ANN%y(nrn)*oop(hddef(i,j)%terms(o,k))**(o - 1)
          nrn = nrn + 1
        end do
      end do
      dhmat(:,j,i) = dhmat(:,i,j)
    end do
  end do
  
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
  real*8 :: oop,de,fnorm(nstates)
  integer :: i,j,info,lwork
  character(4) :: str, str2
  integer :: ptid, nc
  double precision :: mind
  double precision, external :: dnrm2
  
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
  if(info.ne.0) then
    print *, 'Failed to solve eigenvectors!'
    print *, 'info = ', info
    stop
  end if

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

  ! trajectory parsing
  if (parsing) then
     NEval = NEval + 1
     write(str, "(I4)") natoms*3
     write(str2,"(I4)") nstates
     if (GUNIT<0) call openTrajFile(0)
     if (calcmind) then
        call getmind(igeom,mind,ptid)
        NeighborID = ptid
        print *, "mindNN = ", mind
        mind = mind/sqrt(dble(ndcoord))
        print *, "mindNN'= ", mind
     endif
     fnorm = 0d0
     if (isurftraj <= nstates .and. isurftraj > 0) then
        do j = 1, nstates
           if (j == isurftraj) cycle
           fnorm(j) = dnrm2(3*natoms, cg(1,j,isurftraj),1)
        end do
     end if
     write(GUNIT,&
          "(F16.3,',',I9,',',"//trim(str)//"(F12.7,','),"//trim(str2)//"(E16.8,','),F14.8,',',I6,"//&
          trim(str2)//"(',',F16.4),',',I3,"//trim(str2)//"(',',E16.8),',',I4)")&
          timetraj,NEval,geom, e,mind, ptid, e, isurftraj, fnorm, nc
  endif
  return
end subroutine NNEvaluate
!==================================================================================
subroutine NNEvaluate77(mnatm, mnst, gm, e, cg, h, dcg)
  use nnpes, only: natoms, nstates
  implicit none
  integer, intent(in) :: mnst, mnatm
  double precision, intent(in) :: gm(3*mnatm)
  double precision, intent(out):: e(mnst), h(mnst,mnst)
  double precision, intent(out):: cg(3*mnatm,mnst,mnst), dcg(3*mnatm,mnst,mnst)
  ! local packed variables
  double precision :: gmpck(3*natoms), hpck(nstates,nstates)
  double precision ::  cgpck(3*natoms,nstates,nstates)
  double precision :: dcgpck(3*natoms,nstates,nstates)
  if(natoms>mnatm.or.nstates>mnst) then
    print *, "NNEvaluate77: Maximum on atom or state count insufficient."
    print *, "natoms, mnatm = ", natoms,  ", ", mnatm
    print *, "nstates, mnst = ", nstates, ", ", mnst
    stop
  end if
  gmpck = gm(1:3*natoms)
  call NNEvaluate(gmpck, e, cgpck, hpck, dcgpck)
  h(1:nstates,1:nstates) = hpck
  cg(1:3*natoms,1:nstates,1:nstates) = cgpck
  dcg(1:3*natoms,1:nstates,1:nstates)= dcgpck
  return
end subroutine NNEvaluate77
!==================================================================================
! This subroutine is used to set time and surface index data for output 
SUBROUTINE setTrajData(ttime, tisurf)                                   
  use nnpes
  IMPLICIT NONE                                                         
  DOUBLE PRECISION,intent(IN)  :: ttime                                 
  INTEGER, intent(IN)          :: tisurf                                
                                                                        
  timetraj  = ttime                                                     
  isurftraj = tisurf                                                    
END SUBROUTINE setTrajData                                              
!===============================================================================
subroutine getmind(currgeom,mind,ptid)
  use nnpes
  implicit none
  double precision, dimension(ncoord), intent(in) :: currgeom
  double precision, intent(out) :: mind
  integer, intent(out)          :: ptid

  double precision :: dlpt           ! distance to last point
  double precision :: max_mind       ! max(mind)
  integer, dimension(num_rpts) :: wl ! waiting list of candidate points
  integer :: nwl                     ! number of candidates
  double precision :: d ! distance to point with least lower bound
  double precision, parameter :: vsmall = 1d-8 ! resolution of distances
  integer :: i, j

  ! If first point
  if (firstpt) then
     do i = 1, num_rpts
        call getdist(currgeom, refpoints(:,i), d)
        ldbounds(i) = d
        udbounds(i) = d
        if (i .eq. 1 .or. d < mind) then
           mind = d
           ptid = i
        end if
     end do
     lastrgeom = currgeom
     firstpt = .false.
     return
  end if
  
  ! Get distance to last point
  call getdist(currgeom, lastrgeom, dlpt)
  lastrgeom = currgeom
  do i = 1, num_rpts
     if (dlpt > udbounds(i)) then
        ldbounds(i) = dlpt - udbounds(i)
     else if (dlpt < ldbounds(i)) then
        ldbounds(i) = ldbounds(i) - dlpt
     else
        ldbounds(i) = dble(0)
     end if
  end do
  udbounds = udbounds + dlpt

  ! prepare candidate list
  max_mind = minval(udbounds)
  nwl = 0
  do i = 1, num_rpts
     if (ldbounds(i) < max_mind + vsmall) then
        nwl = nwl + 1
        wl(nwl) = i
     end if
  end do

  ! Evaluate these points to find lowest
  mind = max_mind
  ptid = 0
  do i = 1, nwl
     if (ldbounds(wl(i)) > mind) cycle
     call getdist(currgeom, refpoints(:,wl(i)), d)
     ldbounds(wl(i)) = d
     udbounds(wl(i)) = d
     if (d < mind .or. ptid .eq. 0) then
        mind = d
        ptid = wl(i)
     end if
  end do
  return
end subroutine getmind
!============================================================================
subroutine getdist(geom1, geom2, mind)
  use nnpes
  implicit none
  double precision, dimension(ncoord), intent(in) :: geom1, geom2
  double precision,                    intent(out):: mind
  double precision, dimension(ndcoord) :: geomd
  double precision, external :: dnrm2
  geomd = geom1(dcoordls(1:ndcoord)) - geom2(dcoordls(1:ndcoord))
  mind  = dnrm2(ndcoord, geomd, 1)
  return
end subroutine getdist

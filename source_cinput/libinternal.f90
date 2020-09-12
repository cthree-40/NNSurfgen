!=============================================================================
subroutine buildWBmat(natom,ncoords,cgeom,igeom,Bmat)
  use pip
  use hddata
  implicit none
  integer, intent(in) :: natom, ncoords
  real*8, intent(in) :: cgeom(3*natom)
  real*8, intent(out) :: igeom(ncoords)
  real*8, intent(out) :: Bmat(ncoords,3*natoms)
  integer :: i, cc

  real*8 :: dcoord(12)
  real*8 :: ccoord
  
  
  !PIP
  call evaluatepip(cgeom,igeom(1:npip),Bmat(1:npip,:))

  ! User defined coordinates from coord.input file.
  do i = 1, noopc
     call evaluateoopcoord(i, natom, cgeom, ccoord, dcoord)
     igeom((ncoord-noopc) + i) = ccoord
     Bmat((ncoord-noopc) + i,:) = dcoord
  end do

  return
end subroutine buildWBmat
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
  use hddata
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
        
        val = val + v*coord(n)%acoef(i)

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
!=============================================================================
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
!===============================================================================
! FCDP: compute four-center dot product:
!        f = (a1-a2).(a3-a4) / norm(a1-a2) / norm(a3-a4)
subroutine FCDP(a1,a2,a3,a4,coord,dcoord,id)
  implicit none
  real*8, dimension(3), intent(in) :: a1, a2, a3, a4
  integer,              intent(in) :: id
  real*8, dimension(12),intent(out):: dcoord
  real*8,               intent(out):: coord
  
  real*8, dimension(3) :: d12, d34 ! distance vectors: 1-2, 3-4                   
  real*8               :: r12, r34 ! distances:  1-2, 3-4                         
  real*8               :: f ! dot product                                         
  real*8, dimension(12):: dp ! gradients of f             

  dp = 0d0
  dcoord = 0d0
  coord = 0d0
  
  d12 = a1-a2
  d34 = a3-a4
  r12 = sqrt(dot_product(d12,d12))
  if (abs(r12)<1D-30) r12=sign(1D-30,r12)
  r34 = sqrt(dot_product(d34,d34))
  if (abs(r34)<1D-30) r34=sign(1D-30,r34)
  f = dot_product(d12,d34)
  coord = f / r12 / r34
  
  !dp = d12(1)*d34(1)+d12(2)*d34(2)+d12(3)*d34(3)                                 
  if (id .ne. 0) then
          dp(1)  = d34(1)/r12/r34 - f/r12*d12(1)/r12/r34/r12
          dp(2)  = d34(2)/r12/r34 - f/r12*d12(2)/r12/r34/r12
          dp(3)  = d34(3)/r12/r34 - f/r12*d12(3)/r12/r34/r12
          dp(4)  =-d34(1)/r12/r34 + f/r12*d12(1)/r12/r34/r12
          dp(5)  =-d34(2)/r12/r34 + f/r12*d12(2)/r12/r34/r12
          dp(6)  =-d34(3)/r12/r34 + f/r12*d12(3)/r12/r34/r12
          dp(7)  = d12(1)/r34/r12 - f/r34*d34(1)/r34/r12/r34
          dp(8)  = d12(2)/r34/r12 - f/r34*d34(2)/r34/r12/r34
          dp(9)  = d12(3)/r34/r12 - f/r34*d34(3)/r34/r12/r34
          dp(10) =-d12(1)/r34/r12 + f/r34*d34(1)/r34/r12/r34
          dp(11) =-d12(2)/r34/r12 + f/r34*d34(2)/r34/r12/r34
          dp(12) =-d12(3)/r34/r12 + f/r34*d34(3)/r34/r12/r34
  end if
  dcoord = dp
  
  return
end subroutine FCDP

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

!********************************************************************************************
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
  use progdata, only: natoms
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

!************************************************************************************
! SUBROUTINE oop2 
! OOP angle values and derivatives, with the scaling method of distance coordinates
! f= g*Product[ Scale[rij] ]
! where g is the triple product of displacement vectors
!************************************************************************************
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
SUBROUTINE oop2(natoms,atms,cgeom,qval,bval,scale,sctype,factor)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,sctype
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(natoms*3),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: qval
  DOUBLE PRECISION,INTENT(IN)                     :: scale
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)      :: bval
  double precision,intent(IN)                     :: factor 

  INTEGER           :: i, j, k, sgns(6,4), tp
  DOUBLE PRECISION  :: dnorm2(6), dvecs(6,3),   &
                       geom(4,3), ddvec(3,3),   dw(12),       &
                       coef(2), dwdR(6,12), w(6),  denom, ddenrm
! eprod = product of the scaling exponentials
! fval  = value of unscaled oop coordinate (the triple product)
! fbmat = derivatives of unscaled OOP coordinate
  DOUBLE PRECISION  :: eprod, fval, fbmat(12),efactor

! The oriented direction vectors between all 6 atom pairs
  sgns(1,:)=(/ 1, 0, 0,-1/)
  sgns(2,:)=(/ 0, 1, 0,-1/)
  sgns(3,:)=(/ 0, 0, 1,-1/)
  sgns(4,:)=(/ 1, 0,-1, 0/)
  sgns(5,:)=(/ 0, 1,-1, 0/)
  sgns(6,:)=(/ 1,-1, 0, 0/)

! calculate the scaling factors and their derivatives
  coef(1) = scale
  coef(2) = factor
  tp  = sctype
  if(sctype<0) tp = 0
  k = 1
  dwdR=0d0
  do i=1,3
   do j=i+1,4
     call calcwij(tp,atms(i),atms(j),coef,cgeom,w(k),dw)
     dwdR(k,i*3-2:i*3) = dw(1:3)
     dwdR(k,j*3-2:j*3) = dw(4:6)
     k = k+1
   end do !j
  end do


! extract the geometry of the 4 involved atoms
  do i=1,4
    geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do

! calculate displacement vectors between these atoms
  dvecs=matmul(sgns,geom)
  do i=1,6
    dnorm2(i)=dot_product(dvecs(i,:),dvecs(i,:))
  end do

! calculate value of unscaled OOP angle and its derivatives
  fval=det3(dvecs(1,1),6)
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
    end do
  end do
  CALL cross(dvecs(2,1),6,dvecs(3,1),6,fbmat(1),1)
  CALL cross(dvecs(3,1),6,dvecs(1,1),6,fbmat(4),1)
  CALL cross(dvecs(1,1),6,dvecs(2,1),6,fbmat(7),1)
  fbmat(10:12)=-fbmat(1:3)-fbmat(4:6)-fbmat(7:9)

! calculate the scaled oop
  if(sctype<0)then
    denom = (factor/sum(w))**(scale+3)
    qval = fval * denom
    bval = fbmat* denom - qval*(scale+3)/sum(w)*&
              (dwdR(1,:)+dwdR(2,:)+dwdR(3,:)+dwdR(4,:)+dwdR(5,:)+dwdR(6,:))
  else
    bval = 0d0
    eprod = product(w)
    qval = fval * eprod
    do i=1,6
      efactor = 1D0
      do j=1,6
        if(i.ne.j) then
          efactor=efactor*w(j)
        end if
      end do
      bval = bval + dwdR(i,:)*efactor
    end do
    bval = fbmat*eprod + fval*bval
  end if
CONTAINS

  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m,ldm)
    IMPLICIT NONE
    integer :: ldm
    double precision ::  det3, m(ldm,3)
    det3=(m(1,2)*m(2,3)-m(1,3)*m(2,2))*m(3,1)+&
         (m(1,3)*m(2,1)-m(1,1)*m(2,3))*m(3,2)+&
         (m(1,1)*m(2,2)-m(1,2)*m(2,1))*m(3,3)
  END FUNCTION det3

  !this function calculate the cross product of two 3D vectors
  SUBROUTINE cross(v1,iv1,v2,iv2,v3,iv3)
    implicit none
    integer,intent(in)  ::  iv1,iv2,iv3  !increment of v1,v2 and v3
    double precision,dimension(*) :: v1,v2,v3
    integer             ::  i12,i13,i22,i23,i32,i33
    i12=1+iv1
    i13=1+2*iv1
    i22=1+iv2
    i23=1+2*iv2
    i32=1+iv3
    i33=1+2*iv3

    v3(  1)=-v1(i13)*v2(i22)+v1(i12)*v2(i23)
    v3(i32)= v1(i13)*v2(  1)-v1(  1)*v2(i23)
    v3(i33)=-v1(i12)*v2(  1)+v1(  1)*v2(i22)
  END SUBROUTINE cross

END SUBROUTINE oop2

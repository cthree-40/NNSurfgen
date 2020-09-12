program main
  use nnpes
  use opttools, only: readcolgeom, analysegeom, getflunit, analysegeom2
  implicit none
  ! masses = atomic masses array
  ! anum   = atomic number array
  ! aname  = atomic name   array
  double precision, dimension(:), allocatable :: masses, anum
  character(1),     dimension(:), allocatable :: aname
  ! cgeom  = cartesian geometry
  ! hess   = hessian
  ! w      = frequencies
  double precision, dimension(:),  allocatable :: cgeom
  double precision, dimension(:,:),allocatable :: hess
  double precision, dimension(:),  allocatable :: w
  ! gmfile     = geometry file
  ! cmdstr     = command line argument string
  ! new_geomfl = new geometry file
  ! old_geomfl = old geometry file
  character(300) :: gmfile, cmdstr
  character(300) :: old_geomfl, new_geomfl
  ! ios    = input/output status
  ! isurf1 = surface 1 to optimize index
  ! isurf2 = surface 2 to optimize index
  integer :: ios
  integer :: isurf1, isurf2
  logical :: ex
  ! shift = value of shift
  ! gtol  = energy gradient tolerance for convergence
  ! dtol  = energy difference tolerance for convergence
  double precision :: shift
  double precision :: gtol, dtol
  double precision :: evalcutoff, maxd, scale
  integer          :: maxiter
  ! MAXCONS  = maximum number of constraints
  ! cons_val = value of constraints
  ! constype = 1=stretch, 2=bend, 3=OOP
  integer, parameter :: MAXCONS = 20
  integer            :: ncons
  integer,          dimension(MAXCONS)   :: constype
  integer,          dimension(4,MAXCONS) :: cons_atm
  double precision, dimension(MAXCONS)   :: cons_val

  double precision, parameter :: deg2rad=0.01745329251994d0
  integer,          parameter :: fid=326
  integer :: i
  ! Namelist for input
  NAMELIST /MEXOPT/  shift, gtol, dtol, evalcutoff, maxd, maxiter,&
    scale, ncons, cons_atm, cons_val
  
  new_geomfl="new.geom"
  
  print *," ***************************************** "
  print *," *    findmex.x                           * "
  print *," ***************************************** "
  print *," Minimum energy crossing search on Surfgen potential surfaces"
  
  ! Initialize NN surface
  call pesinit()

  shift=1d-1
  gtol =1d-5
  dtol =1d-5
  maxd =1d-1
  evalcutoff=1d-9
  maxiter = 100
  scale = 1d0
  ncons = 0
  cons_atm=0
  cons_val=0d0
  inquire(file="mexopt.in",exist=ex)
  if(ex)then
    open(unit=fid,file="mexopt.in",access='sequential',form='formatted',&
      IOSTAT=ios, position='rewind',action='read',status='old')
    if(ios.eq.0)then
      read(unit=fid,nml=MEXOPT)
      if(ncons>MAXCONS)ncons=MAXCONS
      if(any(cons_atm(1:2,1:ncons)<1) .or. any(cons_atm(3,1:ncons)<0) .or. &
        any(cons_atm(1,1:ncons)==cons_atm(2,1:ncons)).or.  &
        any(cons_atm(1,1:ncons)==cons_atm(3,1:ncons)).or.  &
        any(cons_atm(1,1:ncons)==cons_atm(4,1:ncons)).or.  &
        any(cons_atm(2,1:ncons)==cons_atm(3,1:ncons)).or.  &
        any(cons_atm(2,1:ncons)==cons_atm(4,1:ncons)).or.  &
        any(cons_atm(:,1:ncons)>natoms)) then
        print *,"input atom index for each coordinate constraints:"
        print "(4I4)",cons_atm(:,1:ncons)
        stop "invalid atom index"
      end if
      close(fid)
    end if
  end if
  
  ! allocate arrays
  allocate(masses(natoms))
  allocate(anum(natoms))
  allocate(aname(natoms))
  allocate(hess(natoms*3,natoms*3))
  allocate(w(3*natoms))
  allocate(cgeom(3*natoms))

  print "(A,I6)","Number of Atoms:  ",natoms
  print "(A,I6)","Number of States: ",nstates

  ! process arguments
  ! synposis:    findmex.x geom state1 state2
  ! Default values:
  ! geom        geom
  ! states      1 2
  isurf1 = 1 
  isurf2 = 2
  call get_command_argument(number=1,value=gmfile,status=ios)
  if(ios<-1)  &     ! input argument larger than container
    stop "Filename truncated.  Please use less than 300 characters"
  if(ios>0) then
    print *,"No filename suppied.  Using default."
    write(gmfile,"(A)"), "geom"
  else
    call get_command_argument(number=2,value=cmdstr,status=ios)
    if(ios<0) stop "Cannot get surface number from command line options"
    if(ios==0)then
      read(cmdstr,*,iostat=ios)isurf1
      if(ios/=0) stop "Cannot get surface number from command line options."
      call get_command_argument(number=3,value=cmdstr,status=ios)
      if(ios/=0) stop "Cannot get surface number from command line options."
      read(cmdstr,*,iostat=ios)isurf2
      if(ios/=0) stop "Cannot get surface number from command line options."
    end if
  end if
  ! check values
  if(isurf1<1.or.isurf1>isurf2.or.isurf2>nstates) then
    stop "Invalid state range. Correct usage: findmex.x gmfile s1 s2"
  end if
  if(isurf1.eq.isurf2) print *, "ATTENTION: Single surface optimization."
  print *,"Reading input from input file "//trim(gmfile)
  call readColGeom(gmfile,1,natoms,aname,anum,cgeom,masses)
 
  ! print initial geometry information
  print *,"--------------- Initial Geometries ------------------"
  call analysegeom(natoms,cgeom,aname,anum,masses,2.5d0,.true.)
  
  ! Set constraint types
  do i=1, ncons
    if (cons_atm(3,i) .eq. 0) then
      ! Stretch
      constype(i) = 1
    else if (cons_atm(4,i) .eq. 0) then
      ! Angle
      constype(i) = 2
    else
      ! OOP
      cons_val(i) = cons_val(i) * deg2rad
      constype(i) = 3
    end if
  enddo
  
  
  ! search for intersections
  call findx(natoms,nstates,cgeom,isurf1,isurf2,maxiter,shift,evalcutoff,gtol,&
    dtol,maxd,scale,ncons,cons_atm,cons_val,constype)
  
  ! print final geometry information
  print *,"---------------  Final Geometries  ------------------"
  call analysegeom(natoms,cgeom,aname,anum,masses,2.5d0,.true.)
  
  ! print final geometry
  call analysegeom2(natoms,cgeom,aname,anum,masses,new_geomfl)
  
  ! deallocate arrays
  deallocate(masses)
  deallocate(anum)
  deallocate(aname)
  deallocate(hess)
  deallocate(w)
  deallocate(cgeom)
  
contains
  !*
  ! bend: bending angle with specific scaling
  !*
  SUBROUTINE BEND(natoms,atms,cgeom,qval,bval,scale,sctype,factor)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                              :: natoms,sctype
    INTEGER,DIMENSION(3),INTENT(IN)                 :: atms
    DOUBLE PRECISION,DIMENSION(3,natoms),INTENT(IN) :: cgeom
    DOUBLE PRECISION,INTENT(OUT)                    :: qval
    DOUBLE PRECISION,INTENT(IN)                     :: scale
    DOUBLE PRECISION,DIMENSION(9),INTENT(OUT)       :: bval
    double precision,intent(IN)                     :: factor
    
    double precision ::  v, dv(9), d1(3),d2(3), d12, d22, dd12(9), dd22(9)
    double precision ::  p, dp(9), tmp, pw
    integer,dimension(3,9) ::   dd1,dd2    ! derivatives of d1 and d2
    integer  ::  i
    
    dd1(1,:) =[1,0,0,0,0,0,-1,0,0]
    dd1(2,:) =[0,1,0,0,0,0,0,-1,0]
    dd1(3,:) =[0,0,1,0,0,0,0,0,-1]
    dd2(1,:) =[0,0,0,1,0,0,-1,0,0]
    dd2(2,:) =[0,0,0,0,1,0,0,-1,0]
    dd2(3,:) =[0,0,0,0,0,1,0,0,-1]
    
    ! calculate the dot product of the two borders of the angle and its 
    ! derivative v is the dot product, and dv is the derivative
    
    d1 = cgeom(:,atms(1))-cgeom(:,atms(3))
    d2 = cgeom(:,atms(2))-cgeom(:,atms(3))
    v = dot_product(d1,d2)
    do i=1,9
      dv(i) = dot_product(d1,dd2(:,i))+dot_product(dd1(:,i),d2)
    end do!i=1,9
    
    ! d12 and d22 are the square norm of the two borders, and dd12 and dd22 are
    ! derivatives of these two quantities.
    d12 = dot_product(d1,d1)
    d22 = dot_product(d2,d2)
    do i=1,9
      dd12(i) = 2*dot_product(d1,dd1(:,i))
      dd22(i) = 2*dot_product(d2,dd2(:,i))
    end do!i=1,9
    
    ! calculate the cosine of angle p = v/sqrt(d12*d22)
    p = v/sqrt(d12*d22)
    dp =  ( dv - v/2*(dd12/d12+dd22/d22)  ) / sqrt(d12*d22)
    
    ! calculate final value and derivatives according to scaling method
    select case(sctype)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !       Cos theta, p
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    case ( 1 )
      qval = p
      bval = dp
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !       Raw angle. scaled by `scale` and shifted by `factor`
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  obtained as acos(p).    d acos(p)= - dp / Sqrt(1-p^2)
    case ( 0 )
      qval = (acos(p)-factor)*scale
      bval = -dp/sqrt(1-p**2)*scale
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !       Scaled cosine, p/( 1+exp[c1*(d1^2+d2^2-c2^2)] )
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !   q = p/( 1+exp[c1*(d1^2+d2^2-c2^2)] )
    case  (2)
      pw=scale*(d12+d22-factor**2)
      if(pw>70d0)then
        qval = 0d0
        bval = 0d0
      else
        tmp  = exp(pw)
        qval = p/(1+tmp)
        bval = (  dp  -  scale*tmp*p/(1+tmp)*(dd12+dd22)  ) /(1+tmp)
      end if
    end select!case(sctype)
  END SUBROUTINE BEND
  
  !*
  ! calchess: calculate averaged hessian at a certain geometry
  !*
  subroutine calcHess(natoms,cgeom,nstate,state1,state2,stepsize,hessian,LDH)
    implicit none
    integer, intent(in)          :: natoms, nstate,state1,state2,LDH
    double precision,intent(in)  :: stepsize
    double precision,intent(in)  :: cgeom(3*natoms)
    double precision,intent(out) :: hessian(LDH,3*natoms)
    double precision   ::  mdif
    
    integer           ::   i,  j
    double precision  ::  dispgeom(3*natoms), dgrd(3*natoms)
    real*8    :: h(nstate,nstate),cg(3*natoms,nstate,nstate)
    real*8    :: dcg(3*natoms,nstate,nstate),e(nstate)

    do i=1,3*natoms
      dispgeom=cgeom
      dispgeom(i)=dispgeom(i) - stepsize
      call NNEvaluate(dispgeom,e,cg,h,dcg)
      dgrd  = 0
      do j=state1,state2
        dgrd  = dgrd-cg(:,j,j)
      end do
      dispgeom=cgeom
      dispgeom(i)=dispgeom(i) + stepsize
      call NNEvaluate(dispgeom,e,cg,h,dcg)
      do j=state1,state2
        dgrd = dgrd+cg(:,j,j)
      end do
      hessian(i,1:3*natoms)= dgrd/2/stepsize/(state2-state1+1)
    end do!o=1,3*natoms
    mdif = 0
    do i=1,3*natoms-1
      do j=i+1,3*natoms
        mdif=max(mdif,abs(hessian(i,j)-hessian(j,i)))
        hessian(i,j)=(hessian(i,j)+hessian(j,i))/2
        hessian(j,i)=hessian(i,j)
      end do
    end do
    if(mdif>1d-5)print *,"maximum hermiticity breaking : ",mdif
  end subroutine calcHess

  !*
  ! dvec: compute distance between two atoms
  !*
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
  
  !*
  ! findx: find minimum energy crossing
  !*
  !----search for minimum energy crossings on adiabatic surfaces 
  !This is done by location critical points of the Lagrangian
  ! L=E+sum(l1.O)+sum(l2.dEd)
  ! O are off diagonal elements and dEd are difference of diagonal
  ! elements of the degeneracy state group in an intersection adapted
  ! partially diagonal representation 
  !---------------------
  !After ab initio data are obtained, first construct an intersection
  ! adapted representation.
  !This is then used to construct the Hessian and gradient of the Lagrangian.
  !Newton-Raphson procedure is used to find the critical point.
  !Energy part of the Hessian is done numerically
  subroutine findx(natoms,nstate,cgeom,surf1,surf2,maxiter,shift,evalcutoff,&
    Etol,Stol,maxd,scale,nc,cons_atm,cons_val,constype)
    implicit none
    integer, intent(in) ::  natoms,surf1,surf2,maxiter,nstate,nc,cons_atm(4,nc)
    integer, intent(in) ::  constype(nc)
    double precision,intent(in) :: shift,Etol,Stol
    double precision,intent(in) :: evalcutoff,maxd,scale,cons_val(nc)
    double precision,dimension(3*natoms),intent(inout) :: cgeom

    integer   :: ndeg  ! number of degenerate states
    integer   :: noffd ! number of off-diagonals that has to vanish 
    integer   :: nlag  ! number of Lagrange multipliers 
    integer   :: neq   ! total number of equations 
    real*8    :: h(nstate,nstate),cg(3*natoms,nstate,nstate)
    real*8    :: dcg(3*natoms,nstate,nstate),e(nstate)
    double precision,dimension(3*natoms)          :: grad
    double precision,dimension(:,:),allocatable   :: hess
    double precision,dimension(:,:),allocatable   :: B  
    
    double precision,dimension(:),allocatable     :: WORK
    double precision,dimension(:),allocatable     :: coord_val
    ! lagrange multipliers and gradient of lagrangian
    double precision,dimension(:),allocatable     :: lam, dlag
    ! evals of hessian and intermediate vectors
    double precision,dimension(:),allocatable     :: w,b1,b2
    integer,dimension(:),allocatable  :: IWORK
    double precision                  :: nrmG, nrmD,tmp(1),bval(12)
    double precision                  :: transrot(3*natoms,6),center(3),ovlp(6)
    double precision, external        :: dnrm2
    double precision, parameter       :: amu2au=1.822888484514D3
    double precision, parameter       :: au2cm1=219474.6305d0
    
    integer           :: offs
    integer           :: LIWORK, LWORK, itmp(1),INFO
    integer           :: iter  , i,j,lindex,nskip
    logical,dimension(:),allocatable :: skip
    logical, parameter  :: debugmode=.false.
    
    allocate(coord_val(nc))
    ! allocate arrays
    ndeg = surf2-surf1+1
    noffd= ndeg*(ndeg-1)/2
    nlag = ndeg-1+noffd+nc
    neq  = 3*natoms+nlag
    allocate(hess(neq,neq))
    allocate(lam(nlag))
    allocate(dlag(neq))
    allocate(skip(neq))
    allocate(w(neq))
    allocate(b1(neq))
    allocate(b2(neq))
    allocate(B(3*natoms,nlag))
    ! initialize work spaces
    call DSYEVD('V','U',neq,hess,neq,w,tmp,-1,itmp,-1,INFO)
    if(info/=0)print *,"DSYEVD allocation failed.  info=",info
    LWORK = int(tmp(1))
    LIWORK= itmp(1)
    allocate(WORK(LWORK))
    allocate(IWORK(LIWORK))
    
    print "(3(A,I4),A)","Searching for minimum between surfaces ",surf1," and ",&
      surf2," in ",maxiter," iterations."
    print "(A)","  Convergence Tolerances"
    print "(A,E10.2,A,E10.2)","  Energy Gradient: ",Etol,"   Displacement:",Stol
    lam = 0d0
    transrot=0d0
    do i=1,3
      do j=0,natoms-1
        transrot(j*3+i,i)=1d0
      end do
      transrot(:,i)=transrot(:,i)/dnrm2(3*natoms,transrot(:,i),1)
    end do
    do iter=1,maxiter
      !recenter geometry
      center =0d0
      do i=1,natoms
        center=center+cgeom(i*3-2:i*3)
      end do
      center=center/natoms
      do i=1,natoms
        cgeom(i*3-2:i*3)=cgeom(i*3-2:i*3)-center
      end do
      !get translational vectors
      transrot(:,4:6)=0d0
      do i=1,natoms
        transrot(i*3-1,4)=-cgeom(i*3)   
        transrot(i*3,4)  = cgeom(i*3-1)   
      end do
      do i=1,natoms
        transrot(i*3-2,5)= cgeom(i*3)   
        transrot(i*3,5)  =-cgeom(i*3-2)   
      end do
      do i=1,natoms
        transrot(i*3-2,6)=-cgeom(i*3-1)   
        transrot(i*3-1,6)= cgeom(i*3-2)   
      end do
      do i=4,6
        transrot(:,i)=transrot(:,i)/dnrm2(3*natoms,transrot(:,i),1)
      end do
      call NNEvaluate(cgeom,e,cg,h,dcg,.false.)
      print "(A,I5)","Iteration ",iter
      print "(A,20F20.4)","Energies : ",e*au2cm1 
      print "(A,20F20.4)","Lagrange Multipliers : ",lam
      ! use h and cg to store partially diagonalized representation and their gradients
      h = 0d0
      do i=1,nstate
        h(i,i) = e(i)
      end do
      ! construct Intersection Adapted Partially Diagonalized Representation
      call OrthogonalizeGH(h,cg,surf1,surf2,3*natoms,100,1d-10)
      ! get Hessian matrix of coordinate block
      call calcHess(natoms,cgeom,nstate,surf1,surf2,1D-6,hess,neq)
      ! get Hessian of mixed-coord-Lagrange multiplers block
      lindex=3*natoms+1
      do i=surf1+1,surf2
        hess(lindex,1:3*natoms)=cg(:,i,i)-cg(:,surf1,surf1) 
        hess(1:3*natoms,lindex)= hess(lindex,1:3*natoms)
        lindex=lindex+1
      end do
      ! part 2 : off diagonal element
      do i=surf1,surf2-1
        do j=i+1,surf2
          hess(lindex,1:3*natoms)=cg(:,i,j)
          hess(1:3*natoms,lindex)=cg(:,i,j)
          lindex=lindex+1
        end do
      end do
      ! part 3 : coordinate constraint part
      do i=1,nc
        hess(:,lindex) =  0d0
        if(constype(i) .eq. 3)then
          ! it is an OOP coordinate
          call oop(natoms,cons_atm(:,i),cgeom,coord_val(i),bval,dble(0));
          do j=1,4
            offs=3*(cons_atm(j,i)-1)
            hess(offs+1:offs+3,lindex) = bval(3*j-2:3*j)
          end do!j=1,4
        else if(constype(i) .eq. 2)then
          ! it is an angle. convert to degrees
          call bend(natoms,cons_atm(:,i),cgeom,coord_val(i),bval,&
            5.72957795131d1,0,0d0)
          do j=1,3
            offs = 3*(cons_atm(j,i)-1)
            hess(offs+1:offs+3,lindex) =  bval(3*j-2:3*j)
          end do!j=1,3
        else if (constype(i) .eq. 1) then
          !it is a distance
          call stre(natoms,cons_atm(1,i),cons_atm(2,i),cgeom,coord_val(i),bval)
          do j=1,2
            offs = 3*(cons_atm(j,i)-1)
            hess(offs+1:offs+3,lindex) =  bval(3*j-2:3*j)
          end do!j=1,2
        end if
        hess(lindex,1:3*natoms)=hess(1:3*natoms,lindex)
        lindex=lindex+1
      end do!i=1,nc
      !the LM-LM block is zero
      hess(3*natoms+1:neq,3*natoms+1:neq)=0d0
      ! get average gradient
      dlag = 0d0
      do i=surf1,surf2
        dlag(1:3*natoms) = dlag(1:3*natoms)+cg(:,i,i)
      end do
      dlag=dlag/(surf2-surf1+1)*scale
      if(debugmode)then
        do i=1,neq
          call DGEMV('T',3*natoms,6,1d0,transrot,3*natoms,hess(1:3*natoms,i),&
            1,0d0,ovlp,1)
          print "(A,I3,A,6F10.4)","Hessian column ",i,&
            ",overlap with translation/rotation: ",ovlp
        end do
        call DGEMV('T',3*natoms,6,1d0,transrot,3*natoms,dlag(1:3*natoms),&
          1,0d0,ovlp,1)
        print "(A,6F10.4)",&
          "Original lagrangian grad overlap with translation/rotation: ",ovlp
      end if
      ! calculate optimized lagrange multipliers and update lagrangian gradient
      B = hess(1:3*natoms,3*natoms+1:neq)
      call updateLag(B,3*natoms,nlag,lam,dlag)
      ! calculate gradient of lagragian: coordinate part
      ! part 1 :  difference in diagonal part
      lindex=1
      do i=surf1+1,surf2
        dlag(1:3*natoms)=dlag(1:3*natoms)+&
          lam(lindex)*(cg(:,i,i)-cg(:,surf1,surf1)) 
        lindex=lindex+1
      end do
      ! part 2 : off diagonal element
      do i=surf1,surf2-1
        do j=i+1,surf2
          dlag(1:3*natoms)=dlag(1:3*natoms)+lam(lindex)*cg(:,i,j)
          lindex=lindex+1
        end do
      end do
      ! part 3: coordinate constraint part
      do i=1,nc
        dlag(1:3*natoms)=dlag(1:3*natoms)+&
          lam(lindex)*hess(lindex+3*natoms,1:3*natoms)
        lindex=lindex+1
      end do
      ! calculate gradient of lagrangian: lagrange multipliers part
      lindex=3*natoms+1
      do i=surf1+1,surf2
        dlag(lindex)=h(i,i)-h(surf1,surf1)
        lindex=lindex+1
      end do
      ! part 2 : off diagonal element
      do i=surf1,surf2-1
        do j=i+1,surf2
          dlag(lindex)=h(i,j)
          lindex=lindex+1
        end do
      end do
      print "(A,15E15.5)",&
        "Coordinate values - constraints :",coord_val-cons_val
      do i=1,nc
        dlag(lindex)=coord_val(i)-cons_val(i)
        lindex=lindex+1
      end do
      if(debugmode)then
        print "(A)","Lagrangian Gradients:"
        print "(3E15.5)",dlag
        call DGEMV('T',3*natoms,6,1d0,transrot,3*natoms,dlag(1:3*natoms),&
          1,0d0,ovlp,1)
        print "(A,6F10.4)",&
          "Lagrangian grad overlap with translation/rotation: ",ovlp
      end if
      do i=1,3*natoms
        hess(i,i)=hess(i,i)+shift
      end do
      ! invert the hessian
      call DSYEVD('V','U',neq,hess,neq,w,WORK,LWORK,IWORK,LIWORK,INFO)
      if(info/=0)print *,"DSYEVD failed.  info=",info
      ! hess. w. hess^T = hess_old
      ! so x= hess_old^-1.b = hess. w^-1. hess^T. b
      ! b1= hess^T.g       =>
      call DGEMV('T',neq,neq,1d0,hess,neq,dlag,1,0d0,b1,1)
      ! b1' = w^-1*b1
      nskip = 0
      do i=1,neq
        if(abs(w(i))<evalcutoff)then
          b1(i)=dble(0)
          nskip=nskip+1
        else!
          b1(i)=b1(i)/w(i)
        end if!
      end do!
      ! b2 = hess.b1'
      call DGEMV('N',neq,neq,1d0,hess,neq,b1,1,0d0,b2,1)
      ! remove rotation and translations
      call DGEMV('T',3*natoms,6,1d0,transrot,3*natoms,b2,1,0d0,ovlp,1)
      print "(A,6F10.6)","Overlap with translation/rotation: ",ovlp
      call DGEMV('N',3*natoms,6,-1d0,transrot,3*natoms,ovlp,1,1d0,b2,1)
      call DGEMV('T',3*natoms,6,1d0,transrot,3*natoms,b2,1,0d0,ovlp,1)
      print "(A,6F10.6)","Overlap with translation/rotation: ",ovlp
      ! cgeom' = cgeom - H^-1.g
      nrmD=dnrm2(3*natoms,b2,1)
      if(nrmD>maxD)then
        b2=b2/nrmD*maxD
        nrmD=maxD
      end if
      nrmG=dnrm2(3*natoms,dlag,1)
      print 1000,nrmG,nrmD
      cgeom = cgeom - b2(1:3*natoms)
      lam   = lam   - b2(3*natoms+1:neq)
      if(nrmG<Etol.and.nrmD<Stol)then
        print *,"Optimization converged"
        return
      end if
    end do
1000 format("         |Grd|=",E12.5,", |Disp|=",E12.5)
  end subroutine findx

  !*
  ! oop: compute OOP coordinate
  !*
  SUBROUTINE oop(natoms,atms,cgeom,qval,bval,scale)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                              :: natoms
    INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
    DOUBLE PRECISION,DIMENSION(natoms*3),INTENT(IN) :: cgeom
    DOUBLE PRECISION,INTENT(IN)                     :: scale
    DOUBLE PRECISION,INTENT(OUT)                    :: qval
    DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)      :: bval
    
    INTEGER           :: i, j, sgns(6,4)
    DOUBLE PRECISION  :: denom, dnorm2(6), dvecs(3,6),&
      geom(4,3), ddvec(3,3)
    DOUBLE PRECISION  :: factor  ! normalization
    
    ! The oriented direction vectors between all 6 atom pairs
    sgns(1,:)=(/ 1, 0, 0,-1/)
    sgns(2,:)=(/ 0, 1, 0,-1/)
    sgns(3,:)=(/ 0, 0, 1,-1/)
    sgns(4,:)=(/ 1, 0,-1, 0/)
    sgns(5,:)=(/ 0, 1,-1, 0/)
    sgns(6,:)=(/ 1,-1, 0, 0/)
    
    ! extract the geometry of the 4 involved atoms
    do i=1,4
      geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
    end do
    
    ! calculate displacement vectors between these atoms
    dvecs=transpose(matmul(sgns,geom))
    do i=1,6
      dnorm2(i)=dot_product(dvecs(:,i),dvecs(:,i))
    end do
    
    ! calculate value of scaled OOP angle and its derivatives
    denom=exp(log(product(dnorm2))/2*scale)
    if(abs(denom)<1d-30)stop "OOP angle cannot be defined when atoms coincide."
    qval=det3(dvecs)/denom
    do i=1,3
      do j=1,3
        ddvec(i,j)=sum(sgns(:,i)*dvecs(j,:)/dnorm2(:))
      end do
    end do
    bval(1:3)=cross(dvecs(:,2),dvecs(:,3))/denom-qval*scale*ddvec(1,:)
    bval(4:6)=cross(dvecs(:,3),dvecs(:,1))/denom-qval*scale*ddvec(2,:)
    bval(7:9)=cross(dvecs(:,1),dvecs(:,2))/denom-qval*scale*ddvec(3,:)
    bval(10:12)=-bval(1:3)-bval(4:6)-bval(7:9)
    
    ! calculate the normalization factor
    factor = 2D0 ** (scale*3-3)*3D0**(scale*3)
    qval = qval * factor
    bval = bval * factor

  end SUBROUTINE oop

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
  !******************************************************************

  !*
  ! orthogonalizegh: 
  !
  ! Rotate degenerate ab initio adiabatic data so that g and h vectors
  ! between any pairs of two degenerate states are orthogonalized.
  ! The rotation is made so that total norm of off diagonal gradients
  ! (h vectors) is minimized.
  ! This subroutine is used as interface for external programs. 
  ! dh             : gradients in the current representation
  ! h              : hamiltonian in the current representation
  ! lstate,ustate  : range of degenerate states
  ! nvibs          : vibrational degrees of freedom
  !
  ! maxiter: maximum number of jacobi rotations
  ! toler  : convergence tolerance for rotation angle beta
  ! hasGrad: specifies if gradient/couple data is available for a block
  ! [method]
  ! This subroutine tries orthogonalize all g-h vectors by an algorithm
  ! similar to the Jacobi's method for eigenvalues.  The Jacobi
  ! transformations here extremize the norm of the corresponding coupling
  ! block rather than eliminating them.   
  SUBROUTINE OrthogonalizeGH(h,dh,lstate,ustate,nvibs,maxiter,toler)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nvibs,nstates,nstates) :: dh
    DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nstates,nstates)       :: h
    INTEGER,INTENT(IN)          :: maxiter,lstate,ustate,nvibs
    DOUBLE PRECISION,INTENT(IN) :: toler
    
    integer           ::  i,j,iter
    integer           ::  mi,mj  !location of maximum rotation
    double precision, dimension(nvibs,nstates,nstates)    :: newgrad
    double precision, dimension(nstates,nstates)          :: newener
    ! required rotation 
    double precision, dimension(nstates,nstates)          :: beta
    double precision           :: max_b,t, c,s
    
    beta=dble(0)  
    max_b=-1
    iter=0
    ! build normalized g and h vectors from Hd and compute rotation
    ! angles the first time
    do i=lstate,ustate 
      do j=lstate,i-1
        beta(i,j) = getBeta(i,j,dh,nvibs,nstates)
        if(abs(beta(i,j))>max_b)then
          max_b=abs(beta(i,j))
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_b)
      end do!j=lstate,i-1
    end do!i=lstate,ustate
    do while(iter<maxiter.and.max_b>toler)
      iter=iter+1
      t=beta(mi,mj)
      c = cos(t)
      s = sin(t)
      ! Gnew = J^T.G.J.   Gnew_ij=Jki*Gkl*Jlj
      newgrad = dh 
      newener = h
      do i=1,nstates
        newgrad(:,i,mi)=dh(:,i,mi)*c-dh(:,i,mj)*s
        newgrad(:,i,mj)=dh(:,i,mj)*c+dh(:,i,mi)*s
        newener(i,mi)=h(i,mi)*c-h(i,mj)*s
        newener(i,mj)=h(i,mj)*c+h(i,mi)*s
      end do
      dh     = newgrad
      h      = newener
      do i=1,nstates
        dh(:,mi,i)=newgrad(:,mi,i)*c-newgrad(:,mj,i)*s
        dh(:,mj,i)=newgrad(:,mj,i)*c+newgrad(:,mi,i)*s 
        h(mi,i)=newener(mi,i)*c-newener(mj,i)*s
        h(mj,i)=newener(mj,i)*c+newener(mi,i)*s 
      end do
      ! update rotation angles
      do i=mj+1,ustate
        beta(i,mj)=getBeta(i,mj,dh,nvibs,nstates)
      end do
      do j=lstate,mj-1
        beta(mi,j)=getBeta(mi,j,dh,nvibs,nstates)
      end do
      do j=mj+1,mi-1
        beta(mi,j)=getBeta(mi,j,dh,nvibs,nstates)
      end do
      max_b=-1
      do i=lstate,ustate
        do j=lstate,i-1
          if(abs(beta(i,j))>max_b)then
            max_b=abs(beta(i,j))
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_b)
        end do!j=lstate,i-1
      end do!i=lstate,ustate 
    end do!while(iter<maxiter.and.max_b>toler)do
  end SUBROUTINE OrthogonalizeGH
  
  FUNCTION getBeta(i,j,dh,nvibs,nst) RESULT(beta)
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: i,j
    DOUBLE PRECISION           :: beta
    integer, intent(in)        :: nst, nvibs
    double precision, dimension(nvibs,nst,nst), intent(in) :: dh
    double precision, dimension(nvibs) :: g,h
    
    g=(dh(:,i,i)-dh(:,j,j))/2
    h=dh(:,i,j)
    CALL orthgh(nvibs,g,h,beta)
    
  END FUNCTION getBeta


  SUBROUTINE orthgh(n,g,h,beta)
    IMPLICIT NONE
    INTEGER,intent(IN)          :: n
    DOUBLE PRECISION,intent(IN) :: g(n),h(n)
    DOUBLE PRECISION,intent(OUT):: beta
    double precision :: denom,num
    num=dot_product(g,h)*2
    denom=-dot_product(g,g)+dot_product(h,h)
    if(abs(denom)<1D-10)then
      PRINT *,"DENORM VERY SMALL.  NO ROTATION IS PERFORMED"
      beta=0D0
    else
      beta=atan(num/denom)/4
    end if
  END SUBROUTINE orthgh
  
  !*************************************************************
  
  !*
  ! stre: calculate the bond length between two reference atoms
  !*
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
  
  !*
  ! updateLag: find optimized lagrange multipliers
  !*
  subroutine updateLag(b,nvibs,nlag,lam,grad)
    implicit none
    integer,intent(in)            ::  nvibs,nlag
    double precision,intent(in)   ::  b(nvibs,nlag)
    double precision,intent(inout)::  lam(nlag)
    double precision,intent(in)   ::  grad(nvibs)
    double precision :: btb(nlag,nlag),v(nlag),dl(nlag),w(nlag)
    integer::LWORK 
    double precision :: WORK(3*nlag*nlag+20)
    integer :: INFO,ipiv(nlag),i
    if(nlag.eq.0)return
    LWORK=3*nlag*nlag+20
    !Bt.(g+B.deltaLam)=0
    !BtB.deltaLam=-Bt.g
    btb = matmul(transpose(b),b)
    
    ! v= -Bt.g
    call dgemv('T',nvibs,nlag,-1d0,B,nvibs,grad,1,0d0,v,1)
    
    call dsysv('U',nlag,1,BtB,nlag,ipiv,v,nlag,WORK,LWORK,INFO)
    if(INFO/=0)then
      print *,"INFO=",INFO
      print *,"Bmatrix:"
      do i=1,nvibs
        print "(30F15.7)",b(i,:)
      end do
      STOP "updateLag: DSYSV Failed."
    end if
    lam = v 
  end subroutine updateLag
  
end program main

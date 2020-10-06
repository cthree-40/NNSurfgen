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
  ! dchess = diagonal correction for hessian
  ! w      = frequencies
  ! skip   = skip 
  double precision, dimension(:),  allocatable :: cgeom
  double precision, dimension(:,:),allocatable :: hess
  double precision, dimension(:),  allocatable :: dchess
  double precision, dimension(:),  allocatable :: w
  logical,          dimension(:),  allocatable :: skip
  ! gmfile = geometry file
  ! cmdstr = command line argument string
  character(255) :: gmfile, cmdstr
  ! ios   = input/output status
  ! isurf = surface to optimize index
  integer        :: ios
  integer        :: isurf
  ! niter       = number of max iterations
  ! egrad_tol   = energy gradient tolerance for convergence
  ! shift       = value of shift
  ! disp_tol    = step size tolerance for convergence
  ! grad_scale  = Scaling of gradient
  ! hess_disp   = displacement for calculation of hessian
  ! maxdisp     = size of maximum optimization displacement
  ! sadd_search = Saddle point search? (Y/N)
  ! printlvl    = print level
  ! molden      = molden output
  ! ant_output  = ANT output
  ! geomprint   = print geometries of intermediate points
  ! bndcutoff   = cutoff for printing bond distances
  integer          :: niter
  double precision :: egrad_tol, disp_tol
  double precision :: shift, grad_scale, hess_disp, maxdisp
  integer          :: printlvl
  logical          :: molden, ant_output, geomprint
  character*1      :: sadd_search
  double precision :: bndcutoff
  integer          :: converge_test
  ! new_geomfl = new geometry file
  ! old_geomfl = old geometry file
  character*300    :: new_geomfl, old_geomfl
  integer :: i
  ! Namelist input for inputfile
  namelist /cpsearch/   niter, egrad_tol, shift, disp_tol, grad_scale, hess_disp, &
      maxdisp, sadd_search, printlvl, molden, ant_output, geomprint
  ! If not read in, here are the default values:
  niter      = 100
  egrad_tol  = 1d-9
  shift      = 1d-5    
  disp_tol   = 1d-5 
  grad_scale = 1.0
  hess_disp  = 1d-5
  maxdisp    = 1d-1 
  sadd_search= 'N'
  old_geomfl = "old.geom"
  new_geomfl = "new.geom"
  printlvl   = 0           
  molden     = .false.       
  geomprint  = .false.    
  ant_output = .true.   
  bndcutoff  = 4d0        

  
  print *, "*******************************************"
  print *, "*         findcp.x                        *"
  print *, "*******************************************"
  print *, " Critical point search on nn-surfgen surface."

  ! Initialize NN surface
  call pesinit()

  ! Allocate arrays
  allocate(masses(natoms))
  allocate(anum(natoms))
  allocate(aname(natoms))
  allocate(hess(3*natoms,3*natoms))
  allocate(w(3*natoms))
  allocate(cgeom(3*natoms))
  allocate(skip(natoms))
  allocate(dchess(3*natoms))
  dchess = 0d0
  skip = .false.

  ! Get command line arguments
  call get_command_argument(number=1,value=gmfile,status=ios)
  if (ios .ne. 0) then
    print *, "Usage: findcp.x [geom file] [state]"
    stop
  end if
  print "('Geometry file: ',A)", trim(adjustl(gmfile))
  call get_command_argument(number=2,value=cmdstr,status=ios)
  if (ios .ne. 0) then
    print *, "Usage: findcp.x [geom file] [state]"
    stop
  end if
  read(cmdstr,*,iostat=ios) isurf
  if (ios .ne. 0) then
    print *, "WARNING: Can't get surface number from input. Using default."
    isurf = 1
  end if

  ! Read geometry file
  print *, "Reading input from geometry file "//trim(adjustl(gmfile))
  call readcolgeom(gmfile, int(1), natoms, aname, anum, cgeom, masses)
  print "(/,A)", "---------- Initial Geometry ----------"
  call analysegeom(natoms,cgeom,aname,anum,masses,bndcutoff,.true.)
  print "(/,A)"," Printing original geometry "
  call analysegeom2(natoms,cgeom,aname,anum,masses,old_geomfl)

  print "(/,A)", "---------- Geometry Optimization ----------"
  call findmin(natoms, nstates, cgeom, isurf, niter, shift, egrad_tol, &
    disp_tol, grad_scale, hess_disp, maxdisp, masses, sadd_search,     &
    converge_test, molden, aname, geomprint, anum)

  print "(/,A)", "---------- Final Geometry ----------"
  call analysegeom(natoms,cgeom,aname,anum,masses,bndcutoff,.true.)
  print "(/,A)", "---------- Harmonic Frequencies ----------"
  call calchess(natoms,cgeom,nstates,isurf,hess_disp,hess,.true.,skip)
  call writehess(hess,3*natoms)
  call getfreq(natoms,masses,hess,w,cgeom,aname,printlvl,molden)
  do i = 1, 3*natoms
    print "(I5,F12.2)", i, w(i)
  end do
  call compute_zeropt(natoms, w, 3*natoms)

  ! Print final geometry to file
  call analysegeom2(natoms,cgeom,aname,anum,masses,new_geomfl)
  if (ant_output) then
    call print_ant_output(natoms, cgeom, aname, anum, masses)
  end if

  ! Deallocate arrays
  deallocate(masses)
  deallocate(anum)
  deallocate(aname)
  deallocate(hess)
  deallocate(w)
  deallocate(cgeom)
  deallocate(skip)
  
contains
  !*
  ! calchess: calculate hessian at a certain geometry
  !*
  subroutine calcHess(natoms,cgeom,nstate,istate,stepsize,hessian,centerd,skip)
    implicit none
    integer, intent(in)          :: natoms, nstate,istate
    logical, intent(in),optional :: skip(natoms)
    logical, intent(in)          :: centerd 
    double precision,intent(in)  :: stepsize
    double precision,intent(in)  :: cgeom(3*natoms)
    double precision,intent(out) :: hessian(3*natoms,3*natoms)
    double precision   ::  mdif
    
    integer           ::   i,  j
    logical           ::   skipdisp(natoms*3)
    double precision  ::  dispgeom(3*natoms), dgrd(3*natoms),gref(3*natoms)
    real*8            ::  h(nstate,nstate),cg(3*natoms,nstate,nstate)
    real*8            ::  dcg(3*natoms,nstate,nstate),e(nstate)
    skipdisp=.false.
    if(present(skip))then
      do i=1,natoms
        if(skip(i))skipdisp(i*3-2:i*3)=.true.
      end do
    end if
    ! to perform backward difference, gradient at references is needed
    if(.not.centerd)then
      call NNEvaluate(cgeom,e,cg,h,dcg)
      gref = cg(:,istate,istate)
    end if
    hessian = 0d0
    do i=1,3*natoms
      if(skipdisp(i))cycle
      dispgeom=cgeom
      dispgeom(i)=dispgeom(i) - stepsize
      call NNEvaluate(dispgeom,e,cg,h,dcg)
      dgrd  =-cg(:,istate,istate)
      if(centerd)then ! centered difference
        dispgeom=cgeom
        dispgeom(i)=dispgeom(i) + stepsize
        call NNEvaluate(dispgeom,e,cg,h,dcg)
        dgrd = dgrd+cg(:,istate,istate)
        hessian(i,:)= dgrd/2/stepsize
      else !backward difference
        dgrd = dgrd+gref
        hessian(i,:)= dgrd/stepsize
      end if
    end do!o=1,3*natoms
    do i=1,3*natoms
      if(skipdisp(i))hessian(:,i)=0d0
    end do
    mdif = maxval(abs(hessian-transpose(hessian)))
    if(mdif>1d-5)print *,"maximum hermiticity breaking : ",mdif
    hessian = (hessian+transpose(hessian))/2
  end subroutine calcHess
  !*
  ! compute_zeropt: compute and print zero point energy.
  !*
  subroutine compute_zeropt(na, freqs, nvibs)
    implicit none
    integer, intent(in) :: na, nvibs
    double precision, dimension(nvibs), intent(in) :: freqs
    integer :: i
    double precision :: zp
    zp = 0d0
    do i = 1, nvibs
      if (freqs(i) .gt. 0d0) then
        zp = zp + freqs(i)
      end if
    end do
    zp = zp / 2d0
    print "(A,f10.2)", "Zero point energy = ", zp
    return
  end subroutine compute_zeropt
  !*
  ! findmin: locate minimum/saddlepoint on NN surface
  !*
  subroutine findmin(natom, nstate, cgeom, isurf, maxiter, shift, etol, &
    stol, gscale, hess_disp, maxd, masses, sadd_search, converge_test,  &
    molden, aname, gmprint, anum)
    implicit none

    integer,          intent(in) :: natom, nstate, isurf, maxiter
    double precision, intent(in) :: shift, etol, stol, gscale
    double precision, intent(in) :: maxd, hess_disp
    character(len=1), intent(in) :: sadd_search
    logical,          intent(in) :: molden, gmprint
    double precision, dimension(natom),   intent(in) :: anum, masses
    character(len=3), dimension(natom),   intent(in) :: aname
    
    double precision, dimension(3*natom), intent(inout) :: cgeom
    integer,                              intent(inout) :: converge_test

    double precision, dimension(nstate) :: e
    double precision, dimension(nstate,nstate) :: h
    double precision, dimension(3*natom,nstate,nstate) :: cg, dcg
    double precision, dimension(3*natom) :: grad, b1, b2, w
    double precision, dimension(3*natom,3*natom)   :: hess, hinv
    double precision, dimension(:), allocatable :: WORK
    integer,          dimension(:), allocatable :: IWORK
    integer  :: LIWORK, LWORK, itmp(1),INFO  
    integer  :: iter, i, eig_start, j
    double precision            :: nrmG, nrmD, tmp(1)
    double precision, external  :: dnrm2
    double precision, parameter :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
    double precision, parameter :: au2ang = 0.529
    double precision, dimension(:), allocatable  :: rem_modes
    character*3,      dimension(natom) :: ctmp
    character*25  :: gflname, itstr
    integer :: ios
    
    
    ! Derived type for normal mode removal in eigenvalue decomposition procedure.
    type mode_removed
      double precision  :: freq
      integer           :: removed 
    end type mode_removed
    ! use this type for normal mode removal list
    type(mode_removed),dimension(:),allocatable :: rem_mode
    
    ! Set convergence flag
    converge_test=0
    
    ! initialize work spaces
    call DSYEVD('V','U',3*natom,hess,3*natom,w,tmp,-1,itmp,-1,INFO)
    if(info/=0)print *,"DSYEVD allocation failed.  info=",info
    LWORK = int(tmp(1))
    LIWORK= itmp(1)
    allocate(WORK(LWORK))
    allocate(IWORK(LIWORK))
    
    ! Allocate mode matrix
    ! This array will allow users to see what modes are being removed
    ! in the eigenvalue decomposition procedure.
    allocate(rem_modes(3*natom))
    allocate(rem_mode(3*natom))
    
    ! Molden output, open molden.all file.
    if (molden) then
      open(file="molden.all",unit=88,status="unknown",action="write", &
        position="rewind")
      write(88,"(' [Molden Format]')")
      write(88,"(' [GEOMETRIES] XYZ')")
    end if
    
    print "(A,I4,A,I4,A)","Searching for minimum on surface ",isurf,&
      " in ",maxiter," iterations."
    print "(A)","  Convergence Tolerances"
    print "(A,E10.2,A,E10.2)","  Energy Gradient: ",Etol,"   Displacement:",Stol
    do iter=1,maxiter
      
      if (molden) then
        write(88, "(i2,X,A)") natom, "/coord"
        write(88, "(i2,X,A)") (iter-1), "/current iter"
        do i = 1, natom
          write(88, "(3X,A,3F11.6)") aname(i), cgeom(i*3-2)*au2ang, &
            cgeom(i*3-1)*au2ang, cgeom(i*3)*au2ang
        end do
      end if
      if (gmprint) then
        write(itstr,"(I4)") iter
        gflname = "geom."//trim(adjustl(itstr))    
        open(unit=89,file=gflname,status="unknown",action="write",iostat=ios)
        if (ios.ne.0) then
          print *, "WARNING! Could not open intermediate geometry file."
        else
          do j = 1, natom
            write(89,"(A3,F7.2,4F14.8)") aname(j), anum(j), &
              cgeom(j*3-2),cgeom(j*3-1),cgeom(j*3),masses(j)
          end do
          close(89)
        end if
      end if
      
      call NNEvaluate(cgeom,e,cg,h,dcg)
      grad = cg(:,isurf,isurf)
      nrmG=dnrm2(3*natom,grad,1)
      call calcHess(natom,cgeom,nstate,isurf,hess_disp,hess,.true.)
      if(iter.eq.1 .or. iter.eq.maxiter)then
        call writeHess(hess,3*natom)
        call getFreq(natom,masses,hess,w,cgeom,ctmp,0,.false.)
        do i=1,3*natom
          print "(I5,F12.2)",i,w(i)
        end do
        do i=1,nstate
          print "(2x,A,I2,A,F12.8)", " Energy of state ", i,"= ", e(i)
        end do
      end if
      hinv = hess
      ! invert the hessian
      ! Call DSYEVD
      call DSYEVD('V','U',3*natom,hess,3*natom,w,WORK,LWORK,IWORK,LIWORK,INFO)
      ! hess now contains the orthonormal eigenvectors of hess(old)
      ! w contains the eigenvalues from lowest to highest
      if(info/=0)print *,"DSYEVD failed.  info=",info
      ! [Old Hessian] = [hess][w][hess]^T
      ! Thus, 
      !    x = [Old Hessian]^-1[b] = ([hess]w)^-1[hess]^T[b]
      !
      ! [b1]= [hess]^T[g]       =>
      ! Call DGEMV. perform [hess]^T[grad]=b1
      ! First, scale gradient
      grad=grad*gscale
      call DGEMV('T',3*natom,3*natom,1d0,hess,3*natom,grad,1,0d0,b1,1)
      ! b1' = w^-1*b1
      ! Check if saddle point search
      if( sadd_search .EQ. 'Y' ) then
        eig_start=2      ! First eigenvalue should be large and negative, so we skip it
        b1(1)=b1(1)/w(1)
      else
        eig_start=1      ! Otherwise, we continue as normal
      end if
      !
      rem_modes=0
      do i=eig_start,3*natom
        if(abs(w(i))<shift)then
          b1(i)=b1(i)
          rem_mode(i)%freq=w(i)
          rem_mode(i)%removed=1
        else!
          b1(i)=b1(i)/w(i)
        end if!
      end do!
      ! Print out to what modes are being removed
      do i=1, 3*natom
        if ( rem_mode(i)%removed .eq. 1 ) then      ! If mode was removed, print info
          write(*,1001) i, rem_mode(i)%freq
        end if
      end do
      ! b2 = hess.b1'
      call DGEMV('N',3*natom,3*natom,1d0,hess,3*natom,b1,1,0d0,b2,1)
      ! cgeom' = cgeom - H^-1.g
      nrmD=dnrm2(3*natom,b2,1)
      if(nrmD>maxD)then
        b2=b2/nrmD*maxD
        nrmD=maxD
      end if
      print 1000,iter,e(isurf)*au2cm1,nrmG,nrmD
      cgeom = cgeom - b2
      if(nrmG<Etol.and.nrmD<Stol)then
        print *,"Optimization converged"
        print "(A,10F20.4)","Final energy of all states : ",e*au2cm1
        converge_test=1
        return
      end if
    end do
    close(88)
    return
1000 format("Iteration ",I4,": E=",F20.4,", |Grd|=",E12.5,", |Disp|=",E12.5)
1001 format("Modes Removed: ",I4,F14.8)
  end subroutine findmin

  !*
  ! gen_molden_file: generate molden frequency file
  !*
  subroutine gen_molden_file(cg, evec, eval, natoms, anames)
    implicit none
    integer, intent(in) :: natoms
    character*3,dimension(natoms) :: anames
    real*8, dimension(3*natoms), intent(in) :: cg
    real*8, dimension(3*natoms,3*natoms), intent(in) :: evec
    real*8, dimension(3*natoms), intent(in) :: eval
    real*8 :: vnorm ! vector norm
    real*8, external :: dnrm2
    
    integer :: mu
    character*25 :: mn
    
    integer :: ios, i, j

    mu = getflunit()
    write(mn,"(a)") "molden.freq"
    mn=trim(mn)
    
    open(file=mn,unit=mu,action='write',status='unknown',iostat=ios)
    if (ios .ne. 0) then
      print "(A)", "Could not open molden file. No file generated."
      return
    end if
    
    ! print header
    write(mu,"(1x,A)") "-- > start of molden input"
    write(mu,"(1x,A)") "[Molden Format]"
    ! print frequencies
    write(mu,"(1x,A)") "[FREQ]"
    do i=1,3*natoms
      write(mu,"(f10.2)") eval(i)
    end do
    ! print geometry
    write(mu,"(1x,A)") "[FR-COORD]"
    do i=1, natoms
      write(mu,"(1x,a3,3f13.6)") trim(anames(i)), &
        cg((i-1)*3+1),cg((i-1)*3+2),cg((i-1)*3+3)
    end do
    ! print modes
    write(mu,"(1x,A)") "[FR-NORM-COORD]"
    do i=1, 3*natoms
      write(mu,"(1x,'vibration',i24)") i
      vnorm = dnrm2(3*natoms, evec(1,i), 1)
      do j=1,natoms
        write(mu,"(3f13.5)") evec((j-1)*3+1,i)/vnorm, &
          evec((j-1)*3+2,i)/vnorm, &
          evec((j-1)*3+3,i)/vnorm
      end do
    end do
    write(mu,"(1x,A)") "--> end of molden input"
  end subroutine gen_molden_file
  
  !*
  ! getfreq: calculate harmonic frequencies from hessian matrix
  !*
  subroutine getFreq(natoms,masses,hess,w,cg,anm,pl,mol)
    implicit none
    integer,intent(in)          :: natoms, pl
    double precision,intent(in) :: masses(natoms),hess(3*natoms,3*natoms)
    double precision,intent(in) :: cg(3*natoms)
    logical,intent(in)          :: mol
    character*3, dimension(natoms), intent(in) :: anm
    double precision,intent(out):: w(3*natoms)
    
    double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
    integer  :: i,j
    double precision  :: sqrm,  hmw(3*natoms,3*natoms), tmp(1)
    double precision,dimension(:),allocatable :: WORK
    integer,dimension(:),allocatable :: IWORK
    integer           :: LIWORK, LWORK, itmp(1),INFO
    ! convert hessian into mass weighed coordinates
    hmw = hess/amu2au
    do i=1,natoms
      sqrm = sqrt(masses(i))
      do j=i*3-2,i*3
        hmw(j,:)=hmw(j,:)/sqrm
        hmw(:,j)=hmw(:,j)/sqrm
      end do
    end do
    !calculate eigenvalues of hmw
    call DSYEVD('V','U',3*natoms,hmw,3*natoms,w,tmp,-1,itmp,-1,INFO)
    if(info/=0)print *,"DSYEVD allocation investigation failed.  info=",info
    LWORK = int(tmp(1))
    LIWORK= itmp(1)
    allocate(WORK(LWORK))
    allocate(IWORK(LIWORK))
    
    ! if print level is greater than 0 we want to print the modes.
    if (pl == 0 .and. .not. mol) then
      call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
    else
      call DSYEVD('V','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
    endif
    if(info/=0)print *,"DSYEVD failed.  info=",info
    
    
    
    do i=1,3*natoms
      if(w(i)<0)then
        w(i) = -sqrt(-w(i))*au2cm1
      else 
        w(i) = sqrt(w(i))*au2cm1
      end if
    end do
    
    ! print modes
    if (pl > 0) then
      print *," Modes:"
      do i=1,3*natoms
        print "(2x,'Mode ',i5,5x,'Frequency: ',f12.2)", i, w(i)
        do j=1, natoms
          print "(2x,3f15.8)", &
            hmw((j-1)*3 + 1,i), &
            hmw((j-1)*3 + 2,i), & 
            hmw((j-1)*3 + 3,i)
        end do
      end do
    end if
    
    ! print molden output?
    if (mol) then
      call gen_molden_file(cg, hmw, w, natoms, anm)
    end if
    return
  end subroutine getFreq

  !*
  ! print_ant_output: print geometry in ANT initial geometry format
  !*
  subroutine print_ant_output(na, g, aname, anum, masses)
    implicit none
    integer, intent(in) :: na
    real*8, dimension(3, na),   intent(in) :: g
    real*8, dimension(na),      intent(in) :: anum, masses
    character(3), dimension(na),intent(in) :: aname
    real*8, parameter :: BOHR2ANG = 0.529177
    integer :: i
    print *, ""
    print *, "Final geometry (ANT format):"
    do i = 1, na
      print 100, trim(aname(i)), masses(i), g(1,i) * BOHR2ANG, &
        g(2,i) * BOHR2ANG, g(3,i) * BOHR2ANG
    end do
    return
100 format(a2,4f14.8)
  end subroutine print_ant_output
    
  !*
  ! writehess: write hessian file to disk, columbus format
  !*
  subroutine writehess(hess,nvibs)
    implicit none
    integer, intent(in)         :: nvibs
    double precision,intent(in) :: hess(nvibs,nvibs)
    integer,parameter     :: hessfl= 32175
    integer  :: ios,i
    
    open (unit=hessfl,file="hessian",status='replace',position='rewind',&
      access='sequential', action='write',iostat=ios)
    if(ios/=0)then
      print *,"Cannot open hessian file for write."
      return
    end if
    do i=1,nvibs
      write(hessfl,"(8F13.6)")hess(i,:)
    end do
    close(hessfl)
  end subroutine writeHess
  
end program main

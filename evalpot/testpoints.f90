!testpoints
!------
!Point testing utility for surfgen.
!
! Written by Xiaolei Zhu, May 2013
! Yarkony group, Johns Hopkins Univeristy
!
! For details see README file.
!
!-- Change Log --
! 12-05-2015 - CLM - Added input file routine.
!                    Added printing of s vector.
!                    Added molden output of g/h/s vectors.
! 12-18-2015 - CLM - Molden output is mass-weighted if selected
!                     in input.
!                    Added print level implementation.
! 04-23-2019 - CLM - Added functionality to generate g-h grid of geoms
!----------------
!This program uses surfgen evaluation library for construction and evaluation of
!Hd, as well as geometry input. 
program testpoints
  use nnpes
  implicit none
  integer,parameter  ::  MaxGeoms = 15000
  real*8,  parameter ::  au2ang = 0.529177249
  character(255)     ::  geomfl 
  integer            ::  npts, i,j,k,l,m,ios,ptid

  character(3),dimension(:),allocatable        :: atoms
  double precision,dimension(:),allocatable    :: anums,masses,e,dvec,norms
  double precision,dimension(:,:),allocatable  :: cgeoms,h
  double precision,dimension(:,:,:),allocatable:: cg,dcg
  double precision,external :: dnrm2

  logical :: printm, mweight ! print molden output, print molden as 
                             ! mass-weighted vectors
  logical :: nnout           ! print input for nn initial guess
  logical :: colout          ! print adiabatic data in COL7 formats
  logical :: grid            ! if MEX, generate grid in g and h
  integer :: ngpts           ! number of grid points along one direction
                             ! total points will be ngpts**2
  real*8 :: dsize            ! total displacement size
  integer :: printl          ! print level
  real*8,  dimension(:,:,:),allocatable :: hvec, svec, gvec !h/s/g vectors
  integer, dimension(:),    allocatable :: mex

  integer :: nbond, nangle, noop
  integer, dimension(:,:), allocatable :: bonds, angles, oopang
  integer, dimension(2) :: mexcalc
  double precision :: ezero
  character(255) :: gname
  logical :: csvout

  character(255), dimension(:,:), allocatable :: colgradfl
  integer, dimension(:,:), allocatable        :: colgradun
  integer :: coleun
  character(255) :: colefl
  double precision :: coleshift
  integer :: ii
  
  double precision, dimension(10) :: bondval, angval, oopval
  double precision, dimension(30) :: freqs
  double precision, dimension(4)  :: mexinfo ! ||g||,||h||,sx,sy
  
  print *,"-------------------------------------------------"
  print *,"Entering testpoints, a surfgen point testing utility"
  print *,""
  print *,"  This program is part of the surfgen program"
  print *,"  2013 Yarkony group, Johns Hopkins University"
  print *,"-------------------------------------------------"

  print *,""
  print *,"Checking for input file"
  allocate(mex(2))
  call read_input_file_outnml(printm, mweight, printl, mex, grid, &
          ngpts, dsize, nnout, colout, coleshift)
  allocate(bonds(2,10))
  allocate(angles(3,10))
  allocate(oopang(4,10))
  call read_input_file_csvnml(nbond, bonds, nangle, angles, noop, oopang, &
          mexcalc, ezero, gname)
  if (nbond .ne. 0) csvout = .true.
  if (csvout) print "(a)", " *.csv output to be printed."
  print *,""
  print *,"Initializing potential"

  call pesinit()

  ! allocate arrays
  allocate(atoms(natoms))
  allocate(anums(natoms))
  allocate(masses(natoms))
  allocate(cgeoms(3*natoms,MaxGeoms))
  allocate(dvec(3*natoms))
  allocate(e(nstates))
  allocate(norms(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))
  allocate(hvec(3*natoms,nstates,nstates))
  allocate(svec(3*natoms,nstates,nstates))
  allocate(gvec(3*natoms,nstates,nstates))

  print *," Number of states:",nstates
  print *," Number of natoms:",natoms
 
  call get_command_argument(number=1,value=geomfl,status=ios)
  if(ios/=0)then
    print *,"Cannot get file name from command line.  Using default."
    write(geomfl,"(A)")"geom.all"
  end if
  print "(A,I8)","Processing geometry input. Maximum number of geometries allowed:",MaxGeoms
  print "(A)","  Filename:"//trim(geomfl)
  npts = MaxGeoms
  call readColGeom(geomfl,npts,natoms,atoms,anums,cgeoms,masses)

  print "(A,I8)","  Number of geometries found in file: ",npts

  if (colout) then ! If we are preparing columbus output, ready files
    allocate(colgradfl(nstates,nstates))
    allocate(colgradun(nstates,nstates))
    ! Set unit numbers
    ii = 31
    do i = 1, nstates
      do j = i, nstates
        colgradun(i,j) = ii
        ii = ii + 1
        if (i .ne. j) colgradun(j,i) = colgradun(i,j)
      end do
    end do
    coleun = ii
    ! Set file names
    write(colefl,"(A)") "energy.all"
    do i = 1, nstates
      write(colgradfl(i,i), &
        "('cartgrd.drt1.state',i1,'.all')") i
      print "(A,A)", "Write adiabatic gradients: ", trim(adjustl(colgradfl(i,i)))
      do j = i + 1, nstates
        write(colgradfl(i,j), &
          "('cartgrd_total.drt1.state',i1,'.drt1.state',i1,'.all')") i, j
        colgradfl(j,i) = colgradfl(i,j)
        print "(A,A)","Write couplings to: ",trim(adjustl(colgradfl(i,j)))
      end do
    end do
    ! Open files
    open(file=trim(adjustl(colefl)),unit=coleun,status="unknown",&
      position="rewind",action="write",iostat=ios)
    if (ios .ne. 0) stop "Could not open energy.all to write!"
    do i = 1, nstates
      open(file=trim(adjustl(colgradfl(i,i))),unit=colgradun(i,i),&
        position="rewind",action="write",iostat=ios)
      if (ios .ne. 0) then
        print "(A,i3,i3)", "Could not open gradient file ", i, i
        stop
      end if
      do j = i + 1, nstates
        open(file=trim(adjustl(colgradfl(i,j))),unit=colgradun(i,j),&
          position="rewind",action="write",iostat=ios)
        if (ios .ne. 0) then
          print "(A,i3,i3)", "could not open gradient file ", i, j
          stop
        end if
      end do
    end do
  end if ! if (colout .eq. .true.)

  do i=1,npts
    print *,""
    print *,"Cartesian geometry"
    print "(3F18.10)",cgeoms(:,i)
    print "(2x,A,I8,A)","Hd predictions for point #",i," :"
    call NNEvaluate(cgeoms(1,i),e,cg,h,dcg)
    print *,""
    print *,"Quasi-diabatic Hamiltonian"
    do j=1,nstates
      print "(2x,10F24.15)",h(j,:)
    end do

    if (nnout) then
      write(*,"('>Hd:  ')",advance="no")
      do j = 1, nstates
        do k = j, nstates
          write(*,"(E20.10)",advance="no") h(j,k)
        end do
      end do
      write(*,"('')")
    endif

    print *,"Diabatic energy(cm-1)"
    print "(2x,10F24.15)",(h(j,j)*219474.6305d0,j=1,nstates)
    print *,""
    print *,"Adiabatic energy(a.u.)"
    print "(2x,10F24.15)",e

    ! Columbus-format adiabatic data
    if (colout) then
      do j = 1, nstates
        write(coleun,"(F14.8)",advance="no") e(j)+coleshift
        write(coleun,"(1x)",advance="no")
      end do
      write(coleun,"('')")
    end if
    
    print *,"Adiabatic energy(cm-1)"
    print "(2x,10F24.15)",e*219474.6305d0
    print *,""
    print *,"Norms of g vectors"
    print "(8x,30I17)",(j,j=1,nstates-1)
    do j=2,nstates
      do k=1,j-1
        dvec=cg(:,j,j)-cg(:,k,k)
        gvec(:,j,k)=dvec/2
        norms(k)=dnrm2(3*natoms,dvec,1)/2
        if (mweight) then
          do l=1,natoms
            do m=1,3
              gvec((l-1)*3+m,j,k)=gvec((l-1)*3+m,j,k)/dsqrt(masses(l))
            end do
          end do
        end if
      end do
      print "(I6,30E17.7)",j,norms(1:j-1)
    end do
    print *,""
    print *,"Norms of h vectors"
    print "(8x,30I17)",(j,j=1,nstates-1)
    do j=2,nstates
      do k=1,j-1
        dvec=cg(:,j,k)*(e(j)-e(k))
        hvec(:,j,k)=dvec
        norms(k)=dnrm2(3*natoms,dvec,1)
        if (mweight) then
          do l=1,natoms
            do m=1,3
              hvec((l-1)*3+m,j,k) = hvec((l-1)*3+m,j,k)/dsqrt(masses(l))
              end do
          end do
        end if
      end do
      print "(I6,30E17.7)",j,norms(1:j-1)
    end do
    print *,""
    print *,"Norms of s vectors"
    print "(8x,30I17)",(j,j=1,nstates-1)
    do j=2,nstates
      do k=1,j-1
        dvec=cg(:,j,j)+cg(:,k,k)
        svec(:,j,k)=dvec/2
        norms(k)=dnrm2(3*natoms,dvec,1)/2
        if (mweight) then
          do l=1,natoms
            do m=1,3
              svec((l-1)*3+m,j,k) = svec((l-1)*3+m,j,k)/dsqrt(masses(l))
            end do
          end do
        end if
      end do
      print "(I6,30E17.7)",j,norms(1:j-1)
    end do
    print *,""
    print *,"Cartesian Gradients and Couplings in Adiabatic Representation"
    do j=1,nstates
      do k=1,j
         print "(A,I3,A,I3,A,E13.5)"," block (",j,",",k,"),|F|=",dnrm2(3*natoms,cg(1,j,k),1)
         print "(3F18.10)",cg(:,j,k)

         ! Columbus adiabatic output
         if (colout) then
           if (j .ne. k) then
             write(colgradun(j,k),"(3D15.6)") cg(:,j,k)*abs(e(j)-e(k))
           else
             write(colgradun(j,k),"(3D15.6)") cg(:,j,k)
           end if
         end if ! if (colout .eq. .true.)
         
       end do!k
    end do!j
    print *,""
    print *,"Cartesian Gradients in Diabatic Representation"
    do j=1,nstates
      do k=1,j
         print "(A,I3,A,I3,A)"," block (",j,",",k,")"
         print "(3F18.10)",dcg(:,j,k)
      end do!k
    end do!j
    if (printm .and. i .eq. npts .and. mex(1) .ne. 0) then
            if (mweight) then
                print *, "ALERT: mass-weighted g, h and s being used!"
            end if
            call rot_g_h_vectors(gvec(:,mex(2),mex(1)),hvec(:,mex(2),mex(1)),&
                    natoms)
            call print_molden_output(gvec, hvec, svec, nstates, &
                    natoms, atoms, cgeoms(:,i), mex, printl)    
!            call rot_g_h_vectors(gvec(:,mex(2),mex(1)),hvec(:,mex(2),mex(1)),&
!                natoms) 
            call express_s_in_ghplane(gvec(:,mex(2),mex(1)), &
                    hvec(:,mex(2),mex(1)), svec(:,mex(2),mex(1)), natoms)
            if (grid) then
                    call generate_grid_gh(cgeoms(:,i),gvec(:,mex(2),mex(1)),&
                            hvec(:,mex(2),mex(1)),ngpts,natoms,dsize)
            endif
            print *, "G-vector: "
            print "(3F13.8)", gvec(1:natoms*3,mex(2),mex(1))
            print *, "H-vector: "
            print "(3F13.8)", hvec(1:natoms*3,mex(2),mex(1))
    end if
    if (csvout) then
            call compute_bond_distances(cgeoms(:,i), natoms, nbond, bonds,&
                    bondval)
            call compute_angle_values(cgeoms(:,i), natoms, nangle, angles, &
                    angval)
            call compute_oop_values(cgeoms(:,i), natoms, noop, oopang, oopval)

            if (mexcalc(1) .eq. mexcalc(2) .and. mexcalc(1) .ne. 0) then
                    ! Minimum or saddle point
                    call compute_vibfreq(cgeoms(:,i), natoms, freqs, mexcalc(1),&
                            nstates, masses)
                    call print_minsad_csvout(nbond, bondval, nangle, angval, &
                            noop, oopval, freqs, e, nstates, natoms, gname,  &
                            ezero)
            else
                    ! MEX point
                    call compute_mexinfo(cgeoms(:,i), natoms, &
                            gvec(:,mexcalc(2),mexcalc(1)),    &
                            hvec(:,mexcalc(2),mexcalc(1)),    &
                            svec(:,mexcalc(2),mexcalc(1)), mexinfo)
                    call print_mexinfo_csvout(nbond, bondval, nangle, angval,&
                            noop, oopval, mexinfo, e, nstates, natoms, gname,&
                            ezero)
            end if

    end if
  end do!i

  ! get neighboring point index
  !call getNeighbor(ptid)
  if(ptid/=0) print *,"Index of Closest Data Point : ", ptid

contains
  ! calcHess: calculate hessian at a certain geometry
  subroutine calcHess(natoms,cgeom,nstate,istate,stepsize,hessian,centerd,skip)
    implicit none
    integer, intent(in)          :: natoms, nstate,istate
    logical, intent(in),optional :: skip(natoms)
    logical, intent(in)          :: centerd !whether to do centered difference or only backward difference    
    double precision,intent(in)  :: stepsize
    double precision,intent(in)  :: cgeom(3*natoms)
    double precision,intent(out) :: hessian(3*natoms,3*natoms)
    double precision   ::  mdif
    
    integer   ::   i,  j
    logical   ::   skipdisp(natoms*3)
    double precision  ::  dispgeom(3*natoms), dgrd(3*natoms),gref(3*natoms)
    real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
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
  
  
  ! compute_angle_values: compute the listed angle values for a
  ! geometry.
  ! geom = geometry
  ! natm = number of atoms
  ! na   = number of angles to compute
  ! a    = atoms of angles (center is vertex)
  ! aval = value of angles computed
  subroutine compute_angle_values(geom, natm, na, a, aval)
    implicit none
    integer, intent(in) :: natm, na
    integer, dimension(3,10), intent(in) :: a
    double precision, dimension(3,natm), intent(in) :: geom
    double precision, dimension(10), intent(inout) :: aval
    double precision, external :: dnrm2
    integer, dimension(3) :: aatm
    double precision, dimension(3) :: d1, d2
    integer :: i
    aval = 0d0
    do i = 1, na
            aatm = a(1:3,i)
            d1 = geom(:,aatm(1))-geom(:,aatm(2))
            d1 = d1/dnrm2(3,d1,1)
            d2 = geom(:,aatm(3))-geom(:,aatm(2))
            d2 = d2/dnrm2(3,d2,1)
            aval(i) = 90.0/acos(0d0)*acos(dot_product(d1,d2))
    end do
    return
  end subroutine compute_angle_values
  
  ! compute_bond_distances: compute the listed bond distances for a
  ! geometry.
  ! geom = geometry
  ! natm = number of atoms
  ! nb   = number of bond distances
  ! b    = bonds to measure
  ! bval = value of bonds
  subroutine compute_bond_distances(geom, natm, nb, b, bval)
    implicit none
    double precision, parameter  ::  bohr2ang=0.529177249d0
    integer, intent(in) :: natm, nb
    integer, dimension(2,10), intent(in) :: b
    double precision, dimension(3,natm), intent(in) :: geom
    double precision, dimension(10), intent(inout) :: bval
    double precision, external :: dnrm2
    integer :: i
    bval = 0
    do i = 1, nb
            bval(i) = dnrm2(3, (geom(:,b(1,i))-geom(:,b(2,i))),1)
            bval(i) = bval(i) * bohr2ang
    end do
    return
  end subroutine compute_bond_distances

  ! compute_mexinfo: compute information for a minimum energy crossing
  subroutine compute_mexinfo(geom, natm, gvec, hvec, svec, mexinfo)
    implicit none
    integer, intent(in) :: natm
    double precision, dimension(3,natm), intent(in) :: geom
    double precision, dimension(3,natm), intent(inout) :: gvec, hvec, svec
    double precision, dimension(4), intent(inout) :: mexinfo
    double precision, dimension(3,natm) :: gx, hy
    double precision, external :: ddot, dnrm2
    call rot_g_h_vectors(gvec, hvec, natm)
    gx = gvec / dnrm2(3*natm, gvec, 1)
    hy = hvec / dnrm2(3*natm, hvec, 1)
    mexinfo(1) = dnrm2(3*natm, gvec, 1)
    mexinfo(2) = dnrm2(3*natm, hvec, 1)
    mexinfo(3) = ddot(3*natm,svec,1,gx,1)
    mexinfo(4) = ddot(3*natm,svec,1,hy,1)
    return
  end subroutine compute_mexinfo
  
  ! compute_oop_values: compute selected out-of-plane angle values
  ! for a geometry
  ! geom = geometry
  ! natm = number of atoms
  ! no   = number of OOP angles
  ! o    = OOP angle to measure
  ! oval = value of OOP angle
  subroutine compute_oop_values(geom, natm, no, o, oval)
    implicit none
    integer, intent(in) :: natm, no
    integer, dimension(2,10), intent(in) :: o
    double precision, dimension(3,natm), intent(in) :: geom
    double precision, dimension(10), intent(inout) :: oval
    double precision, external :: dnrm2
    integer :: i
    integer, dimension(4) :: oatm
    double precision, dimension(3) :: d1, d2, d3, p1, p3
    oval = 0d0
    do i = 1, no
            oatm = o(1:4,i)
            d1 = geom(:,oatm(1))-geom(:,oatm(2))
            d1 = d1/dnrm2(3,d1,1)
            d2 = geom(:,oatm(3))-geom(:,oatm(2))
            d2 = d2/dnrm2(3,d2,1)
            d3 = geom(:,oatm(4))-geom(:,oatm(3))
            d3 = d3/dnrm2(3,d3,1)

            p1 = d1-dot_product(d2,d1)*d2
            p1 = p1/dnrm2(3,p1,1)
            p3 = d3-dot_product(d2,d3)*d2
            p3 = p3/dnrm2(3,p3,1)

            oval(i) = 90/acos(0d0)*acos(dot_product(p1,p3))

    end do
    return
  end subroutine compute_oop_values

  ! compute_vibfreqs: compute vibrational frequencies for a minimum or
  ! saddle point.
  subroutine compute_vibfreq(geom, natm, freqs, istate, nstate, masses)
    implicit none
    integer, intent(in) :: natm, istate, nstate
    double precision, dimension(3,natm), intent(in) :: geom
    double precision, dimension(natm),   intent(in) :: masses
    double precision, dimension(30),  intent(inout) :: freqs

    double precision, dimension(3*natm,3*natm)  :: hessian
    double precision, parameter :: stepsize = 1.0d-5
    logical, parameter :: centerd = .true.
    logical, dimension(3*natm) :: skip
    skip = .false.
    
    call calcHess(natm,geom,nstate,istate,stepsize,hessian,centerd,skip)
    call getFreq(natm,masses,hessian,freqs,geom)
    
    return
  end subroutine compute_vibfreq
  
  ! generate_grid_gh: generate grid of points in g and h
  subroutine generate_grid_gh(geom, gv, hv, ndisps, na, dsize)
    implicit none
    integer, intent(in) :: na, ndisps
    real*8 , intent(in) :: dsize
    real*8, dimension(3*na), intent(in) :: geom, gv, hv
    real*8, dimension(3*na) :: x, y, gdisp
    real*8 :: gnrm, hnrm, dsz
    integer :: flun = 31, ios
    character(25) :: flname = "gh_disp.dat"
    integer :: i, j
    real*8, dimension(nstates) :: e
    real*8, dimension(nstates,nstates) :: h
    real*8, dimension(3*na,nstates,nstates) :: cg
    real*8, dimension(3*na,nstates,nstates) :: dcg
    real*8, external :: ddot
    
    gnrm = ddot(na*3,gv,1,gv,1)
    hnrm = ddot(na*3,hv,1,hv,1)
    x = gv / gnrm
    y = hv / hnrm
    dsz = dsize / ndisps
    
    open(unit=flun,file=trim(adjustl(flname)),status="replace", &
            action="write",iostat=ios)
    if (ios .ne. 0) stop "Error! Could not open gh_disp.dat!"

    do i = 0 - ndisps, ndisps
            do j = 0 - ndisps, ndisps
                    gdisp = geom
                    gdisp = gdisp + dsz * i * x
                    gdisp = gdisp + dsz * j * y
                    call NNEvaluate(gdisp,e,cg,h,dcg)
                    write(flun,"(f10.5,',',f10.5,',')",advance="no") &
                            (dsz * i), (dsz * j)
                    do k = 1, nstates
                            write(flun,"(f20.15,',')",advance="no") e(k)
                    end do
                    write(flun,"")
            end do
    end do

    close(unit=flun)
    
    
    
  end subroutine generate_grid_gh

  ! getFreq: get harmonic frequencies from hessian matrix
  subroutine getFreq(natoms,masses,hess,w,cg)
    implicit none
    integer,intent(in)          :: natoms
    double precision,intent(in) :: masses(natoms),hess(3*natoms,3*natoms)
    double precision,intent(in) :: cg(3*natoms)
    
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
    call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
    if(info/=0)print *,"DSYEVD failed.  info=",info
    
    
    
    do i=1,3*natoms
            if(w(i)<0)then
                    w(i) = -sqrt(-w(i))*au2cm1
            else
                    w(i) = sqrt(w(i))*au2cm1
            end if
    end do
    
    return
  end subroutine getFreq

  ! rot_g_h_vectors: rotate g and h vectors
  subroutine rot_g_h_vectors(g, h, na)
    implicit none
    integer,intent(in) :: na
    real*8, dimension(3*na), intent(inout) :: g, h
    real*8, dimension(3*na) :: tmp1, tmp2
    real*8 :: gh, hh, gg, atin, theta
    real*8, external :: ddot
    
    gh = ddot(na*3,g,1,h,1)
    hh = ddot(na*3,h,1,h,1)
    gg = ddot(na*3,g,1,g,1)
    atin = (2 * gh) / (hh - gg)
    theta= datan(atin) / 4.0
    
    tmp1 = g
    tmp2 = h
    print "(A,f13.5)"," g.h = ", gh
    print "(A,f8.3)"," Rotating by ", theta
    g = dcos(2 * theta) * tmp1 - dsin(2 * theta) * tmp2
    h = dsin(2 * theta) * tmp1 + dcos(2 * theta) * tmp2
    gh = ddot(na*3,g,1,h,1)
    print "(A,f13.5)"," g.h = ", gh
    
    return
  end subroutine rot_g_h_vectors
  
  ! express_s_in_ghplane: gives components of s in g/h plane.
  subroutine express_s_in_ghplane(g, h, s, na)
          implicit none
          integer, intent(in) :: na
          real*8, dimension(3*na), intent(in) :: g, h, s
          real*8, dimension(3*na) :: gx, hy
          real*8, external :: ddot
          real*8, dimension(2) :: sxy
          
          gx = g / dsqrt(ddot(3*na,g,1,g,1))
          hy = h / dsqrt(ddot(3*na,h,1,h,1))

          sxy(1) = ddot(3*na,s,1,gx,1)
          sxy(2) = ddot(3*na,s,1,hy,1)
          
          print "('||s|| = ',f15.8)", dnrm2(3*na,s,1)
          print "('S vector in g-h plane:')"
          print "('S(x,y) = (',f15.8,',',f15.8,')')", sxy(1), sxy(2)
        
          return
  end subroutine express_s_in_ghplane

  ! read_input_file_csvnml: reads input file (testpoints.in) csvoutput
  ! namelist for testpoints.x
  ! nbond (bonds) = number of bonds (bonds)
  ! nangle (angles) = number of angles (angles)
  ! noop (oopang) = number of out of plane angles (oop angles)
  ! mexcalc = states that form mex
  ! ezero = shift to energy output
  ! gname = name for geometry
  subroutine read_input_file_csvnml(nbond, bonds, nangle, angles, noop, &
          oopang, mexcalc, ezero, gname)
    implicit none
    integer, intent(out) :: nbond, nangle, noop
    integer, dimension(2,10), intent(inout) :: bonds
    integer, dimension(3,10), intent(inout) :: angles
    integer, dimension(4,10), intent(inout) :: oopang
    integer, dimension(2),    intent(inout) :: mexcalc
    double precision, intent(out) :: ezero
    character(255), intent(out) :: gname
    
    character(25) :: flname = 'testpoints.in'
    integer       :: flunit = 21, ios
    integer, parameter :: fl_not_found = 29
    
    namelist /csvoutput/ nbond, bonds, nangle, angles, noop, oopang, &
            mexcalc, ezero, gname
    ! Defaults
    mexcalc = 0
    ezero   = 0d0
    oopang  = 0
    noop    = 0
    angles  = 0
    nangle  = 0
    bonds   = 0
    nbond   = 0
    gname   = ""

    ! Open input file. If not found return, doing nothing.
    open(file = flname, unit = flunit, status = 'old', &
            action = 'read', iostat = ios)
    if (ios .eq. 0) then
            read(unit = flunit, nml = csvoutput, iostat = ios)
            if (ios .lt. 0) then
                    print "(a,i5)", "No input file found. ios = ", ios
                    return
            end if
    else
            print "(a,i5,a)", " ** Warning: ", ios, &
                    " occurred opening input file! **"
            return
    end if
    print *, "Input file found. Reading csvoutput namelist."
    return
    
  end subroutine read_input_file_csvnml
  
  ! read_input_file: reads input file (testpoints.in) testoutput namelist
  !  for testpoints.x
  subroutine read_input_file_outnml(printm, mweight, printl, mex, grid, &
          ngpts, dsize, nnout, colout, coleshift)
          implicit none
          logical, intent(out) :: printm, mweight, grid, nnout, colout
          integer, intent(out) :: printl, ngpts
          real*8,  intent(out) :: dsize, coleshift
          integer, dimension(2), intent(inout) :: mex

          character(25) :: flname = 'testpoints.in'
          integer       :: flunit = 21, ios
          integer, parameter :: fl_not_found = 29

          namelist /testoutput/ printm, mweight, printl, mex, grid, ngpts, &
                  dsize, nnout, colout, coleshift
          printm = .false.
          mweight= .false.
          grid=.false.
          nnout=.false.
          colout=.false.
          ngpts=0
          dsize = 0.01
          printl = 0
          mex(1) = 0
          mex(2) = 0
          coleshift = 0.0d0

          ! Open input file. If file is not found, set printm = .false.
          ! and return. Else, read the file.
          open(file = flname, unit = flunit, status = 'old', &
                  action = 'read', iostat = ios)
          if (ios .eq. 0) then
                  read(unit = flunit, nml = testoutput)
                  return
          else
                  print "(a,i5,a)", "  ** Warning ", ios, &
                          " occurred opening input file! **"
                  return
          end if
          print *, "Input file found."    
          return
  end subroutine read_input_file_outnml

  ! read COLUMBUS geom file and obtain geometry and atom info
  SUBROUTINE readColGeom(gfile,ngeoms,na,atoms,anums,cgeom,masses)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)                           :: gfile
    INTEGER,INTENT(IN)                                    :: na
    INTEGER,INTENT(INOUT)                                 :: ngeoms
    CHARACTER(3),dimension(na),INTENT(INOUT)              :: atoms
    DOUBLE PRECISION,dimension(na),INTENT(INOUT)          :: anums,masses
    DOUBLE PRECISION,dimension(3*na,ngeoms),INTENT(INOUT) :: cgeom
    INTEGER                                               :: i,j,k,GUNIT,ios
    
    GUNIT=get_unit()
    print "(A,I4,2A)","unit of geom file = ",GUNIT,", filename=",trim(adjustl(gfile))
    open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted',&
        status='old',action='read',position='rewind',iostat=ios)
    if(ios/=0)then
        print *,"Failed to open file [", trim(adjustl(gfile)),"]"
        ngeoms = 0
        close(GUNIT)
        return
    end if
    do i = 1,ngeoms
        do j = 1,na
            read(GUNIT,*,iostat=ios)atoms(j),anums(j),(cgeom(3*(j-1)+k,i),k=1,3),masses(j)
            if(ios/=0)then
                print *,"  End of file encountered. i=",i,",j=",j
                ngeoms=i-1
                close(GUNIT)
                return
            end if
        enddo
    enddo
    close(GUNIT)
    
    return
  END SUBROUTINE readColGeom
  !-------------------------------------------------------------------
  ! get_unit: get available file unit
  function get_unit () result(u)
    implicit none
    integer :: i
    integer :: u
    logical :: uexist, uopen
    u = 0
    do i = 15, 9999
            inquire(unit=i,exist=uexist,opened=uopen)
            if (uexist .and. .not. uopen) then
                    u=i
                    exit
            end if
    end do
    if (u .eq. 0) stop "No unit available."
  end function get_unit

  
  ! print csv output for minimum/saddle point 
  ! Input:
  !  nb   = number of bonds
  !  na   = number of angles
  !  no   = number of out of plane angles
  !  nst  = number of states
  !  natm = number of atoms
  !  nm   = file name
  subroutine print_minsad_csvout(nb, bv, na, av, no, ov, w, e, nst, natm,&
          nm, ez)
    implicit none
    integer, intent(in) :: nb, na, no, nst, natm
    double precision, dimension(10),  intent(in) :: bv, av, ov
    double precision, dimension(30),  intent(in) :: w
    double precision, dimension(nst), intent(in) :: e
    double precision, intent(in) :: ez
    character(255), intent(in) :: nm
    character(255), parameter  :: flname = "geomdata.csv"
    integer, parameter :: flunit = 25
    double precision, parameter :: au2cm1=219474.6305d0
    integer :: ios
    integer :: i, f

    open(file="geomdata.csv",unit=flunit,status="unknown",&
            action="write",position="append",iostat=ios)
    if (ios .ne. 0) return
    ! Geom name
    write(unit=flunit,fmt=10,advance="no") trim(adjustl(nm))
    ! Bond values
    do i = 1, nb
            write(unit=flunit,fmt=11,advance="no") bv(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! Angle values
    do i = 1, na
            write(unit=flunit,fmt=12,advance="no") av(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! OOp values
    do i = 1, no
            write(unit=flunit,fmt=12,advance="no") ov(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! Frequencies
    f=1
    do i = 30, 1, -1
            if (abs(w(i)) .gt. 10d0) then
                    write(unit=flunit,fmt=13,advance="no") w(i)
                    f=f+1
            end if
    end do
    do i = f, 3*natm-6
            write(unit=flunit,fmt="('',',')",advance="no")
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! Energies
    do i = 1, nst
            write(unit=flunit,fmt=14,advance="no") (e(i) * au2cm1) - ez
    end do
    close(unit=flunit)
    return
10  format(1x,A,',')  ! Name
11  format(f10.4,',') ! Bonds
12  format(f10.3,',') ! Angles
13  format(f10.1,',') ! Frequencies
14  format(f10.1,',') ! Energies
  end subroutine print_minsad_csvout

  ! print_mexinfo_csvout: print mex information in csv format
  subroutine print_mexinfo_csvout(nb, bv, na, av, no, ov, mi, e, nst, natm,&
          nm, ez)
    implicit none
    integer, intent(in) :: nb, na, no, nst, natm
    double precision, dimension(10), intent(in) :: bv, av, ov
    double precision, dimension(4),  intent(in) :: mi
    double precision, dimension(nst),intent(in) :: e
    double precision, intent(in) :: ez
    character(255), intent(in) :: nm
    character(255), parameter  :: flname = "geomdata.csv"
    integer, parameter :: flunit = 25
    double precision, parameter :: au2cm1=219474.6305d0
    integer :: ios
    integer :: i

    open(file="geomdata.csv",unit=flunit,status="unknown",&
            action="write",position="append",iostat=ios)
    if (ios .ne. 0) return
    ! Geom name
    write(unit=flunit,fmt=10,advance="no") trim(adjustl(nm))
    ! Bond values
    do i = 1, nb
            write(unit=flunit,fmt=11,advance="no") bv(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! Angle values
    do i = 1, na
            write(unit=flunit,fmt=12,advance="no") av(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! OOp values
    do i = 1, no
            write(unit=flunit,fmt=12,advance="no") ov(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! mexinfo
    do i = 1, 4
            write(unit=flunit,fmt=13,advance="no") mi(i)
    end do
    write(unit=flunit,fmt="('',',')",advance="no")
    ! Energies
    do i = 1, nst
            write(unit=flunit,fmt=14,advance="no") (e(i) * au2cm1) - ez
    end do
    close(unit=flunit)
    return
10  format(1x,A,',')  ! Name
11  format(f10.4,',') ! Bonds
12  format(f10.3,',') ! Angles
13  format(f10.6,',') ! g, h, sx, sy
14  format(f10.1,',') ! Energies
  end subroutine print_mexinfo_csvout
  
  ! print_molden_output: print g/h/s vectors in molden format
  subroutine print_molden_output(gv, hv, sv, ns, na, a, g, mex, pl)
          implicit none
          integer, intent(in) :: ns, na, pl
          integer, dimension(2),          intent(in) :: mex
          character(3), dimension(na),    intent(in) :: a
          real*8,  dimension(3*na,ns,ns), intent(in) :: gv, hv, sv
          real*8,  dimension(3*na),       intent(in) :: g
          real*8,  dimension(:),         allocatable :: fq
          character(25) :: flnm = 'vectors.molden'
          integer       :: flun = 22, ios, i, j, nfq
          
          ! Allocate dummy frequencies. There will be n(n+1)/2, were 
          ! n = states. Open file. Print vectors in order: g, h, s. 
          !Lower triangle is printed.
          
          nfq = 3 * ns * (ns + 1) / 2
          allocate(fq(nfq))
          fq = 1d0

          open(file = flnm, unit = flun, status = 'unknown', action = 'write',&
                  iostat = ios)
          if (ios .ne. 0) then
                  print *, "**Error ", ios, " occurred opening vectors.molden!"
                  return
          end if        
          write(flun,"(1x,'[Molden Format]')")

          call print_molden_freq(flun, fq, 3)
          call print_molden_geom(flun, g, a, na)
          call print_molden_vecs(flun, gv, hv, sv, ns, na, nfq, mex, pl)

          close(flun)
          return
  end subroutine print_molden_output

  ! print_molden_vecs: print vectors in molden format. Expects g/h/s
  subroutine print_molden_vecs(u, gv, hv, sv, ns, na, nfq, mex, pl)
          implicit none
          integer, intent(in) :: nfq, na, ns, u, pl
          integer, dimension(2), intent(in) :: mex
          real*8,  dimension(3*na,ns,ns), intent(in) :: gv, hv, sv
          real*8 :: ngv, nhv, nsv
          integer :: i, j, p
          
          j = mex(1)
          i = mex(2)
          ngv = dnrm2(3*na,gv(:,i,j),1)
          nhv = dnrm2(3*na,hv(:,i,j),1)
          nsv = dnrm2(3*na,sv(:,i,j),1)

          write(u,"(1x,'[FR-NORM-COORD]')")
          p = 1
          write(u,"(1x,'vibration',i24)") p
          if (pl > 0) then
                  print "(a,2i2,a)", &
                          " g-vector for ", mex(1), mex(2)," intersection:"
                  print "(a,f12.5)", "  ||g|| = ", ngv
                  print "(3f12.5)", gv(:,i,j)
          end if
          write(u,"(3f12.5)") gv(:,i,j)/ngv
          p = 2
          write(u,"(1x,'vibration',i24)") p
          if (pl > 0) then
                  print "(a,2i2,a)", &
                          " h-vector for ", mex(1), mex(2)," intersection:"
                  print "(a,f12.5)", "  ||h|| = ", nhv
                  print "(3f12.5)", hv(:,i,j)
          end if
          write(u,"(3f12.5)") hv(:,i,j)/nhv
          p = 3
          write(u,"(1x,'vibration',i24)") p
          write(u,"(3f12.5)") sv(:,i,j)/nsv
          return
  end subroutine print_molden_vecs


  ! print_molden_geom: print geometry in molden format, file is open.
  subroutine print_molden_geom(u, g, a, na)
          implicit none
          integer, intent(in) :: u, na
          character(3), dimension(na),intent(in) :: a
          real*8,  dimension(3 * na), intent(in) :: g
          integer :: i

          write(u,"(1x,'[FR-COORD]')")
          do i = 1, na
                write(u,"(1x,a3,3f13.6)") a(i), g((i-1)*3+1), g((i-1)*3+2),&
                        g((i-1)*3+3)
          end do
          return
  end subroutine print_molden_geom

  ! print_molden_freq: print frequencies to molden file, file is opened.
  subroutine print_molden_freq(u, f, nf)
          implicit none
          integer, intent(in) :: u, nf
          real*8,  dimension(nf), intent(in) :: f
          integer :: i

          write(u,"(1x,'[FREQ]')")
          do i = 1, nf
                write(u,"(f10.2)") f(i)
          end do

          return
  end subroutine print_molden_freq
end program testpoints

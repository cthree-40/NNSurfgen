!===========================================================================
!read COLUMBUS geom file and obtain geometry and atom info
SUBROUTINE readColGeom(gfile,ngeoms,na,atoms,anums,cgeom,masses)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)                           :: gfile
  INTEGER,INTENT(IN)                                    :: na
  INTEGER,INTENT(INOUT)                                 :: ngeoms
  CHARACTER(3),dimension(na),INTENT(INOUT)              :: atoms
  DOUBLE PRECISION,dimension(na),INTENT(INOUT)          :: anums,masses
  DOUBLE PRECISION,dimension(3*na,ngeoms),INTENT(INOUT) :: cgeom
  INTEGER                                               :: i,j,k,GUNIT,ios

  call FLUnit(GUNIT)
  !print "(A,I4,2A)","unit of geom file = ",GUNIT,", filename=",trim(adjustl(gfile))
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
    read(GUNIT,*,iostat=ios) atoms(j),anums(j),(cgeom(3*(j-1)+k,i),k=1,3),masses(j)
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
!===============================================================================
! efile     pathname of geometry file to be read in
! ngeoms    on input: maximum number of data entries that will be read in
!           on output: number of data entries that are read in from file
! nstates   total number of energies (states) per data. could be overriden with
!           state bounds stored inside the file.
! eners     energies read from the file
! st1       lower bound of energy data state index.  retrieved from file, with 
!           0 as the default when file does not specify states.
! st2       upper bound of energy data state index.  retrieved from file. the
!           default is nstates
SUBROUTINE readEner(efile,ngeoms,nstates,eners,st1,st2)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)                           :: efile
  INTEGER,INTENT(IN)                                    :: nstates
  INTEGER,INTENT(INOUT)                                 :: ngeoms
  INTEGER,INTENT(OUT)                                   :: st1,st2
  DOUBLE PRECISION,dimension(nstates,ngeoms),INTENT(INOUT)   :: eners
  INTEGER                                               :: i,j,EUNIT,ios, ne
  CHARACTER(6)  :: card

  call FLUnit(EUNIT)
  open(unit=EUNIT,file=trim(adjustl(efile)),access='sequential',form='formatted',&
    position='rewind',action='read',status='old',iostat=ios)
  if(ios/=0)then
    print *,'WARNING: readEner: cannot open file for read'
    ngeoms=0
    close(EUNIT)
    return
  end if
  read(EUNIT,"(6A)",advance='no') card
  if(card=='STATES') then
    read(EUNIT,*) st1,st2
    if(st1<1) st1=1
    if(st2>nstates) st2=nstates
  else
    rewind(EUNIT)
    st1 = 1
    st2 = nstates
  end if
  ne = st2-st1+1
  do i = 1,ngeoms
    read(EUNIT,*,IOSTAT=ios) (eners(j,i),j=st1,st2)
    if(ios/=0)then
        ngeoms = i-1
        close(EUNIT)
        return
    end if
  enddo
  close(EUNIT)
  return
END SUBROUTINE readEner
!===================================================================================
! efile     pathname of geometry file to be read in
! ngeoms    on input: maximum number of data entries that will be read in
!           on output: number of data entries that are read in from file
! nstates   total number of energies (states) per data. could be overriden with
!           state bounds stored inside the file.
! hd        energies read from the file
SUBROUTINE readHd(efile,ngeoms,nstates,hd)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)                           :: efile
  INTEGER,INTENT(INOUT)                                 :: ngeoms
  INTEGER,INTENT(IN)                                    :: nstates
  real*8,INTENT(INOUT)   :: hd(nstates*(nstates+1)/2,ngeoms)
  INTEGER                                               :: i,j,k,EUNIT,ios

  call FLUnit(EUNIT)
  open(unit=EUNIT,file=trim(adjustl(efile)),access='sequential',form='formatted',&
    position='rewind',action='read',status='old',iostat=ios)
  if(ios.ne.0) then
     print *,'WARNING: readHd: cannot open file for read'
     ngeoms=0
  end if

  k=nstates*(nstates+1)/2
  do i = 1,ngeoms
    read(EUNIT,*,IOSTAT=ios) (hd(j,i),j=1,k)
    if(ios.ne.0) then
        ngeoms = i-1
        close(EUNIT)
        return
    end if
  end do
  close(EUNIT)

  return
END SUBROUTINE readHd
!===================================================================================
SUBROUTINE readGrads(gfile,ngrads,na,cgrads)
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(IN)                           :: gfile
  INTEGER,INTENT(INOUT)                                 :: ngrads
  INTEGER,INTENT(IN)                                    :: na
  DOUBLE PRECISION,dimension(3*na,ngrads),INTENT(INOUT) :: cgrads
  INTEGER                                               :: i,j,k,GUNIT,ios

  call FLUnit(GUNIT)
  open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted',&
   action='read',position='rewind',status='old',iostat=ios)
  if(ios/=0)then
    ngrads = 0
    close(GUNIT)
    return
  end if

  do i = 1,ngrads
   do j = 1,na
    read(GUNIT,*,IOSTAT=ios)(cgrads(3*(j-1)+k,i),k=1,3)
    if(ios/=0)then
        ngrads=i-1
        close(GUNIT)
        return
    end if
   enddo
  enddo
  close(GUNIT)

  return
END SUBROUTINE readGrads
!===================================================================================
FUNCTION filename(s1,s2,grdptn,cpptn)
  IMPLICIT NONE
  INTEGER,INTENT(IN)          :: s1,s2
  CHARACTER(LEN=*),INTENT(IN) :: grdptn,cpptn
  CHARACTER(255)              :: filename
  CHARACTER(1)                :: st1,st2
  integer :: i

  write(st1,'(i1)')s1
  write(st2,'(i1)')s2

  if(s1.eq.s2)then
    filename=grdptn
    i=index(filename,'$')
    if(i>0)filename(i:i) = st1
    if(i==0)filename=''
  else
    filename=cpptn
    i=index(filename,'$')
    if(i>0)filename(i:i) = st1
    i=index(filename,'$')
    if(i>0)filename(i:i) = st2
    if(i==0)filename=''
  endif

  filename = trim(adjustl(filename))
END FUNCTION filename
!===================================================================================

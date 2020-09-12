module opttools
  implicit none
contains
  !*
  ! analysegeom: print out geometry information: bond lengths, angles, torsions
  ! Arguments:
  ! INTEGER			[in] natoms
  !				Number of atoms in the molecule
  ! DOUBLE PRECISION(3*natoms)	[in] geom
  !				geometry input
  ! CHARACTER*3(3*natoms)		[in] aname
  !				name of elements
  ! DOUBLE PRECISION(3*natoms)	[in] anum
  !				atomic number of each atom
  ! DOUBLE PRECISION(3*natoms)	[in] masses
  !				Nuclear mass of each atom
  ! DOUBLE PRECISION		[in] TLen
  !				Threshold for bond-length
  ! LOGICAL			[in] ShortList
  !				Whether to generate only a short list of
  !                             coordinates, instead those from all connected
  !                             atoms
  subroutine analysegeom(natoms,geom,aname,anum,masses,TLen,ShortList)
    implicit none
    integer, intent(in)          ::  natoms
    character*3,intent(in)       ::  aname(natoms)
    double precision,intent(in)  ::  anum(natoms),masses(natoms)
    double precision,intent(in)  ::  geom(3,natoms)
    double precision,intent(in)  ::  TLen
    logical,intent(in)           ::  ShortList
    double precision, parameter  ::  bohr2ang=0.529177249d0
    integer   ::  i,j,k,l
    double precision  ::  distmat(natoms,natoms)
    double precision  :: d,d1(3),d2(3),d3(3),cpd(3),p1(3),p3(3),ang
    double precision,external  ::   dnrm2
    logical :: hasOOP(natoms)
    print *,"Cartesian Geometries in Atomic Units"
    do i=1,natoms
      print "(x,a3,1x,f4.1,3F14.8,F14.8)",aname(i),anum(i),geom(:,i),masses(i)
    end do
    distmat = 0d0
    print "(/,A)","   Atom1   Atom2   Bond Length(Ang)"
    do i=1,natoms
      do j=i+1,natoms
        d = dnrm2(3,geom(:,i)-geom(:,j),1)*bohr2ang
        distmat(i,j) = d
        distmat(j,i) = d
        if(d<TLen)print "(2x,I5,3x,I5,3x,F12.8)",i,j,d
      end do
    end do
    
    print "(/,A)","   Atom1   Atom2   Atom3   Bond Angle (Degrees)"
    do i=1,natoms
      do j=1,natoms
        if(j==i .or. distmat(i,j)>TLen)cycle
        d1 = geom(:,j)-geom(:,i)
        d1 = d1/dnrm2(3,d1,1)
        do k=j+1,natoms
          if(k==i .or. distmat(i,k)>TLen)cycle
          d2 = geom(:,k)-geom(:,i)
          d2 = d2/dnrm2(3,d2,1)
          print "(2x,3(I5,3x),F12.4)",J,I,K, 90/Acos(0d0)* &
            ACOS(dot_product(d1,d2))
        end do!k
      end do !j
    end do   !i    
    
    print "(/,A)","   Atom1*  Atom2   Atom3   Atom4   OOP Angle (Degrees)"
    do i=1,natoms
      do j=1,natoms-2
        if(j==i)cycle
        d1 = geom(:,j)-geom(:,i)
        d1 = d1/dnrm2(3,d1,1)
        do k=j+1,natoms-1
          if(k==i .or. distmat(i,k)>TLen .or. distmat(i,j)>TLen  )cycle
          d2 = geom(:,k)-geom(:,i)
          d2 = d2/dnrm2(3,d2,1)
          do l=k+1,natoms
            if(l==i .or. distmat(i,l)>TLen )cycle
            d3 = geom(:,l)-geom(:,j)
            d3 = d3/dnrm2(3,d3,1)
            cpd(1) = d1(3)*(d2(2)-d3(2))+d2(3)*d3(2)-d2(2)*d3(3)+d1(2)*(d3(3)-d2(3))
            cpd(2) =-d2(3)*d3(1)+d1(3)*(d3(1)-d2(1))+d1(1)*(d2(3)-d3(3)) +d2(1)*d3(3)
            cpd(3) = d1(2)*(d2(1)-d3(1))+d2(2)*d3(1)-d2(1)*d3(2)+d1(1)*(d3(2)-d2(2))
            print "(2x,4(I5,3x),F12.4)",I,J,K,L, 90/Acos(0d0)* &
              asin((-d1(3)*d2(2)*d3(1)+d1(2)*d2(3)*d3(1)+d1(3)*d2(1)*d3(2)  &
              -d1(1)*d2(3)*d3(2)-d1(2)*d2(1)*d3(3)+d1(1)*d2(2)*d3(3))/      &
              dnrm2(3,cpd,1))
          end do! l
        end do!k
      end do !j
    end do   !i    
    
    print "(/,A)","   Atom1   Atom2   Atom3   Atom4   Torsion Angle (Degrees)"
    do i=1,natoms-1
      do j=1,natoms
        if(j==i .or. distmat(i,j)>TLen)cycle
        d1 = geom(:,i)-geom(:,j)
        d1 = d1/dnrm2(3,d1,1)
        do k=1,natoms
          if(k==i .or. k==j .or. distmat(j,k)>TLen  )cycle
          d2 = geom(:,k)-geom(:,j)
          d2 = d2/dnrm2(3,d2,1)
          p1=d1-dot_product(d2,d1)*d2
          p1=p1/dnrm2(3,p1,1)
          do l=i+1,natoms
            if(l==j .or. l==k .or. distmat(k,l)>TLen )cycle
            d3 = geom(:,l)-geom(:,k)
            d3 = d3/dnrm2(3,d3,1)
            p3=d3-dot_product(d2,d3)*d2
            p3=p3/dnrm2(3,p3,1)
            ang= 90/Acos(0d0)* acos(dot_product(p1,p3))
            print "(2x,4(I5,3x),F12.4)",I,J,K,L, ang 
          end do! l
        end do!k
      end do !j
    end do   !i    
    
  end subroutine analysegeom

  !*
  ! analysegeom2: this subroutine prints a geometry to a file
  !*
  subroutine analysegeom2(natoms,geom,aname,anum,masses,gflname)
    implicit none
    integer, intent(in)          ::  natoms
    character*3,intent(in)       ::  aname(natoms)
    character(len=300),intent(IN)      ::  gflname
    double precision,intent(in)  ::  anum(natoms),masses(natoms)
    double precision,intent(in)  ::  geom(3,natoms)
    double precision, parameter  ::  bohr2ang=0.529177249d0
    integer   ::  i, flun, ios
    
    ! Open new file
    flun = getflunit()
    open(unit=flun,file=gflname,status="unknown",position="rewind",iostat=ios)
    if (ios .ne. 0 ) then ! If file cannot be opened
      print "(1x,A)", "Could not open geom file. (analysegeom2)"
    end if
    do i=1,natoms
      write(unit=flun, fmt="(x,a3,1x,f4.1,3F14.8,F14.8)") &
        aname(i),anum(i),geom(:,i),masses(i)
    end do
    close(unit=flun)
    return
  end subroutine analysegeom2
  
  !*
  ! get an available UNIT index for input/output
  !*
  FUNCTION getFLUnit() RESULT(UNIT)
    IMPLICIT NONE
    INTEGER         :: UNIT
    integer         :: i
    logical         :: unitex,unitop
    UNIT=0
    do i=15,99999
      inquire(UNIT=i,EXIST=unitex,OPENED=unitop)
      if(unitex .and. .not. unitop)then
        UNIT=i
        exit
      end if!(unitex .and. .not. unitop)
    end do!i=10,99
    if(UNIT==0)stop "getFLUnit:  failed to find an available unit."
  END FUNCTION getFLUnit

  !*
  ! read COLUMBUS geometry file and obtain geometry and atom info
  !*
  SUBROUTINE readcolgeom(gfile,ngeoms,na,atoms,anums,cgeom,masses)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)                           :: gfile
    INTEGER,INTENT(IN)                                    :: na
    INTEGER,INTENT(IN)                                    :: ngeoms
    CHARACTER(3),dimension(na),INTENT(INOUT)              :: atoms
    DOUBLE PRECISION,dimension(na),INTENT(INOUT)          :: anums,masses
    DOUBLE PRECISION,dimension(3*na,ngeoms),INTENT(INOUT) :: cgeom
    INTEGER                                               :: i,j,k,GUNIT,ios
    INTEGER                                               :: ng
    
    GUNIT=getFLUnit()
    print "(A,I4,2A)","unit of geom file = ",GUNIT,", filename=",trim(adjustl(gfile))
    open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted',&
      status='old',action='read',position='rewind',iostat=ios)
    if(ios/=0)then
      print *,"Failed to open file [", trim(adjustl(gfile)),"]"
      close(GUNIT)
      return
    end if
    do i = 1,ngeoms
      ng = i
      do j = 1,na
        read(GUNIT,*,iostat=ios)atoms(j),anums(j),&
          (cgeom(3*(j-1)+k,i),k=1,3),masses(j)
        if(ios/=0)then
          print *,"  End of file encountered. i=",i,",j=",j
          ng=i-1
          close(GUNIT)
          return
        end if
      enddo
    enddo
    close(GUNIT)
    if (ng .ne. ngeoms) then
      print *, "***Warning: Fewer geometry entries found than requested.***"
    end if
    return
  end SUBROUTINE readColGeom
end module opttools

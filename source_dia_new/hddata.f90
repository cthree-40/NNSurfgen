!------------------------------------------------------------------------
MODULE HdDATA
  use nn_class
  IMPLICIT NONE
  type(nn) :: ANN !Feed-Forward Neural Network
  INTEGER :: ncoord   !total number of internal coords
  INTEGER :: nstates  !number of electronic states

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

CONTAINS
!------------------------------------------------------------------------
SUBROUTINE EvaluateHd(dispgeom,hmat)
      use progdata, only: abpoint, natoms
      IMPLICIT NONE
      TYPE(abpoint),INTENT(IN) :: dispgeom
      DOUBLE PRECISION,DIMENSION(nstates,nstates),INTENT(OUT) :: hmat
      double precision, allocatable :: xin(:)
      real*8 :: oop(noopc) ! user-defined oop coordinates
      integer :: R, i, j, k, o
      integer :: nrn ! neuron
      
      R=ANN%RX
      allocate(xin(R))

      !hmat
      xin=dispgeom%igeom(1:R)
      do i = 1, noopc
         oop(i)=dispgeom%igeom((ncoord-noopc)+i)
      end do

      call ANN%output(xin)

      hmat = 0d0
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
      
!      hmat(1,1) = ANN%y(1)
!      hmat(1,2) = ANN%y(2)*oop(1) + ANN%y(3)*oop(1)**3
!      hmat(2,1) = hmat(1,2)
!      hmat(1,3) = ANN%y(4)*oop(2)**2 + ANN%y(5)*oop(2)**4
!      hmat(3,1) = hmat(1,3)
!      hmat(2,2) = ANN%y(6)
!      hmat(2,3) = ANN%y(7)*oop(3) + ANN%y(8)*oop(3)**3
!      hmat(3,2) = hmat(2,3)
!      hmat(3,3) = ANN%y(9)

      return
END SUBROUTINE EvaluateHd
!------------------------------------------------------------------------
subroutine copy_wb
      call ANN%wb_copy()
      return
end subroutine copy_wb
!------------------------------------------------------------------------
subroutine rec_wb
      call ANN%wb_rec()
      return
end subroutine rec_wb
!------------------------------------------------------------------------
subroutine init_wb
      call ANN%wb_init(-1.d0,1.d0)
      return
end subroutine init_wb
!------------------------------------------------------------------------
subroutine savenet(outfl)
  implicit none
  character(len=*) :: outfl
  call ANN%savenn(outfl,1)
  return
end subroutine savenet
!------------------------------------------------------------------------
subroutine readnet(wfile)
  implicit none
  character(len=*) :: wfile
  call ANN%savenn(wfile,0)
  return

end subroutine readnet
!-----------------------------------------------------------------------
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
        print "(A,4i3,A,i3)", " Atoms: ", coord(i)%alist(:,j),&
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
  
1000 format(i10)         ! first line
1001 format(A)           ! coordinate titles
1002 format(3i6)         ! ctype, scaling type, perms
1003 format(4i6,i6)      ! atom list defs, coefficient
1004 format(2f6.3)       ! scaling coeficients
end subroutine readcoordinput
!-----------------------------------------------------------------------



END MODULE HdDATA
!=======================================================================

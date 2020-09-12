module hddata
  implicit none
  type, public :: cdef
          integer :: ctype
          integer :: nperm
          integer, dimension(4,3) :: alist
  end type cdef
  type, public :: hdform
          integer :: nc
          integer, dimension(4,10) :: terms
  end type hdform


  type(cdef), dimension(:),     allocatable :: coord

  type(hdform), dimension(:,:), allocatable :: hddef

  integer :: nstates
  integer :: lastlayer
  integer :: noopc
  
contains
  subroutine readcoordinput ()
    implicit none
    integer :: ncoords
    character(255), dimension(:), allocatable :: title
    integer, dimension(:,:), allocatable :: nc
    integer :: i, j, k, ios
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
            read(unit=11,fmt=1002,iostat=ios) coord(i)%ctype, coord(i)%nperm
            if (coord(i)%ctype .eq. -4) then
                    print "(A,i2,A)", "Coordinate ", i, " is a product coordinate."
            end if
            do j = 1, coord(i)%nperm
                    read(unit=11,fmt=1003,iostat=ios) coord(i)%alist(:,j)
            end do
    end do
    print *, "Coordinate definitions: "
    do i = 1, ncoords
            print "(A)", trim(adjustl(title(i)))
            print *, coord(i)%ctype, coord(i)%nperm
            do j = 1, coord(i)%nperm
                    print *, " Atoms: ", coord(i)%alist(:,j)
            end do
    end do
    ! Get state info (for block number)
    read(unit=11,fmt=1000,iostat=ios) nstates
    print *, " NSTATES = ", nstates
    if (.not. allocated(hddef)) allocate(hddef(nstates,nstates))
    allocate(nc(nstates,nstates))
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
    


    
    close(unit=11)
    return
    
1000 format(i10) ! first line
1001 format(A)   ! coordinate titles
1002 format(i5,i5) ! ctype, perms
1003 format(4i6) ! atom list defs
  end subroutine readcoordinput

  subroutine evaluatehd
    implicit none
    integer :: i, j, k
    integer :: order

    print *, "Evaluating hd form:"
    ! Loop over blocks of hd
    do i = 1, nstates
       do j = i, nstates
          print "('Block :  ', i2,',',i2)", j, i
          if (hddef(j,i)%nc .lt. 0) then
             do k = 1, abs(hddef(j,i)%nc)
                print *, "PIP"
             end do
          end if
          do k = 1, hddef(j,i)%nc
             do order = 1, 4
                if (hddef(j,i)%terms(order,k) .eq. 0) cycle
                print "(i2,'^',i2)", hddef(j,i)%terms(order,k), order
             end do
          end do
          
       end do
    end do
            
    return
  end subroutine evaluatehd
end module hddata
program main
  use hddata
  implicit none
  call readcoordinput
  call evaluatehd
end program main

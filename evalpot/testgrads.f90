!==========================================================================
program main
  use nnpes
  implicit none

  real*8, parameter :: au2cm=219474.63d0
  integer,parameter :: max_geom = 15000
  
  character(255) :: gfile
  
  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: h
  real*8, dimension(:), allocatable :: geom
  real*8, dimension(:), allocatable :: dispgeom, egrad
  real*8, dimension(:,:,:), allocatable :: cg, dcg

  integer :: i,j,k,ios

  real*8, parameter :: stepsize = 0.001
  real*8, external  :: dnrm2
  
  ! Get geometry filename from command line
  call get_command_argument(number=1,value=gfile,status=ios)
  if (ios .ne. 0) then
    print *, "Cannot get filename from command line. Exiting..."
    stop
  end if

  !initialize NN PES
  call pesinit()
  
  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))
  allocate(geom(3*natoms))

  ! Read in geometry
  call read_colgeom(gfile, geom, natoms)
  
  !NN PES evaluation
  call NNEvaluate(geom(:),e,cg,h,dcg)

  ! Print output
  print *,""
  print *,"Adiabatic energy(cm-1):"
  print "(5F24.15)", (e(i)*au2cm,i=1,nstates)
  print *,"Quasi-diabatic Hamiltonian energy(cm-1):"
  do i = 1, nstates
    print "(5F24.15)", h(i,:)*au2cm
  end do
  print *, ""
  do i = 1, nstates
    do j = i, nstates
      print "('Adiabatic gradients of states: ',2i3,E13.5)", i, j,&
        dnrm2(3*natoms,cg(1,i,j),1)
      print "(3f25.15)", cg(:,i,j)
      print *, ""
      print "('Diabatic gradients of states: ',2i3,E13.5)", i, j,&
        dnrm2(3*natoms,dcg(1,i,j),1)
      print "(3f25.15)", dcg(:,i,j)
      print *, ""
    end do
  end do

  ! Calculate gradient via centered difference
  allocate(dispgeom(3*natoms))
  allocate(egrad(3*natoms))
  print *, "Adiabatic gradients"
  do j = 1, nstates
      ! For state j
      egrad = 0d0
      do i = 1, 3*natoms
          dispgeom = geom
          dispgeom(i) = dispgeom(i) - stepsize
          call NNevaluate(dispgeom, e, cg, h, dcg)
          egrad(i) = egrad(i) - e(j) 
          dispgeom = geom
          dispgeom(i) = dispgeom(i) + stepsize
          call NNevaluate(dispgeom, e, cg, h, dcg)
          egrad(i) = egrad(i) + e(j)
          egrad(i) = egrad(i) / 2 / stepsize
      end do
      print *, "Energy gradient for state ", j
      print "(3f25.15)", egrad
      print *, ""
  end do
  print *, "Diabatic gradients"
  do j = 1, nstates
      do k = j, nstates
          ! For block j, k
          egrad = 0d0
          do i = 1, 3*natoms
              dispgeom = geom
              dispgeom(i) = dispgeom(i) - stepsize
              call NNevaluate(dispgeom, e, cg, h, dcg)
              egrad(i) = egrad(i) - h(j,k) 
              dispgeom = geom
              dispgeom(i) = dispgeom(i) + stepsize
              call NNevaluate(dispgeom, e, cg, h, dcg)
              egrad(i) = egrad(i) + h(j,k)
              egrad(i) = egrad(i) / 2 / stepsize
          end do
          
          print *, "Energy gradient for block ", j, k
          print "(3f25.15)", egrad
          print *, ""
      end do
  end do
  
  deallocate(egrad)
  deallocate(dispgeom)
  deallocate(e)
  deallocate(h)
  deallocate(cg)
  deallocate(dcg)
  deallocate(geom)
  
contains

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

  !-------------------------------------------------------------------
  ! read_colgeom: read geometry in columbus format
  subroutine read_colgeom (flnm, gm, na)
    implicit none
    character(255), intent(in) :: flnm
    integer, intent(in) :: na
    double precision, dimension(3,na), intent(out) :: gm
    integer :: i, flun, ios
    
    ! Open file
    flun = get_unit()
    open(unit=flun,file=flnm,status="old",iostat=ios)
    if (ios .ne. 0) stop "*** Can't open input geometry file. ***"
    ! Read geometry
    do i = 1, na
            read(flun,fmt="(11x,3f14.8,14x)") gm(1:3,i)
    end do
    close(flun)
    return
  end subroutine read_colgeom

end program main

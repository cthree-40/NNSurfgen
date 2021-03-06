!=================================================================================
MODULE makesurfdata
  use progdata, only: abpoint
  implicit none
  CHARACTER(255) :: gmfptn,enfptn,hdfptn ! naming pattern of input files.
  DOUBLE PRECISION :: eshift    ! uniform shift on ab initio energies

  !General options
  INTEGER :: npoints
  TYPE(abpoint),dimension(:),allocatable       :: dispgeoms

  ! weights for energy at every data point.
  double precision, dimension(:,:),allocatable :: eWeight
  double precision, dimension(10) :: energyT
  double precision, dimension(10) :: highEScale
  
! The following parameters control the construction of local coordinates
! Translational and rotational coordinates of the whole system will be removed,
! as well as those of the dissociated pieces if threshold and decaying parameter for
! the scaling functions are set properly.

! Eigenvectors and fitted energies
  real*8, allocatable :: ckl(:,:,:)
  DOUBLE PRECISION,dimension(:,:,:),allocatable     :: fitE

! Variables for Levenberg-Marquart
! Jacobian
integer :: npara, nse
real*8, allocatable :: Jac(:,:)
real*8, allocatable :: err(:)
real*8, allocatable :: jtjmui(:,:)
real*8, allocatable :: jtj(:,:)
real*8, allocatable :: jte(:)
real*8, allocatable :: BL(:) !for cholesky decomposition
real*8 :: mu
real*8, allocatable :: rmse(:)
real*8, allocatable :: rmsh(:)
integer :: epmax
integer :: nc
integer, allocatable :: trainings(:) ! which trainings to run
real*8 :: rhog

logical :: restart ! restart from WB-restart.txt

!Use GPU to compute dsyrk & cholesky decomposition
logical :: useGPU

CONTAINS
!----------------------------------------------------------------------------------
subroutine calc_jtj
    implicit none
    if(useGPU) then
      stop "** NO GPU **"
          !call gpu_dsyrk('U','N',npara,nse,1.d0,jac,npara,0.d0,jtj,npara)
    else
          call dsyrk('U','N',npara,nse,1.d0,jac,npara,0.d0,jtj,npara)
    end if
    return
end subroutine calc_jtj
!----------------------------------------------------------------------------------
subroutine calc_jte
    implicit none
    call dgemv('N',npara,nse,1.d0,jac,npara,err,1,0.d0,jte,1)
    return
end subroutine calc_jte
!----------------------------------------------------------------------------------
subroutine calc_jtjmui
      implicit none
      integer :: i
      jtjmui=jtj
      forall(i=1:npara) jtjmui(i,i)=jtjmui(i,i)+mu
      return
end subroutine calc_jtjmui
!----------------------------------------------------------------------------------
!*
! printerrorE: print energy error for each point at the end of a calculation.
!*
subroutine printerrorE
  use progdata, only: ev2cm
  use hddata
  implicit none
  integer :: i, j
  
  print *, ""
  print *, "Energy Error: Ab-initio vs NN-Surfgen adiabatic energy (cm-1)"
  print *, ""
  do i = 1, npoints
    write(*,"('PT:',i10)",advance="no") i
    do j = 1, nstates
      write(*,"(F13.2,F13.2)",advance="no") dispgeoms(i)%energy(j,j)*ev2cm, &
        (dispgeoms(i)%energy(j,j)-fitE(i,j,j))*ev2cm
    end do
    write(*,"('')")
  end do
  print *, ""
  return
end subroutine printerrorE
!----------------------------------------------------------------------------------
end module makesurfdata
!----------------------------------------------------------------------------------
subroutine update_wb(info)
      use hddata, only: ANN 
      use makesurfdata, only:npara,jtjmui,jte,BL,useGPU
      implicit none
      integer, intent(out) :: info

      BL=jte
      if(useGPU) then
        stop "*** No GPU ***"
            !call cuda_cholesky(npara,jtjmui,BL,info)
      else
            call cholesky(npara,jtjmui,BL,info)
      end if

      if(info .ne. 0) then
            print*, 'Error in Cholesky decomposition!'
            return
      end if

      call ANN%wb_update(BL)

      return
end subroutine update_wb
!----------------------------------------------------------------------------------
subroutine cholesky(L,jtjmui,BL,info)
      implicit none
      integer, intent(in) :: L
      real*8, intent(inout) :: jtjmui(L,L)
      real*8, intent(inout) :: BL(L)
      integer, intent(out) :: info

      call dpotrf('U',L,jtjmui,L,info)
      if(info .eq. 0) call dpotrs('U',L,1,jtjmui,L,BL,L,info)

      return
end
!----------------------------------------------------------------------------------
!Read in geometry and ab initio data for a set of points which will be used
!in the Hd fitting procedure.
SUBROUTINE readdisps()
  use hddata
  use progdata,only:natoms,printlvl,au2cm,au2ev
  use makesurfdata
  IMPLICIT NONE
  INTEGER :: i,j,k,l,s,ios
  DOUBLE PRECISION,allocatable,dimension(:,:) :: hd,cgeoms,eners,bmat
  character(255),external :: filename
  character(255) :: infile
  character(3),dimension(natoms) :: atoms
  double precision,dimension(natoms) :: anums,masses
  !ptinfile : number of points in a geometry input file.
  !lb, ub   : energies for states between lb and ub are given in energy file
  integer :: ptinfile,lb,ub
  double precision :: NaN
  double precision,external :: dnrm2

  if (.not. allocated(eWeight)) allocate(eWeight(npoints,nstates))
  allocate(cgeoms(3*natoms,npoints))
  allocate(eners(nstates,npoints))
  allocate(hd(nstates*(nstates+1)/2,npoints))
  allocate(bmat(ncoord,3*natoms))
  NaN = 0d0
  NaN = NaN/NaN

  if(printlvl>0) print '(4X,A,I5,A)','Reading',npoints,' displacements'
  if(allocated(dispgeoms)) deallocate(dispgeoms)
  allocate(dispgeoms(npoints))
  do j=1,npoints
     allocate(dispgeoms(j)%cgeom(3*natoms))
     allocate(dispgeoms(j)%igeom(ncoord))
     allocate(dispgeoms(j)%energy(nstates,nstates))
     allocate(dispgeoms(j)%hd(nstates,nstates))
  end do

  !read geometries
  infile = trim(adjustl(gmfptn))
  if(printlvl>1)print 1001,'- Reading geometry input <'//trim(infile)//'>'
  ptinfile=npoints
  call readColGeom(infile,ptinfile,natoms,atoms,anums,cgeoms,masses)
  if(printlvl>1)print 1002,"found ",ptinfile, " geometries."
  if(npoints .ne. ptinfile) stop 'npoints .ne. ptinfile in readColGeom!'
 
  do j=1,ptinfile
      dispgeoms(j)%id=j
      dispgeoms(j)%cgeom(1:3*natoms)=cgeoms(1:3*natoms,j)
  enddo!j=1,ptinfile

  !read energy data
  infile = trim(adjustl(enfptn))
  if(printlvl>1)print 1001,'- Reading energy input <'//trim(infile)//'>'
  ptinfile=npoints
  call readEner(infile,ptinfile,nstates,eners,lb,ub)
  if(npoints .ne. ptinfile) stop 'npoints .ne. ptinfile in readEner!'
  eners=(eners+eshift)*au2ev !to eV
  print 1002,"found ",ptinfile," energy data entries."

  do j=1,ptinfile
      dispgeoms(j)%energy = NaN
      dispgeoms(j)%energy(lb:ub,lb:ub) = 0d0 
      do k=lb,ub
          dispgeoms(j)%energy(k,k) = eners(k,j)
      end do
      dispgeoms(j)%lb = lb
      dispgeoms(j)%ub = ub
      print *, "Point, lb, ub = ", j, lb, ub
  enddo

  ! CLM: Print out energy information
  print *, " Energy information: "
  write(*, "(2x,'Point',6x)", advance="no")
  do j = 1, nstates
          write(*, "(5x,'State ',i1,3x)", advance="no") j
  end do
  print *, ""
  do j = 1, ptinfile
          write(*, "(2x,i6,2x)", advance="no") j
          do k = 1, nstates
                  write(*, "(F15.6)", advance="no") eners(k,j)
          end do
          print *, ""
  end do
  print *, ""
  !read in diabats data
  infile = trim(adjustl(hdfptn))
  if(printlvl>1)print 1001,'- Reading diabats input <'//trim(infile)//'>'
  ptinfile=npoints
  call readHd(infile,ptinfile,nstates,hd)
  if(npoints .ne. ptinfile) stop 'npoints .ne. ptinfile in readHd!'
  hd=hd*au2ev  !to eV

  do j=1,ptinfile
      s=0
      do k=1,nstates
        do l=k,nstates
          s=s+1
          dispgeoms(j)%hd(k,l)=hd(s,j)
          if(k.ne.l) dispgeoms(j)%hd(l,k)=dispgeoms(j)%hd(k,l)
        end do
      end do
  enddo

  !process the data that have been read in
  do l = 1,npoints
    if(printlvl>0) print *,""
    if(printlvl>0) print *,"   Constructing Internal Coordinates at point ",l," :"
    call buildWBmat(natoms,ncoord,dispgeoms(l)%cgeom,dispgeoms(l)%igeom,Bmat)
    if(printlvl>0)then
        print *,"  Internal Coordinates"
        print "(10E11.3)",dispgeoms(l)%igeom
    end if
  end do!l

  !process data point energy data for eWeight
  eWeight = 1d0
  do l = 1, npoints
      write(*,"('Point ',i5, ' weights: ')",advance="no") l
      do i = 1, nstates
          k = 0
          do while(dispgeoms(l)%energy(i,i) .gt. energyT(k+1))
              k=k+1
              if (k .eq. 10) exit
          end do
          if (k .gt. 0) eWeight(l,i) = highEScale(k)
          write(*,"(f12.5,'(',f5.3,')')",advance="no") &
              dispgeoms(l)%energy(i,i), eWeight(l,i)
      end do
      write(*,"('')")
  end do
  
1001 format(4X,A)
1002 format(10X,A,I6,A)

  return
end SUBROUTINE readdisps
!==================================================================================
subroutine initmakesurf
  use hddata, only:nstates,ANN
  use makesurfdata
  implicit none

  !For eigenvectors
  if(allocated(ckl)) deallocate(ckl)
  allocate(ckl(npoints,nstates,nstates))
  if(allocated(fitE)) deallocate(fitE)
  allocate(fitE(npoints,nstates,nstates))

  !For L-M
  npara=ANN%L
  allocate(jtj(npara,npara),jtjmui(npara,npara),jte(npara),BL(npara))

  return
end
!=================================================================================

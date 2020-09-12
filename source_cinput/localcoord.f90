  !-------------------------------------------------------------------------------
  ! Construct a unitary transformation L that transforms cartesian
  ! coordinate into an orthongal system defined by the eigenvectors
  ! of matrix B^T.B, where B is Wilson's B matrix.
  ! Eigenvectors will be ordered with descending eigenvalues.  Coordinates
  ! with an eigenvalue larger than intGT will be considered internal.
  ! Number of internal degrees of freedom will be stored in field nvibs 
  !---------------------------------------------------------------------
  ! This subroutine is used to construct a local coordinate system that
  ! is used for the fitting.  All gradients and couplings will be fitted
  ! with this coordinate, which is made form linear combinations of the 
  ! 3*N cartesian coordinates, but has translation/rotation and dissappearing
  ! vibrations from dissociation separated from vibrational coordinates.
  ! The coordinates will also reflect the symmetry of the point, because 
  ! right singular vectors of the B matrix is used, which carries the full
  ! symmetry of the point as long as the internal coordinates are symmetric
  ! (which is always the case).
SUBROUTINE makeLocalIntCoord(ptdata,nstates,intGT,intGS,nvibs)
    USE hddata, ONLY: ncoord
    USE progdata, ONLY: abpoint,natoms,printlvl
    IMPLICIT NONE
    double precision, parameter  :: degcap=1d-3 ! threshold for degeneracy
    TYPE(abpoint),intent(inout)  :: ptdata
    INTEGER,intent(IN)           :: nstates,nvibs
    DOUBLE PRECISION,intent(IN)  :: intGT,intGS
    integer                      :: i,n1,n2
    double precision,dimension(3*natoms,3*natoms) :: btb      !B^T.B
    double precision,dimension(3*natoms,3*natoms) :: btb0
    double precision,dimension(3*natoms)          :: ev,grad  !eigenvalues of BtB
    double precision,dimension(45*natoms*natoms)  :: scr
    double precision             ::  nrmerr, alpha
    integer                      ::  INFO
    character(3)                 ::  str
    !Number of degeneracies among coordinates
    integer     :: ndeg
    !The list of coordinate ranges that are degenerate.  
    integer     :: degcoords(nvibs,2) 
    !The starting index of the degenerate block. Coordinates from this index to
    !index i-1 are all degenerate.
    integer     ::  degstart
    !bmatrix with new coordinates
    double precision  :: bnew(ncoord,3*natoms)
    double precision,external :: dnrm2

    !Construct matrix B^T.B and get its eigenvalues
    call DSYRK('U','T',3*natoms,ncoord,1d0,ptdata%cbmat,ncoord,0d0,btb,3*natoms)
    call DSYEV('V','U',3*natoms,btb,3*natoms,ev,scr,45*natoms*natoms,INFO)
    if(INFO .ne. 0) stop 'DSYEV FAILED IN MAKELOCALINTCOORD!'
    !reorder eigenvalues and eigenvectors in descending order
    btb0=btb
    do i=1,3*natoms
      !ptdata%eval(i)=ev(3*natoms+1-i)
      btb(:,i)=btb0(:,3*natoms+1-i)
    end do
     
    ndeg=0
    degstart=1
    ptdata%bmat=matmul(ptdata%cbmat,btb)
    ptdata%lmat=btb
    do i=1,3*natoms
      ptdata%eval(i)=ev(3*natoms+1-i)
    end do

    if(printlvl>1)then
      print *,"  Eigenvalues of B^T.B"
      print "(10E11.3)",ptdata%eval
    end if
    if (printlvl>3)then
      print *,"  Eigenvectors of B^T.B"
      do i=1,3*natoms
         print "(12E11.3)",btb(:,i)
      end do
    end if
    ptdata%scale=dble(1)
    do i=1,3*natoms
      if(ev(3*natoms+1-i)>intGT) then
        !update the list of degenerate modes
        if(i>1)then
          if(abs(ptdata%eval(i)-ptdata%eval(i-1))>degcap) then
          !this mode is not degenerate with the previous  
             if(degstart<i-1)then !there are more than one degenerate coordiante
               ndeg=ndeg+1
               degcoords(ndeg,1)=degstart
               degcoords(ndeg,2)=i-1
             end if
             degstart=i
          end if
        end if
        ptdata%nvibs=i
        if(i>nvibs) then
        !Degrees of freedom greater than 3N-6.  Something must have gone
        !wrong if it ever gets here.
          print *,"PT",ptdata%id,"NVIBS(PT) = ",i,", natoms = ",natoms
          print *, "GEOM:"
          print "(3F11.7)",ptdata%CGEOM

          print *, "B Matrix "
          do n1=1,ncoord
            print "(12E11.2)",ptdata%bmat(n1,:)
          end do
          print *,"EV:",ev
          stop "Error: Number of degrees of freedom exceeded theoretical maximum"
        end if
        if(ev(3*natoms+1-i)<intGS) ptdata%scale(i)=ev(3*natoms+1-i)/intGS
      else
        ptdata%scale(i:)=dble(0)
        if(degstart<i-1)then
          ndeg=ndeg+1
          degcoords(ndeg,1)=degstart
          degcoords(ndeg,2)=i-1
        end if
        exit
      end if
    end do !i=1,3*natoms

    if(printlvl>1) then
      print "(A,I7)","     Number of internal degrees of freedom : ",ptdata%nvibs
      if(ndeg>0) print "(A,I7)","     Number of degenerate internal motions : ",ndeg
    end if

    !construct local coordinates
    do n1=1,nstates
      do n2=1,nstates
        do i=1,3*natoms
          grad(i)=dot_product(btb(:,i),ptdata%grads(:,n1,n2)) 
        end do
        ptdata%grads(:,n1,n2)=grad
        nrmerr=dnrm2(3*natoms-ptdata%nvibs,grad(ptdata%nvibs+1:),1)
        if(nrmerr*0d0 .eq. 0d0 .and. nrmerr.eq.nrmerr)then
            if(printlvl>1.and.nrmerr>1D-10.or.printlvl>2) &
                  print 1000,"residual norm of external gradients for block",&
                  n1,n2,dnrm2(3*natoms-ptdata%nvibs,grad(ptdata%nvibs+1:),1)
            if(printlvl>0.and.nrmerr>1D-5) print "(4X,12F9.4)",grad(ptdata%nvibs+1:)
        end if
      end do
    end do !n1=1,nstates

    !Make linear combinations of degenerate internal coordinates to symmetrize them
    if(ndeg>0) call dgemm('n','n',ncoord,3*natoms,3*natoms,1d0,ptdata%bmat,ncoord, &
                 btb,3*natoms,0d0,bnew,ncoord) 
    do i=1,ndeg
      if(printlvl>2) print "(3X,2(A,I4))",&
            "Symmetrizing coordinate range ",degcoords(i,1)," to ",degcoords(i,2)
      if(printlvl>3)then
        do n1=degcoords(i,1),degcoords(i,2)
          print *,"  Internal contribution of coordinate ",n1
          print "(10(x,F10.5))",bnew(:,n1)
        end do
      end if
    end do

    !Print out information about the final coordinate system
    if(ptdata%nvibs<3*natoms-6.or.printlvl>1)then
      write(str,'(I3)') 3*natoms
      print "(7X,A,"//trim(str)//"F7.3)",&
            "scaling factors of local coordinates: ",ptdata%scale
      if(ptdata%nvibs<3*natoms-6) then
        print *,"      Local coordinate system has reduced dimensionalities."
        if(printlvl>2)then
          print *,"   Cartesian vectors of non-vanishing coordinates:"
          do n1=1,ptdata%nvibs
             print "(A,I3,A,F9.5)","    Fitting coordinate ",n1,&
                   " ,eigenvalue=",ptdata%eval(n1)
             print "(3F12.7)",btb(:,n1)
          end do
          print *,"   Cartesian vectors of vanishing coordinates:"
          do n1=ptdata%nvibs+1,3*natoms
             print "(A,I3,A,F9.5)","    Fitting coordinate ",n1,&
                   " ,eigenvalue=",ptdata%eval(n1)
             print "(3F12.7)",btb(:,n1)
          end do
          print *,"Overlaps of non-vanishing coordinates with nascent internals:"
          do n1=1,ptdata%nvibs
             print "(A,I3,A,F9.5)","    Fitting coordinate ",n1,&
                   " ,eigenvalue=",ptdata%eval(n1)
             print "(12F9.4)",(dot_product(btb(:,n1),ptdata%bmat(n2,:)),n2=1,ncoord)
          end do
          print *,"Overlaps of vanishing coordinates with nascent internals:"
          do n1=ptdata%nvibs+1,3*natoms
             print "(A,I3,A,F9.5)","    Fitting coordinate ",n1,&
                   " ,eigenvalue=",ptdata%eval(n1)
             print "(12F9.4)",(dot_product(btb(:,n1),ptdata%bmat(n2,:)),n2=1,ncoord)
          end do
        end if
      end if
    end if

  1000 format(8x,A,I3,",",I3," : ",E11.4)

  return
END SUBROUTINE makeLocalIntCoord
!====================================================================================

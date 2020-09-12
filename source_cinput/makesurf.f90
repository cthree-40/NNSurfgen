!=================================================================================
MODULE makesurfdata
  use progdata, only: abpoint
  implicit none
!derived type to identify a certain block of a point
  TYPE tblk
    INTEGER :: point
    INTEGER :: i,j
  END TYPE
!derived type for eq lists
  TYPE TEqList
   type(tblk),dimension(:),allocatable         :: List
   INTEGER                                     :: length
  END TYPE
!specifies the presence of data at each point
  LOGICAL,dimension(:,:),allocatable :: hasEner
  !specifies the presence of data at each point
  LOGICAL,dimension(:,:,:),allocatable :: hasGrad
  CHARACTER(255) :: gmfptn,enfptn,grdfptn,cpfptn  ! naming pattern of input files.

!General options
  DOUBLE PRECISION :: eshift    ! uniform shift on ab initio energies
  INTEGER :: npoints
  TYPE(abpoint),dimension(:),allocatable :: dispgeoms

  type(TEqList) :: exclEner, exclGrad

! number of least squares fitting equations
  integer :: neqs
  integer :: ncons

!Weights for least squars equations
  DOUBLE PRECISION,dimension(:),allocatable         :: weight
!Weights for each gradient and coupling blocks, caused by energy differences
!and high energy.
  DOUBLE PRECISION,dimension(:,:,:),allocatable     :: gWeight
  !Weights for each energy at every data point, caused by high energy thresholds
  DOUBLE PRECISION,dimension(:,:,:),allocatable     :: eWeight
  DOUBLE PRECISION,dimension(:),allocatable :: ptWeights
  real*8 :: ptwthreshold

!These parameters control the automatic downscaling of high energy points
!EnergyT      Threshold energy above which gradients/energies will be scaled down
!HighEScale   The scaling factor for these high energy data 
  DOUBLE PRECISION,dimension(10) :: energyT, highEScale
  DOUBLE PRECISION,dimension(10) :: energyT_en, highEScale_en

  ! energy above which gradients will be automatically removed from fitting
! equations
  DOUBLE PRECISION  :: gradcutoff
! energy above which couplings with this state will be removed from fitting
! equations.  when <=0, value of gradcutoff is used instead
  DOUBLE PRECISION  :: cpcutoff

! weights for equations of reproducing nonadiabatic coupling
  DOUBLE PRECISION :: nrmediff, ediffcutoff
  
!weights for equations of reproduction of energy, energy gradients,
!and non-adiabatic coupling
  DOUBLE PRECISION :: w_energy,w_grad,w_fij
  double precision :: rho_energy,rho_grad,rho_fij

!scaling factors to exact equations
  logical,dimension(:,:,:),allocatable :: incgrad
  logical,dimension(:,:,:),allocatable :: incener

!map of equations
  integer,dimension(:,:),allocatable :: EqMap

!threshold below which intersection adapted coordinates will be used
  real*8 :: deg_cap
  real*8 :: gorder
  real*8 :: gcutoff   !gradient cutoff
  logical :: deggrdbinding

! The following parameters control the construction of local coordinates
! Translational and rotational coordinates of the whole system will be removed,
! as well as those of the dissociated pieces if threshold and decaying parameter for
! the scaling functions are set properly.

! useInternal  :   transform ab initio gradients to an orthogonal local internal
! coordinates defined as eigenvectors of B^T.B, where B is Wilson's B matrix.
! intGradT     :   Threshold for eigenvalue of B^T.B higher than which
! the coordinate is considered internal
! intGradS     :   Threshold for eigenvalue of B^T.B lower than which the coordinate 
! will be scaled with a factor eval/intGrad, where eval is the eigenvalue
  DOUBLE PRECISION :: intGradT,intGradS 

! nvibs:  number of coordinates along which gradients will be taken.
! 3*natoms when useIntGrad=.FALSE.,  3*natoms-5 when useIntGrad=.TRUE. 
! Individual data point will have a number of internal degree of freedom
! equal to or lower than nvibs, stored in the field abpoint%nvibs
  INTEGER :: nvibs 

! Eigenvectors and fitted energies, gradients and couplings
  real*8, allocatable :: ckl(:,:,:)
  !Hd predictions for all data points
  DOUBLE PRECISION,dimension(:,:,:,:),allocatable   :: fitG
  DOUBLE PRECISION,dimension(:,:,:),allocatable     :: fitE
! Phase List for fixing phase
  integer, allocatable :: phaseList(:,:)
! Permutation List
  integer, dimension(:,:), allocatable :: pmtList

! Variables for Levenberg-Marquart
real*8, allocatable :: Jac(:,:)
real*8, allocatable :: err(:)
real*8, allocatable :: jtjmui(:,:)
real*8, allocatable :: jtj(:,:)
real*8, allocatable :: jte(:)
real*8, allocatable :: BL(:) !for cholesky decomposition
real*8, allocatable :: coef(:), coefc(:) !nonlinear parameters
real*8 :: mu
real*8 :: lambda !regularization
integer :: epmax
real*8 :: rmsee,mue,rmsec,rmseg

!Use GPU to compute dsyrk & cholesky decomposition
logical :: useGPU

CONTAINS
!----------------------------------------------------------------------------------
!determine if each equation will be include / excluded / fitted exactly
SUBROUTINE getPtList()
  use hddata, only: nstates
  use progdata, only: printlvl
  IMPLICIT NONE
  INTEGER i,j,s1,s2,k,l

  if(printlvl>0)print *,"Processing data inclusion/exclusion/exact fit data"
  if(allocated(incgrad)) deallocate(incgrad)
  if(allocated(incener)) deallocate(incener)

  allocate(incgrad(npoints,nstates,nstates))
  allocate(incener(npoints,nstates,nstates))

  !Generate destiny table for each piece of ab initio data
  incgrad=.false.
  incener=.false.
  do i=1,npoints
    if (ptWeights(i)>1D-8) then
      do s1=1,nstates
        do s2=1,s1
          if ((dispgeoms(i)%energy(s1,s1)<gradcutoff) .and. &
            ((s1 .eq. s2) .or. (dispgeoms(i)%energy(s1,s1) .lt. cpcutoff))) then
            incgrad(i,s1,s2)=hasGrad(i,s1,s2)
          end if
        end do
        incener(i,s1,s1)=hasEner(i,s1)
      end do!s1=1,nstates
      if (i == 1) then
        print *, "INCENER=", incener(i,:,:)
        print *, "INCGRAD=", incgrad(i,:,:)
        print *, "HASENER=", hasEner(i,:)
        print *, "HASGRAD=", hasGrad(i,:,:)
      end if
      !add off-diagonal elements between states within degeneracy groups
      do j=1,dispgeoms(i)%ndeggrp
        s1 = dispgeoms(i)%deg_groups(j,1)
        s2 = s1-1+dispgeoms(i)%deg_groups(j,2)
        !If all energies are available in the degeneracy group, add off-diagonal
        !elements to the fitting set because these states are subject to rotations
        !to form intersection adapted coordinates
        if(all(hasEner(i,s1:s2))) incener(i,s1:s2,s1:s2)=.true.
      end do !j
    end if!ptWeights(i)>1D-8
  end do!i=1,npoints

  do i=1,exclEner%length
    s1=exclEner%List(i)%i
    s2=exclEner%List(i)%j
    incener(exclEner%List(i)%point,s1:s2,s1:s2)=.false.
  end do

  do i=1,exclGrad%length
    incgrad(exclGrad%List(i)%point,exclGrad%List(i)%i,exclGrad%List(i)%j)=.false.
    incgrad(exclGrad%List(i)%point,exclGrad%List(i)%j,exclGrad%List(i)%i)=.false.
  end do

  !exclude all gradient and couplings if not all the data are present
  !in a degenerate block if deggrdbinding is set to true
  if(deggrdbinding) then
    do i=1,npoints
      do j=1,dispgeoms(i)%ndeggrp
        s1 = dispgeoms(i)%deg_groups(j,1)
        s2 = s1-1+dispgeoms(i)%deg_groups(j,2)
stloop: do k = s1,s2
          do l = k,s2
            if(.not.(incgrad(i,k,l).or.incgrad(i,l,k))) then
              incgrad(i,:,s1:s2)=.false.
              incgrad(i,s1:s2,:)=.false.
              exit stloop
            end if
          end do !l
        end do stloop!k
      end do !j
    end do! i=1,npoints
  end if

  return
end subroutine getPtList
!----------------------------------------------------------------------------------
!This subroutine reads in an array of points, then output a list of equations.
!In the list are the point ID/block ID/gradient direction of each equation.
!This can be used to generate either exact or least square equation list
!maxEqs          :  Largest possible number of equations
!nEqs            :  Number of equations actually generated
!lseMap(maxEqs,4):  The first nEqs rows record the (/point ID,state1,state2,grad ID/)
!                     of the equations generated. gradID=0 for energies
SUBROUTINE makeEqMap(maxEqs,lseMap,wvec)
  use hddata, only: nstates
  use progdata, only: au2ev
  implicit none
  integer,intent(in) :: maxEqs
  integer,dimension(MaxEqs,4),intent(out) :: lseMap
  double precision,dimension(MaxEqs),intent(out) :: wvec
  integer :: i,j,k,s1,s2,ng
  character(4),dimension(0:1) :: fatesymb
  character(4),dimension(nstates) :: enerfate
  character(7),dimension(nstates*(nstates+1)/2) :: gradfate
  character(3) :: str1, str2  ! # of ener and grad data
  double precision :: ediff
  
  !Output point destiny table
  write(str1,"(I3)")  nstates
  write(str2,"(I3)")  nstates*(nstates+1)/2

  print *, "Destiny of ab initio data at each data point"
  print *, "  LS = Least Squares Fit"
  print *, "  EX = Excluded from Fit"
  ng = 1
  do s1 = 1,nstates
    write(enerfate(s1),"(X,I2,X)") s1
    do s2= 1,s1
      write(gradfate(ng),"(x,I2,',',I2,x)") s1,s2
      ng= ng + 1
    end do !s2
  end do !s1
  print *," POINT ",enerfate,gradfate

  fatesymb(0) = " LS "
  fatesymb(1) = " EX "
  do i=1, npoints
    ng = 1
    do s1 = 1,nstates
      enerfate(s1) = fatesymb(0)
      if(.not.incener(i,s1,s1))  enerfate(s1) = fatesymb(1)
      do s2=1,s1
        gradfate(ng) = " "//fatesymb(0)//"  "
        if(.not.incgrad(i,s1,s2)) gradfate(ng) = " "//fatesymb(1)//"  "
        ng = ng+1
      end do !s2
    end do !s1
    print "(2X,I5,2X,"//trim(adjustl(str1))//"A4,"//trim(adjustl(str2))//"A7)",&
           i,enerfate,gradfate
  end do!i

  !make equations
  nEqs=0
  if(allocated(gWeight)) deallocate(gWeight)
  allocate(gWeight(npoints,nstates,nstates))
  gWeight=1.d0
  if(allocated(eWeight)) deallocate(eWeight)
  allocate(eWeight(npoints,nstates,nstates))
  eWeight=w_energy
  do i=1,npoints
    !calculate weights of gradient and coupling blocks caused by energy difference
    do s1=1,nstates
      ! Weights for energy gradients
      gWeight(i,s1,s1)=w_grad
      !determine the bracket of current energy
      k=0
      do while(dispgeoms(i)%energy(s1,s1) .gt. energyT(k+1))
        k=k+1
        if(k .eq. 10) exit
      end do
      gWeight(i,s1,s1) = gWeight(i,s1,s1)*highEScale(k+1)
      !Weights for derivative couplings
      do s2=1,s1-1
        gWeight(i,s1,s2)=w_fij
        !determine the bracket of ab initio energy
        k=0
        do while(dispgeoms(i)%energy(s1,s1)+dispgeoms(i)%energy(s2,s2) .gt. &
                 2.d0*energyT(k+1))
          k=k+1
          if(k .eq. 10) exit
        end do
        gWeight(i,s1,s2)=gWeight(i,s1,s2)*highEScale(k+1)
        if (nrmediff>0) then
          ediff=abs(dispgeoms(i)%energy(s1,s1)-dispgeoms(i)%energy(s2,s2))
          ediff=(ediff+ediffcutoff)/nrmediff
          gWeight(i,s1,s2)=gWeight(i,s1,s2)/ediff
        end if
        gWeight(i,s2,s1)=gWeight(i,s1,s2)
      end do!s2
    end do!s1

    !calculate the weights of energy equations
    do s1=1,nstates
      k=0
      do while(dispgeoms(i)%energy(s1,s1) .gt. energyT_en(k+1))
        k=k+1
        if(k .eq. 10) exit
      end do
      eWeight(i,s1,s1)=w_energy*highEScale_en(k+1)
      do s2=s1+1,nstates
        eWeight(i,s1,s2)=eWeight(i,s1,s1)
        eWeight(i,s2,s1)=eWeight(i,s1,s1)
      end do !s2=1,nstates
    end do!s1

    !create gradient and coupling component of weight vector
    if(ptWeights(i) .lt. ptwthreshold) cycle

    !Create equations for gradients and derivative couplings
    do s1=1,nstates
      do s2=1,nstates
        if(incgrad(i,s1,s2)) then
          do j=1,dispgeoms(i)%nvibs
            nEqs=nEqs+1
            if(nEqs>maxEqs) stop 'makeEqMap: neqs>maxeqs'
            lseMap(nEqs,:)=(/i,s1,s2,j/)
            wvec(nEqs)=ptWeights(i)*gWeight(i,s1,s2)
          end do!j=1,dispgeom()%nvibs
        end if!(incgrad(i,s1,s2))
      end do !s2=1,s1
    end do !s1=1,nstates

  end do !i=1,npoints

  !create the component of weight vector for energy equations
  do i=1,npoints
    if(ptWeights(i) .lt. ptwthreshold) cycle
    do s1=1,nstates
      do s2=1,nstates
        if(incener(i,s1,s2)) then
          nEqs=nEqs+1
          if(nEqs>maxEqs)stop 'makeEqMap: neqs>maxeqs'
          lseMap(nEqs,:)=(/i,s1,s2,0/)
          wvec(nEqs)=ptWeights(i)*eWeight(i,s1,s2)
        end if
      end do !s2=1,nstates
    end do !s1=1,nstates
  end do!i=1,npoints

  1000 format(8X,"Point",I5,", St",2I3,", dir",I3,", grad=",E10.3)
  1001 format(8X,3F14.10)

END SUBROUTINE makeEqMap
!----------------------------------------------------------------------------------
subroutine calc_jtj
    implicit none
    if(useGPU) then
      stop "*** No GPU on NERSC! Please set useGPU=.false. ***"
          !call gpu_dsyrk('U','N',ncons,neqs,1.d0,jac,ncons,0.d0,jtj,ncons)
    else
          call dsyrk('U','N',ncons,neqs,1.d0,jac,ncons,0.d0,jtj,ncons)
    end if
    return
end subroutine calc_jtj
!----------------------------------------------------------------------------------
subroutine calc_jte
    implicit none
    jte=lambda*coef
    call dgemv('N',ncons,neqs,1.d0,jac,ncons,err,1,1.d0,jte,1)
    return
end subroutine calc_jte
!----------------------------------------------------------------------------------
subroutine calc_jtjmui
      implicit none
      integer :: i
      jtjmui=jtj
      forall(i=1:ncons) jtjmui(i,i)=jtjmui(i,i)+mu+lambda
      return
end subroutine calc_jtjmui
!----------------------------------------------------------------------------------
! Get the Dij at a point of degeneracy by solving linear equations
! The subroutine suppose that the d^I.H.d^J is already correctly stored in fitE
! and the fit gradients in fitG.  
! WIJ and VIJ for the point for the coefficient also needs to be provided.
SUBROUTINE getDegDij(pt,WIJ,VIJ,DIJ)
  use hddata, only:  nstates
  IMPLICIT NONE
  INTEGER,intent(IN) :: pt
  DOUBLE PRECISION,intent(IN) :: WIJ(nstates,nstates),VIJ(nvibs,nstates,nstates)
  DOUBLE PRECISION,intent(OUT) :: DIJ(nstates,nstates)
  integer :: nDij
  double precision :: EMat(nstates*(nstates-1)/2,nstates*(nstates-1)/2)
  double precision,dimension(nstates*(nstates-1)/2) :: RVec
  integer,dimension(nstates*(nstates-1)/2) :: jpvt
  integer :: i,j,k,l,lb,ub
  !specifies which degeneracy group a state belongs to
  integer :: grpind(nstates)
  !mapping between dij and state indices 
  integer :: smap(nstates,nstates),dmap(nstates*(nstates-1)/2,2)
  integer :: sgn(nstates,nstates)
  double precision,dimension (nvibs) :: gvec, hvec
  double precision,dimension(:),allocatable :: WORK
  integer :: INFO,LWORK,rank

  nDij = nstates*(nstates-1)/2
  LWORK = 1000+nDij*nDij*100
  allocate(WORK(LWORK))

  grpind = 0
  !construct state groupings index table
  do i=1,dispgeoms(pt)%ndeggrpHd
    lb = dispgeoms(pt)%deg_groupsHd(i,1)
    ub = dispgeoms(pt)%deg_groupsHd(i,2)+lb-1
    grpind(lb:ub) = i
  end do!i=1,dispgeoms(pt)%ndeggrp

  !construct maps between dij and state indices
  smap = 0
  sgn  = 0
  k    = 0
  do i=1,nstates
    do j=i+1,nstates
      k=k+1
      dmap(k,1) = i
      dmap(k,2) = j
      sgn(i,j)  = 1
      smap(i,j) = k
      sgn(j,i)  = -1
      smap(j,i) = k
    end do!j
  end do!i
  if(k .ne. nDij) stop 'k .ne. nDij in getDegDij!'

  !construct linear equations for Dij
  EMat = 0.d0
  do k=1,nDij
    i=dmap(k,1)
    j=dmap(k,2)
    if(grpind(i)==grpind(j).and.grpind(i).ne.0) then
      hvec = fitG(pt,:,I,J)
      gvec = fitG(pt,:,J,J)-fitG(pt,:,I,I)
      rvec(k) = dot_product(VIJ(:,I,J),gvec) + &
                dot_product(VIJ(:,J,J)-VIJ(:,I,I),hvec)
      do l=1,nstates
        if(grpind(i).ne.grpind(l)) cycle
        if(l/=i) EMat(k,smap(l,i)) = EMat(k,smap(l,i)) + sgn(l,i)*&
          (-dot_product(fitG(pt,:,l,j),gvec)+2.d0*dot_product(fitG(pt,:,l,i),hvec))
        if(l/=j) EMat(k,smap(l,j)) = EMat(k,smap(l,j)) + sgn(l,j)*&
          (-dot_product(fitG(pt,:,l,i),gvec)-2.d0*dot_product(fitG(pt,:,l,j),hvec))
      end do!l
    else!i and j are not in same degeneracy group
      rvec(k) = -WIJ(i,j)
      do l=1,nstates
         if(grpind(l)==grpind(j))&
            EMat(k,smap(l,i)) = EMat(k,smap(l,i)) + fitE(pt,l,j)*sgn(l,i)
         if(grpind(l)==grpind(i))&
            EMat(k,smap(l,j)) = EMat(k,smap(l,j)) + fitE(pt,l,i)*sgn(l,j)
      end do
    end if!i and j are in same degeneracy group
  end do!k

  !solve the Dij values
  jpvt=0
  CALL DGELSY(nDij,nDij,1,EMat,nDij,rvec,nDij,jpvt,1.d-9,rank,WORK,LWORK,INFO)
  !if(rank<nDij)print "(A,I4,A,I5)","Reduced rank at point ",pt,", rank=",rank
  if(INFO/=0) then
     print *,"Failed to solve linear equations in getDegDij"
     print *,"Point index:",pt
     print *,"ckl at point :"
     do i=1,nstates
       print "(10F10.5)",ckl(pt,:,i)
     end do
     print *,"Block map: "
     print "(10(2I4,4x))",(dmap(i,:),i=1,nDij)
     print *,"EMat  :"
     do i=1,nstates*(nstates-1)/2
       print "(20F10.5)",EMat(:,i)
     end do
     print *,"INFO=",info
     print "(A,20I4)","Group indexing for states:  ",grpind
     print *,"fitE at point"
     do i=1,nstates
       print "(20F10.5)",fitE(pt,i,:)
     end do
     print *,"fitG at point"
     do i=1,nstates
       do j=i,nstates
         print "(A,2I4)","Block",i,j
         print "(20F10.5)",fitG(pt,:,i,j)
       end do
     end do
     stop
  end if

  !pack the solution into Dij matrix
  DIJ = 0.d0
  do k=1,nDij
    i=dmap(k,1)
    j=dmap(k,2)
    DIJ(i,j) = rvec(k)
    DIJ(j,i) = -rvec(k)
  end do!k 

  return
END SUBROUTINE getDegDij
!----------------------------------------------------------------------------------
!calculate error of fit
subroutine getError
  use progdata, only: printlvl,ev2cm
  use hddata, only: nstates
  implicit none
  integer :: j,k,l,nvpt,inc_e,inc_cp,inc_grad
  real*8 :: nrmener,avgener,nrmgrad,avggrad
  real*8 :: gnrm,dgrd,dcp,nrmdcp,nrmcp,ncp,de1,de2
  real*8, external :: dnrm2
  integer :: maxe_pt ! Max E error point index
  real*8  :: maxe    ! Max E error

  maxe_pt = 1
  maxe    = 0.d0
  
  nrmener = 0.d0
  nrmgrad = 0.d0
  avgener = 0.d0
  avggrad = 0.d0
  nrmdcp = 0.d0
  nrmcp  = 0.d0
  inc_grad = 0
  inc_e   = 0
  inc_cp  = 0

  do j=1,npoints
    nvpt=dispgeoms(j)%nvibs
    if(ptWeights(j) .gt. ptwthreshold) then
      do k=1,nstates

        do l=k+1,nstates
          if(hasGrad(j,k,l).and.(incgrad(j,k,l).or.incgrad(j,l,k))) then
            inc_cp=inc_cp+1
            de1 = max(abs(dispgeoms(j)%energy(k,k)-dispgeoms(j)%energy(l,l)),1.d-5)
            de2 = max(abs(fitE(j,k,k)-fitE(j,l,l)),1.d-5)
            dcp = dot_product(&
                  dispgeoms(j)%grads(:nvpt,k,l)/de1-fitG(j,:nvpt,k,l)/de2,&
                  dispgeoms(j)%grads(:nvpt,k,l)/de1-fitG(j,:nvpt,k,l)/de2)
            ncp = dot_product(dispgeoms(j)%grads(:nvpt,k,l),&
                              dispgeoms(j)%grads(:nvpt,k,l))/de1**2
            nrmdcp = nrmdcp+dcp/ncp

            dcp = dot_product(dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l),&
                              dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l))
            ncp = dot_product(dispgeoms(j)%grads(:nvpt,k,l),&
                              dispgeoms(j)%grads(:nvpt,k,l))
            nrmcp = nrmcp+dcp/ncp
            if(dcp>4d-2*ncp .and. printlvl>2) print "(4x,A,I5,A,2I2,A,F9.2)",&
              "Large coupling error at pt",j," bkl",k,l,": ",sqrt(dcp/ncp)
          end if !coupling included
        end do

        if(dispgeoms(j)%energy(k,k) .lt. min(energyT_en(1),energyT(1))) then
          if(hasEner(j,k)) then
            nrmener = nrmener + (dispgeoms(j)%energy(k,k)-fitE(j,k,k))**2
            avgener = avgener + abs(dispgeoms(j)%energy(k,k)-fitE(j,k,k))
            inc_e=inc_e+1
          end if
          if(hasGrad(j,k,k)) then
            dgrd = dot_product(dispgeoms(j)%grads(:nvpt,k,k)-fitG(j,:nvpt,k,k),&
                               dispgeoms(j)%grads(:nvpt,k,k)-fitG(j,:nvpt,k,k))
            gnrm = dot_product(dispgeoms(j)%grads(:nvpt,k,k),&
                               dispgeoms(j)%grads(:nvpt,k,k))
            dgrd = dgrd/gnrm
            if(gnrm.gt.1.d0) then
              nrmgrad = nrmgrad + dgrd
              avggrad = avggrad + sqrt(dgrd)
              inc_grad = inc_grad + 1
            else
              if(printlvl>1) print "(5x,A,I5,A,E14.4,A,E14.4)",&
                                   "Small gradients excluded at point ",j,&
                                   ".  Norm of grad =",sqrt(gnrm),&
                                   ",Norm of error = ",sqrt(dgrd*gnrm)
            end if !(gnrm>1D-4)
          end if!hasGrad
        end if!dispgeoms(j)%energy(k,k)

      end do!k
    end if!ptWeights(j) .gt. ptwthreshold
  end do!j

  if(printlvl>1)&
    print *,"    ",inc_grad," point/states included in gradient RMS analysis"
  if(printlvl>1)&
    print *,"    ",inc_e  ," point/states included in enery RMS analysis"
  if(printlvl>1)&
    print *,"    ",inc_cp  ," point/states included in coupling RMS analysis"

  nrmener = sqrt(nrmener/inc_e)
  nrmgrad = sqrt(nrmgrad/inc_grad)
  avgener = avgener/inc_e
  avggrad = avggrad/inc_grad

  rmsee=nrmener*ev2cm
  mue=avgener*ev2cm
  rmseg=nrmgrad
  rmsec=sqrt(nrmcp/inc_cp)

  return
end
!----------------------------------------------------------------------------------
!*
! printerrorE: print energy error for each point at the end of a calculation.
!*
subroutine printerrorE
  use progdata, only: ev2cm
  use hddata,   only: nstates
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

end module makesurfdata
!----------------------------------------------------------------------------------
subroutine update(info)
      use makesurfdata, only:ncons,coef,jtjmui,jte,BL,useGPU
      implicit none
      integer, intent(out) :: info

      BL=jte
      if(useGPU) then
        stop "*** No GPU on NERSC! Please set useGPU=.false. ***"
        !call cuda_cholesky(ncons,jtjmui,BL,info)
      else
            call cholesky(ncons,jtjmui,BL,info)
      end if

      if(info .ne. 0) then
            print *, 'Error in Cholesky decomposition!', info
            return
      end if

      coef=coef-BL

      return
end subroutine update
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
!determined the phases of Wavefunctions to best reproduce ab initio couplings
subroutine fixphase(nvibs,fitgrad,abgrad,ckl,phaseList,incGrad,cpw)
  use hddata, only: nstates
  implicit none
  integer, intent(in) :: nvibs
  real*8,dimension(nvibs,nstates,nstates),intent(inout) :: fitgrad
  real*8,dimension(nvibs,nstates,nstates),intent(in) :: abgrad
  real*8,dimension(nstates,nstates),intent(inout) :: ckl
  integer,dimension(2**(nstates-1),nstates),intent(in) :: phaseList
  LOGICAL,dimension(nstates,nstates),intent(in) :: incGrad
  DOUBLE PRECISION,dimension(nstates,nstates),INTENT(IN) :: cpw
  integer :: j,k,l,m,minID
  integer, allocatable :: phaseMat(:,:)
  real*8, allocatable :: diff(:,:)
  real*8 :: errNorm,minerr
  LOGICAL,dimension(nstates,nstates) :: hasGrad

  allocate(phaseMat(nstates,nstates),diff(nstates,nstates))
  hasGrad=incGrad.or.transpose(incGrad)

  minID=0
  do j=1,2**(nstates-1)
    do k=1,nstates
      phaseMat(k,1:nstates)=phaseList(j,k)*phaseList(j,1:nstates)
    end do
    errNorm=0.d0
    do k=1,nvibs
      diff=phaseMat*fitgrad(k,:,:)-abgrad(k,:,:)
      do l=2,nstates
        do m=1,l-1
          !errNorm = errNorm + diff(l,m)**2
          if(hasGrad(l,m))  errNorm = errNorm + (cpw(l,m)*diff(l,m))**2
        end do!m
      end do!l
    end do !k=1,nvibs
    if(minID.eq.0 .or. errNorm.lt.minerr) then
      minID=j
      minerr=errNorm
    end if!(minID==0.or.errNorm<minerr)
  end do !j=1,2**nstates

  do k=1,nstates
    phaseMat(k,1:nstates)=phaseList(minID,k)*phaseList(minID,1:nstates)
    ckl(:,k)=phaseList(minID,k)*ckl(:,k)
  end do
  do k=1,nvibs
    fitgrad(k,:,:)=phaseMat*fitgrad(k,:,:)
  end do!k=1,nvibs

  return
end SUBROUTINE fixphase
!----------------------------------------------------------------------------------
!Read in geometry and ab initio data for a set of points which will be used
!in the Hd fitting procedure.
SUBROUTINE readdisps()
  use hddata
  use progdata,only:natoms,printlvl,au2cm,au2ev
  use makesurfdata
  IMPLICIT NONE
  INTEGER :: i,j,k,l,ios,fid
  DOUBLE PRECISION, allocatable, dimension(:,:) :: cgrads,cgeoms
  DOUBLE PRECISION, allocatable, dimension(:,:) :: eners
  character(255), external :: filename
  character(255) :: infile
  character(3), dimension(natoms) :: atoms
  double precision, dimension(natoms) :: anums,masses
  !ptinfile : number of points in a geometry input file.
  !lb, ub   : energies for states between lb and ub are given in energy file
  integer :: ptinfile,lb,ub
  double precision :: NaN
  double precision,external :: dnrm2

  allocate(cgrads(3*natoms,npoints),cgeoms(3*natoms,npoints))
  allocate(eners(nstates,npoints))
  NaN = 0d0
  NaN = NaN/NaN

  if(printlvl>0) print '(4X,A,I8,A)','Reading',npoints,' displacements'
  if(allocated(dispgeoms)) deallocate(dispgeoms)
  allocate(dispgeoms(npoints))
  do j = 1,npoints
     allocate(dispgeoms(j)%cgeom(3*natoms))
     allocate(dispgeoms(j)%igeom(ncoord))
     allocate(dispgeoms(j)%energy(nstates,nstates))
     allocate(dispgeoms(j)%grads(3*natoms,nstates,nstates))
     dispgeoms(j)%grads=NaN
     allocate(dispgeoms(j)%bmat(ncoord,3*natoms))
     allocate(dispgeoms(j)%cbmat(ncoord,3*natoms))
     allocate(dispgeoms(j)%lmat(3*natoms,3*natoms))
     allocate(dispgeoms(j)%scale(3*natoms))
     allocate(dispgeoms(j)%eval(3*natoms))
  end do
  allocate(hasEner(npoints,nstates))
  allocate(hasGrad(npoints,nstates,nstates))

  hasEner = .false.
  hasGrad = .false.

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
      call genEnerGroups(dispgeoms(j),deg_cap)
      hasEner(j,lb:ub)=.true.
  enddo

  !read in gradients and couplings data
  if(printlvl>1) print 1001,'- Reading gradients and couplings'
  do j=lb,ub
    do k=j,ub
      infile = filename(j,k,grdfptn,cpfptn)
      if(printlvl>3) print 1000,"Searching for gradients in <"//trim(infile)//">"
      ptinfile=npoints
      call readGrads(infile,ptinfile,natoms,cgrads)
      if(npoints .ne. ptinfile) then
            write(*,"('npoints .ne. ptinfile in readGrads: ', 2i3)") j,k
            stop
      end if
      print 1002,"found ",ptinfile, " gradient data entries." 
      do l=1,ptinfile
        if(j/=k) then
          !interstate couplings
          !dispgeoms(l)%grads(:,j,k)=cgrads(:,l)*(eners(k,l)-eners(j,l))
          ! CLM: input is NAD in au/bohr, so we convert E(au) to E(eV)
          dispgeoms(l)%grads(:,j,k)=cgrads(:,l)*au2ev

        else
          !energy gradients
          dispgeoms(l)%grads(:,j,k)=cgrads(:,l)*au2ev
        end if
        if(j/=k) dispgeoms(l)%grads(:,k,j)=dispgeoms(l)%grads(:,j,k)
        hasGrad(l,j,k)=.true.
        hasGrad(l,k,j)=.true.
      enddo!l
    enddo!k
  enddo!j

  !process the data that have been read in
  if(printlvl>0) print 1000,"Generating displacement Wilson's B-Matrices"
  if(printlvl>0) &
      print 1000,'Local internal coordinates will be constructed from B-matrices'
  do l = 1,npoints

    if(printlvl>0) print *,""
    if(printlvl>0) print *,"      Constructing Wilson B matrix at point ",l," :"
    call buildWBMat(natoms,ncoord,dispgeoms(l)%cgeom,dispgeoms(l)%igeom,&
                                  dispgeoms(l)%cbmat)
    if(printlvl>0)then
        print *,"  Internal Coordinates"
        print "(10E11.3)",dispgeoms(l)%igeom
    end if

    if(printlvl>3)then
        print *,"  Wilson B Matrix in nascent coordinate system"
        do j=1,ncoord
            print "(10E11.3)",dispgeoms(l)%cbmat(j,:)
        end do
    end if

    if(printlvl>1)then
        print *,"  Contribution of each nascent coordinate to B matrix"
        print "(10E11.3)",(dnrm2(3*natoms,dispgeoms(l)%cbmat(j,1),ncoord),j=1,ncoord)
    end if

    if(printlvl>0)print *,"      Constructing local coordinate system for point ",l
    CALL makeLocalIntCoord(dispgeoms(l),nstates,intGradT,intGradS,nvibs)

    if(printlvl>0)then
       print *,"  Wilson B Matrix in fitting coordinate system"
       do j=1,ncoord
         print "(10E11.3)",dispgeoms(l)%bmat(j,:)
       end do
    end if

    if(printlvl>1)then
       print *,"  Energy gradients in fitting coordinate system:"
       do j=1,nstates
         if(.not.hasGrad(l,j,j))cycle
         write(unit=*,fmt="(A,I3,A)",advance='no')"   state ",j," : "
         print "(10E11.3)",dispgeoms(l)%grads(:,j,j)
       end do
       print *,"  Derivative couplings in fitting coordinate system:"
       do j=1,nstates-1
         do k=j+1,nstates
           if(.not.hasGrad(l,j,k))cycle
           write(unit=*,fmt="(A,I2,A,I2,A)",advance='no')"   block(",j,",",k,") : "
           print "(10E11.3)",dispgeoms(l)%grads(:,j,k)
         end do!k
       end do!j
    end if!printlvl>1

  end do!l

  print *, "Calling readDispOptions.."
  call readDispOptions(exclEner,exclGrad)
  
  
  call FLUnit(fid)
  open(fid,file='grd.out')
  do i=1,npoints
    do j=1,nstates
      l=dispgeoms(i)%nvibs
      write(fid,"(2i4,e15.6)") i,j,&
      dot_product(dispgeoms(i)%grads(1:l,j,j),dispgeoms(i)%grads(1:l,j,j))
    end do
  end do
  close(fid)

  call FLUnit(fid)
  open(fid,file='coup.out')
  do i=1,npoints
    l=dispgeoms(i)%nvibs
    do j=1,nstates-1
      do k=j+1,nstates
        write(fid,"(3i4,e15.6)") i,j,k,&
        dot_product(dispgeoms(i)%grads(1:l,j,k),dispgeoms(i)%grads(1:l,j,k))
      end do
    end do
  end do
  close(fid)

1000 format(7X,A)
1001 format(4X,A)
1002 format(10X,A,I10,A)

  return
end SUBROUTINE readdisps
!==================================================================================
!load pointwise options
subroutine readDispOptions(exclEner,exclGrad)
  use hddata, only:nstates
  use progdata, only: printlvl
  use makesurfdata, only: TEqList,npoints,ptWeights
  implicit none
  TYPE(TEqList) :: exclEner,exclGrad
  character(255) :: comment
  integer :: PTFL,ios
  integer :: ptid,i,j,tmp,lsID
  character(2) :: job
  integer, allocatable :: tpLists(:,:,:) 
  integer :: tpCount(2)
  real*8 :: wt

  allocate(tpLists(2,npoints*nstates*(nstates+1)/2,3)) !EE, EG
  if(allocated(ptWeights)) deallocate(ptWeights)
  allocate(ptWeights(npoints))
  ptWeights=1.d0

  if(printlvl>0) print *,"   Reading point specific options."
  call FLUnit(PTFL)
  open(unit=PTFL,file='points.in',access='sequential',form='formatted',&
       status='old',action='read',position='rewind',iostat=ios)
  if(ios .ne. 0) then
    print *,"Cannot open points.in.  Skipping..."
    exclEner%Length=0
    exclGrad%Length=0
    return
  end if !ios .ne. 0
  read(PTFL,1000,IOSTAT=ios) comment
  if (ios .ne. 0) stop "Error reading points.in file! Line 1"
  read(PTFL,1000,IOSTAT=ios) comment
  if (ios .ne. 0) stop "Error reading points.in file! Line 2"
  tpCount=0

  print *, "Reading through data in points.in."
  do while(ios==0)
    read(PTFL,1001,IOSTAT=ios) job,ptid,i,j
    if(ios/=0) then
      print *, "IOS = ", ios
      exit
    end if
    if(i>j) then
      Print *,"Warning: i>j in points.in.  i and j switched"
      tmp=i
      i=j
      j=tmp
    end if

    if(job=='WT'.or.job=='wt'.or.job=='Wt'.or.job=='wT') then
      if(ptid .lt. 0) then
        wt=-1.d0/dble(ptid)
      else
        wt=dble(ptid)
      end if
      if(printlvl>0) &
        print "(3A,2I7,A,F12.7)","    Job:",job,", Point Range",i,j,", wt=",wt
      if(i<1.or.j>npoints) then
        print *,"Point index out of range. "
        i=max(i,1)
        j=min(j,npoints)
      end if
      ptWeights(i:j)=wt
      cycle
    end if

    if(printlvl>0) print "(3A,I7,A,2I4)","    Job:",job,", Point",ptid,",states",i,j
    if(ptid<1.or.ptid>npoints) then
      print *,"Point index out of range."
      cycle
    end if
    if(i>nstates.or.j<1) then
      print *,"Warning: state id out of range in points.in. input removed"
      cycle
    end if

    select case (job)
      case ('EG','eg','eG','Eg')
        lsID=1
      case ('EE','ee','eE','Ee')
        lsID=2
      case default
        print *,"Warning: control characters not recognized in points.in"
        cycle
    end select
    tpCount(lsID)=tpCount(lsID)+1
    tpLists(lsID,tpCount(lsID),1)=ptid
    tpLists(lsID,tpCount(lsID),2)=i
    tpLists(lsID,tpCount(lsID),3)=j
  end do

  call fillout(exclGrad,1)
  call fillout(exclEner,2)
  close(unit=PTFL)

  !deallocate(tpLists)

1000 format(a255)
1001 format(a2,x,i5,i5,i5)

CONTAINS
  SUBROUTINE fillout(list,job)
    IMPLICIT NONE
    TYPE(TEqList),INTENT(OUT) :: list
    INTEGER,INTENT(IN) :: job
    CHARACTER(40),DIMENSION(2) :: titles
    integer   :: i

    titles(1)="Gradients Excluded"
    titles(2)="Energies Excluded"

    list%length=tpCount(job)
    if(allocated(list%List)) deallocate(list%List)
    allocate(list%List(tpCount(job)))
    if(tpCount(job)>0) then
      write(*,*) trim(titles(job))
      write (*,1004)
      write (*,1002)
      write (*,1004)
      do i=1,tpCount(job)
        list%List(i)%point=tpLists(job,i,1)
        list%List(i)%i    =tpLists(job,i,2)
        list%List(i)%j    =tpLists(job,i,3)
        write(*,1003) tpLists(job,i,:)
      end do
      write (*,1004)
      write (*,'(/)')
    end if
  1002 format(4X,"POINT  STATE1  STATE2")
  1003 format(4X,I4,3X,I4,4X,I4)
  1004 format(4X,"---------------------")
  END SUBROUTINE fillout
end
!==================================================================================
subroutine initmakesurf
  use progdata, only:printlvl
  use combinatorial
  use hddata, only:nstates,ANN
  use makesurfdata
  implicit none
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: tmpEqMap
  DOUBLE PRECISION,dimension(:),allocatable :: tmpW
  integer :: i,j,k,maxeqs
  integer :: ieq,g,s1,s2,pt
  real*8 :: t1,t2,t3

  call getPtList()

  !For eigenvectors
  if(allocated(ckl)) deallocate(ckl)
  allocate(ckl(npoints,nstates,nstates))
  if(allocated(fitE)) deallocate(fitE)
  if(allocated(fitG)) deallocate(fitG)
  allocate(fitE(npoints,nstates,nstates))
  allocate(fitG(npoints,nvibs,nstates,nstates))

  !Generate list of all possible phase factors for fixphase
  if(allocated(phaseList)) deallocate(phaseList)
  allocate(phaseList(2**(nstates-1),nstates))
  phaseList(1,1:nstates)=1
  do i=2,2**(nstates-1)
    phaseList(i,:)=phaseList(i-1,:)
    j=2
    do while(phaseList(i,j)==-1)
      phaseList(i,j)=1
      j=j+1
    end do
    phaseList(i,j)=-1
  end do

  !Permutation List
  if(allocated(pmtList)) deallocate(pmtList)
  allocate(pmtList(factl(nstates),nstates))
  pmtList=Permutation(nstates,factl(nstates))

  !make equations
  maxeqs=npoints*nstates*(nstates+1)/2*(nvibs+1)
  allocate(tmpEqMap(maxeqs,4))
  allocate(tmpW(maxeqs))
  call makeEqMap(maxEqs,tmpEqMap,tmpW)

  !connecting weights
  t1=-1.d0; t2=-1.d0; t3=-1.d0
  do ieq=1,neqs
    pt=tmpEqMap(ieq,1)
    s1=tmpEqMap(ieq,2)
    s2=tmpEqMap(ieq,3)
    g =tmpEqMap(ieq,4)
    if(g.eq.0) then !energy
      t1=max(t1,(dispgeoms(pt)%energy(s1,s2))**2)
    else 
      if(s1.eq.s2) then !gradient
        t2=max(t2,(dispgeoms(pt)%grads(g,s1,s2))**2)
      else !interstate coupling
        t3=max(t3,(dispgeoms(pt)%grads(g,s1,s2))**2)
      end if
    end if
  end do
  rho_energy=1.d0
  rho_grad=sqrt(t1/t2)
  rho_fij=sqrt(t1/t3)

  write(*,"('Connecting weights for gradients and interstate couplings: ')")
  write(*,"(2e15.6)") rho_grad,rho_fij

  !remake equations
  w_energy=w_energy*rho_energy
  w_grad=w_grad*rho_grad
  w_fij=w_fij*rho_fij
  call makeEqMap(maxEqs,tmpEqMap,tmpW)
  if(allocated(EqMap)) deallocate(EqMap)
  allocate(EqMap(neqs,4))
  if(allocated(weight)) deallocate(weight)
  allocate(weight(neqs))
  weight=tmpW(1:neqs)
  eqmap(1:neqs,:)=tmpEqMap(1:neqs,:)

  !print equations
  if(printlvl>0) print*,"    Number of Least Squares Equations:",neqs
  ncons=ANN%L
  if(neqs<ncons.and.printlvl>0) &
    print *,"    Warning: neqs<ncons, data set might be inadequate."

  print *,"Information of Equations:"
  write(*,"('       No.','   Point','   i','   j','   g','     Weight')")
  do ieq=1,neqs
    write(*,"(i10,i8,3i4,e11.3)") ieq,eqmap(ieq,:),weight(ieq)
  end do

  !For L-M
  allocate(coef(ncons),coefc(ncons))
  allocate(jtj(ncons,ncons),jtjmui(ncons,ncons),jte(ncons),BL(ncons))
  allocate(Jac(ncons,neqs),err(neqs))

  return
end
!-------------------------------------------------------------------------------
!reorder degenerate eigenvectors to best fit ab initio gradients
!both gradients and couplings are considered when doing the ordering
!-------------------------------------------------------------------------------
SUBROUTINE gradOrder(ptid,fitpt,abpt,ckl,pmtList,LDP,w_en,w_grd,w_cp)
  use combinatorial
  use hddata, only: nstates
  use progdata, only: abpoint,printlvl
  use makesurfdata, only:hasEner,hasGrad
  IMPLICIT NONE
  type(abpoint),INTENT(IN) :: abpt
  type(abpoint),INTENT(INOUT) :: fitpt
  DOUBLE PRECISION,DIMENSION(nstates,nstates),intent(INOUT) :: ckl
  DOUBLE PRECISION,INTENT(IN) :: w_en,w_grd,w_cp
  INTEGER,INTENT(IN) :: ptid,LDP
  INTEGER,DIMENSION(LDP,nstates),INTENT(IN) :: pmtList

  integer :: i, j, ldeg, udeg, ndeg, pmt, best_p, m
  double precision :: min_err, err, err1
  DOUBLE PRECISION,DIMENSION(nstates,nstates) :: ckl_new
  DOUBLE PRECISION,DIMENSION(abpt%nvibs) :: diff,diff2
  DOUBLE PRECISION,DIMENSION(abpt%nvibs,nstates,nstates) :: gradnew
  double precision,dimension(nstates,nstates) :: en_new
  double precision :: ndiff
  integer :: reord(nstates)

  if(printlvl>3.and.fitpt%ndeggrp>0) print "(A,I4)","Reordering point",ptid
  min_err = 0.d0
  err1 = 0.d0
  best_p = 0
  do i=1,fitpt%ndeggrp
    ldeg=fitpt%deg_groups(i,1)
    ndeg=fitpt%deg_groups(i,2)
    udeg=ldeg+ndeg-1
    if(printlvl>3) then
      print "(A,I3,A,2I3)","Permutation group",i,", range ",ldeg,udeg
    end if

    do pmt=1,factl(ndeg)
      err=dble(0)
      ! establish new ordering and store in array reord
      do j=1,nstates
        reord(j)=j
      end do
      do j=ldeg,udeg
        reord(j)=pmtList(pmt,j-ldeg+1)+ldeg-1
      end do
      if(printlvl>3) print "(A,10I3)","  New order: ",reord
      do j=ldeg,udeg
      ! calculate error contributions from gradients 
        if(hasGrad(ptid,j,j)) then
          diff=fitpt%grads(:abpt%nvibs,reord(j),reord(j))-abpt%grads(:abpt%nvibs,j,j)
          err=err+dot_product(diff,diff)*w_grd**2
        end if
      ! calculate error contributions from energies 
        if(hasEner(ptid,j))  &
          err=err+(fitpt%energy(reord(j),reord(j))-abpt%energy(j,j))**2*w_en**2
      ! add coupling contrubutions
        do m=1,nstates
           if(m==j) cycle
           if(hasGrad(ptid,m,j)) then
             diff = fitpt%grads(:abpt%nvibs,reord(m),reord(j))&
                   -abpt%grads(:abpt%nvibs,m,j)
             diff2 = fitpt%grads(:abpt%nvibs,reord(m),reord(j))&
                    +abpt%grads(:abpt%nvibs,m,j)
             ndiff = min(dot_product(diff,diff),dot_product(diff2,diff2))
             err=err+ndiff*w_cp**2
           end if
        end do
      end do!j=ldeg,udeg
      if(printlvl>3) print "(A,I3,A,E12.5)","pmt",pmt,",err=",sqrt(err)
      if(pmt==1)then
        min_err=err
        best_p=pmt
        err1=err
      else if (err<min_err) then
        min_err=err
        best_p=pmt
      end if
    end do!pmt=1,factl(ndeg)

    if(printlvl>3)  print "(A,I3)","best permutation:",best_p
    ckl_new=ckl
    en_new =fitpt%energy
    gradnew=fitpt%grads(:abpt%nvibs,:,:)
    do j=ldeg,udeg
      m=pmtList(best_p,j-ldeg+1)+ldeg-1
      ckl_new(:,j) = ckl(:,m)
      en_new(:,j) = fitpt%energy(:,m)
      gradnew(:,:,j) = fitpt%grads(:abpt%nvibs,:,m)
    end do
    if(printlvl>0.and.best_p/=1) then 
      print 1000,ptid,ldeg,ldeg+ndeg-1,sqrt(err1),sqrt(min_err)
    end if
    fitpt%grads(:abpt%nvibs,:,:)=gradnew
    fitpt%energy=en_new
    do j=ldeg,udeg
      m=pmtList(best_p,j-ldeg+1)+ldeg-1
      gradnew(:,j,:)=fitpt%grads(:abpt%nvibs,m,:)
      en_new(j,:)=fitpt%energy(m,:)
    end do
    ckl=ckl_new
    fitpt%grads(:abpt%nvibs,:,:) = gradnew

  end do !i=1,fitpt%ndeggrp

1000 format(6X,"Point ",I4," States ",I2," to ",I2,&
          " ordered by gradients. Error:",E11.4,"->",E11.4)
END SUBROUTINE gradOrder
!-------------------------------------------------------------------------------
subroutine copy_coef
  use makesurfdata, only: coef,coefc
  coefc=coef
  return
end
!-------------------------------------------------------------------------------
subroutine rec_coef
  use makesurfdata, only: coef,coefc
  coef=coefc
  return
end
!-------------------------------------------------------------------------------
subroutine init_coef
  use makesurfdata, only: coef,coefc
  call random_number(coef)
  coef=2.d0*coef-1.d0
  return
end
!-------------------------------------------------------------------------------
subroutine savecoef(outfl)
  use makesurfdata, only: coef,ncons
  implicit none
  character(len=*) :: outfl
  integer :: i,fid,ios

  call FLUnit(fid)
  open(unit=fid,file=trim(adjustl(outfl)),iostat=ios)
  if(ios.ne.0) stop 'Error writing parameters to file!'
  do i=1,ncons
    write(fid,*) coef(i)
  end do
  close(fid)

  return
end
!-------------------------------------------------------------------------------
subroutine readcoef(wfile)
  use makesurfdata, only: coef,ncons
  implicit none
  character(len=*) :: wfile
  integer :: i,fid,ios

  call FLUnit(fid)
  open(unit=fid,file=trim(adjustl(wfile)),access='sequential',&
  form='formatted',STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  if(ios.ne.0) stop 'Error opening parameter file!'
  do i=1,ncons
    read(fid,*,iostat=ios) coef(i)
    if(ios.ne.0) stop 'Error reading NN parameter file!'
  end do 
  close(fid)

  return
end
!-------------------------------------------------------------------------------
subroutine idx(i,j,nstates,id)
  !get array index for the matrix element (i,j)
  !in the order of upper triangle matrix
  implicit none
  integer, intent(in) :: i,j,nstates
  integer, intent(out) :: id
  integer :: ii,jj,k,l

  ii=min(i,j)
  jj=max(i,j)

  id=0
  do k=1,ii-1
    id=id+nstates-k+1
  end do
  id=id+jj-ii+1

  return
end
!-------------------------------------------------------------------------------

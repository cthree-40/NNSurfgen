!===============================================================================
subroutine calc_perf(P)
  use combinatorial
  use progdata, only: abpoint,printlvl,ev2cm
  use hddata, only: nstates,EvaluateHd,ANN
  use makesurfdata
  implicit none
  real*8, intent(out) :: P
  type(abpoint) :: fitpt
  real*8, allocatable :: hmat(:,:), dhmat(:,:,:)
  real*8, allocatable :: work(:),vtmp(:)
  !variables that specify equation
  integer :: pt, s1, s2, g, ieq
  integer :: lwork,info
  integer :: i,j,k,l

  lwork=999
  allocate(hmat(nstates,nstates),dhmat(nvibs,nstates,nstates),vtmp(nstates))
  allocate(work(lwork))
  allocate(fitpt%energy(nstates,nstates))
  allocate(fitpt%grads(nvibs,nstates,nstates))

  !neural neiwork parameters
  call ANN%wb(0,coef)

  !update eigenvectors
  do i=1,npoints
    !fitpt
    fitpt%lb = dispgeoms(i)%lb
    fitpt%ub = dispgeoms(i)%ub
    fitpt%id = i
    !eigenvector
    call EvaluateHd(dispgeoms(i),nvibs,hmat,dhmat)
    ckl(i,:,:)=hmat
    call dsyev('V','U',nstates,ckl(i,:,:),nstates,vtmp,work,lwork,info)
    if(info.ne.0) then
        print *,"Failed to solve eigenvectors at point",i
        print *,"INFO=",info
        stop 'error solving eigenvalue equations'
    end if

    !fitE
    fitE(i,:,:)=matmul(transpose(ckl(i,:,:)),matmul(hmat,ckl(i,:,:)))
    fitpt%energy=fitE(i,:,:)
    call genEnerGroups(fitpt,gorder)
    if(allocated(dispgeoms(i)%deg_groupsHd)) deallocate(dispgeoms(i)%deg_groupsHd)
    allocate(dispgeoms(i)%deg_groupsHd(fitpt%ndeggrp,2))
    dispgeoms(i)%deg_groupsHd=fitpt%deg_groups
    dispgeoms(i)%ndeggrpHd=fitpt%ndeggrp
    fitpt%nvibs=dispgeoms(i)%nvibs

    !partial diagonalization
    if(fitpt%ndeggrp .gt. 0) then
      call OrthGH_Hd(fitpt,dhmat(1:fitpt%nvibs,:,:),ckl(i,:,:),100,gcutoff)
    end if

    !gradients before fixphase
    do k=1,dispgeoms(i)%nvibs
      fitpt%grads(k,:,:)=&
      matmul(transpose(ckl(i,:,:)),matmul(dhmat(k,:,:),ckl(i,:,:)))
    end do!k=1,dispgeoms(i)%nvibs

    !fixphase
    call fixphase(dispgeoms(i)%nvibs,fitpt%grads(1:dispgeoms(i)%nvibs,:,:),&
                  dispgeoms(i)%grads(1:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),&
                  phaseList,incGrad(i,:,:),gWeight(i,:,:))
    fitpt%grads(dispgeoms(i)%nvibs+1:,:,:)=0.d0
    fitG(i,:,:,:)=fitpt%grads
    fitE(i,:,:)=matmul(transpose(ckl(i,:,:)),matmul(hmat,ckl(i,:,:)))
  end do !i

  !calculate difference vectors
  do ieq = 1, neqs
    pt=eqmap(ieq,1)
    s1=eqmap(ieq,2)
    s2=eqmap(ieq,3)
    g =eqmap(ieq,4)
    if(g==0) then !is an energy
      if(s2>0) then !its adiabatic energy or adiabatic off diagonal terms(0)
        err(ieq) = fitE(pt,s1,s2) - dispgeoms(pt)%energy(s1,s2)
      else !(s2>0)
        stop 's2.le.0 in clac_perf!'
      end if !(s2>0)
    else !(g==0) is gradient
      err(ieq) = fitG(pt,g,s1,s2) - dispgeoms(pt)%grads(g,s1,s2)
    end if !(g==0)
  end do !ieq
  err=err*weight
  P = 0.5d0*(dot_product(err,err)+lambda*dot_product(coef,coef)) !regularization

  !calculate RMS Error
  call getError

  return
end
!===============================================================================
subroutine calc_jac
  use progdata, only:natoms
  use hddata
  use makesurfdata
  implicit none
  real*8, allocatable :: dpdh(:,:,:)
  real*8, allocatable :: dpdhdx(:,:,:,:)   !w.r.t. local coordinates
  real*8, allocatable :: dpdhdr(:,:,:,:)   !w.r.t. internal coordinates
  real*8, allocatable :: dtmp(:)
  !W^IJ_i =<f^I|W_i|f^J>
  !D^IJ_i = W^IJ_i / (E^J-E^I)
  !V^IJ_i =<f^I|del W_i|f^J>
  real*8, allocatable :: WIJ(:,:,:), DIJ(:,:,:)
  real*8, allocatable :: VIJ(:,:,:,:)
  real*8, allocatable :: cc(:)
  real*8, allocatable :: xin(:)
  real*8 :: oop(noopc)
  integer :: nrn, blk
  integer :: ii, o
  integer :: i,j,k,l,m,n,R,ico,igd,ipt,nt,id
  !variables that specify equation
  integer :: ieq,pt,s1,s2,g

  nt=nstates*(nstates+1)/2
  allocate(dpdh(npoints,ncons,nt))
  allocate(dpdhdx(npoints,ncons,nvibs,nt))
  allocate(dpdhdr(ncons,ncoord,nstates,nstates))
  allocate(dtmp(3*natoms))
  allocate(WIJ(npoints,nstates,nstates))
  allocate(VIJ(npoints,nvibs,nstates,nstates))
  allocate(DIJ(npoints,nstates,nstates))
  allocate(cc(npoints))
  R=ANN%RX
  allocate(xin(R))

  do i=1,npoints
    xin=dispgeoms(i)%igeom(1:R)
    do j = 1, noopc
       oop(j)=dispgeoms(i)%igeom((ncoord-noopc)+j)
    end do
    call ANN%output(xin)
    call ANN%cal_dydn
    call ANN%cal_dydp
    call ANN%cal_dadx
    call ANN%cal_dxdydn
    call ANN%cal_dpdydx

    ! Loop over blocks
    dpdh(i,:,:) = 0.0
    nrn = 1
    blk = 1
    do ii = 1, nstates
       ! PIP on diagonal
       do j = 1, abs(hddef(ii,ii)%nc)
          dpdh(i,:,blk) = dpdh(i,:,blk) + ANN%dydp(:,nrn)
          nrn = nrn + 1
       end do
       blk = blk + 1
       do j = ii+1, nstates
          ! PIP * OOP off diagonal
          do k = 1, hddef(ii,j)%nc
             do o = 1, 4 ! order
                if (hddef(ii,j)%terms(o,k) .eq. 0) cycle
                dpdh(i,:,blk) = dpdh(i,:,blk) + &
                     ANN%dydp(:,nrn)*oop(hddef(ii,j)%terms(o,k))**o
                nrn = nrn + 1
             end do
          end do
          dpdh(i,:,blk) = dpdh(i,:,blk)
          blk = blk + 1
       end do
    end do

    !dpdh(i,:,1) = ANN%dydp(:,1)  !(1, 1)
    !dpdh(i,:,2) = ANN%dydp(:,2)*oop + ANN%dydp(:,3)*oop**2 + &   ! (1, 2)
    !     ANN%dydp(:,4)*oop**3 + ANN%dydp(:,5)*oop**4             !
    !dpdh(i,:,3) = ANN%dydp(:,6)*oop + ANN%dydp(:,7)*oop**2 + &   ! (1, 3)
    !     ANN%dydp(:,8)*oop**3 + ANN%dydp(:,9)*oop**4             !
    !dpdh(i,:,4) = ANN%dydp(:,10)  !(2,2)
    !dpdh(i,:,5) = ANN%dydp(:,11)*oop + ANN%dydp(:,12)*oop**2 + & !(2,3)
    !     ANN%dydp(:,13)*oop**3 + ANN%dydp(:,14)*oop**4 + &       !
    !     ANN%dydp(:,15)*oop    + ANN%dydp(:,16)*oop**2 + &
    !     ANN%dydp(:,17)*oop**3 + ANN%dydp(:,18)*oop**4
    !dpdh(i,:,6) = ANN%dydp(:,19)  !(3,3)

     ! Loop over blocks
     nrn = 1
     dpdhdr = 0.0d0
     do ii = 1, nstates
        ! PIP on diagonal
        do j = 1, abs(hddef(ii,ii)%nc)
           ! R part of derivative
           dpdhdr(:,1:R,ii,ii) = dpdhdr(:,1:R,ii,ii) + ANN%dpdydx(:,1:R,nrn)
           ! OOP part of derivative
           dpdhdr(:,R+1:ncoord,ii,ii) = 0.0d0
           nrn = nrn + 1
        end do
        do j = ii + 1, nstates
           ! PIP*OOP off diagonal
           do k = 1, hddef(ii,j)%nc
              do o = 1, 4 ! order
                 if (hddef(ii,j)%terms(o,k) .eq. 0) cycle
                 ! R part of derivative
                 dpdhdr(:,1:R,ii,j) = dpdhdr(:,1:R,ii,j) + &
                      ANN%dpdydx(:,1:R,nrn)*oop(hddef(ii,j)%terms(o,k))**o
                 ! OOP part of derivative
                 dpdhdr(:,(R + hddef(ii,j)%terms(o,k)),ii,j) = &
                      dpdhdr(:,(R + hddef(ii,j)%terms(o,k)),ii,j) + &
                      o*ANN%dydp(:,nrn)*oop(hddef(ii,j)%terms(o,k))**(o - 1)
                 nrn = nrn + 1
              end do
           end do
           dpdhdr(:,1:R,j,ii) = dpdhdr(:,1:R,ii,j)
        end do
     end do

     
    !(1,1)
    !dpdhdr(:,1:R,1,1)=ANN%dpdydx(:,1:R,1)
    !dpdhdr(:,ncoord,1,1)=0.d0
    !(1,2)
    !dpdhdr(:,1:R,1,2)=ANN%dpdydx(:,1:R,2)*oop + ANN%dpdydx(:,1:R,3)*oop**2 +&
    !     ANN%dpdydx(:,1:R,4)*oop**3 + ANN%dpdydx(:,1:R,5)*oop**4
    !dpdhdr(:,ncoord,1,2) = ANN%dydp(:,2) + 2.d0*oop*ANN%dydp(:,3) + &
    !     3.d0*oop**2*ANN%dydp(:,4) + 4.d0*oop**3*ANN%dydp(:,5)
    !dpdhdr(:,:,2,1)=dpdhdr(:,:,1,2)

    !(1,3)
    !dpdhdr(:,1:R,1,3)=ANN%dpdydx(:,1:R,6)*oop + ANN%dpdydx(:,1:R,7)*oop**2 +&
    !     ANN%dpdydx(:,1:R,8)*oop**3 + ANN%dpdydx(:,1:R,9)*oop**4
    !dpdhdr(:,ncoord,1,3) = ANN%dydp(:,6) + 2.d0*oop*ANN%dydp(:,7) + &
    !     3.d0*oop**2*ANN%dydp(:,8) + 4.d0*oop**3*ANN%dydp(:,9)
    !dpdhdr(:,:,3,1)=dpdhdr(:,:,1,3)

    !(2,2)
    !dpdhdr(:,1:R,2,2)=ANN%dpdydx(:,1:R,10)
    !dpdhdr(:,ncoord,2,2)=0.d0

    !(2,3)
    !dpdhdr(:,1:R,2,3)=ANN%dpdydx(:,1:R,11)*oop + ANN%dpdydx(:,1:R,12)*oop**2 +&
    !     ANN%dpdydx(:,1:R,13)*oop**3 + ANN%dpdydx(:,1:R,14)*oop**4 + &
    !     ANN%dpdydx(:,1:R,15)*oop + ANN%dpdydx(:,1:R,16)*oop**2 + &
    !     ANN%dpdydx(:,1:R,17)*oop**3 + ANN%dpdydx(:,1:R,18)*oop**4
    !dpdhdr(:,ncoord,2,3) = ANN%dydp(:,11) + 2.d0*oop*ANN%dydp(:,12) + &
    !     3.d0*oop**2*ANN%dydp(:,13) + 4.d0*oop**3*ANN%dydp(:,14) + &
    !     ANN%dydp(:,15) + 2.d0*oop*ANN%dydp(:,16) + &
    !     3.d0*oop**2*ANN%dydp(:,17) + 4.d0*oop**3*ANN%dydp(:,18)
    !dpdhdr(:,:,3,2)=dpdhdr(:,:,2,3)

    !(3,3)
    !dpdhdr(:,1:R,3,3)=ANN%dpdydx(:,1:R,19)
    !dpdhdr(:,ncoord,3,3)=0.d0

    do l=1,ncons
      id=0
      do j=1,nstates
        do k=j,nstates
          id=id+1
          call dgemv('T',ncoord,3*natoms,1.d0,dispgeoms(i)%bmat,&
                         ncoord,dpdhdr(l,:,j,k),1,0.d0,dtmp,1)
          dpdhdx(i,l,1:nvibs,id)=dtmp(1:nvibs)
        end do
      end do
      if(id.ne.nt) stop 'id.ne.nt in calc_perf0!'
    end do

  end do !i

  !WIJ,VIJ,DIJ
  do ico=1,ncons
    WIJ=0.d0
    VIJ=0.d0
    DIJ=0.d0
    !construct vectors W,V and D for the current coefficient
    do i=1,nstates
      do j=i,nstates
        do l=1,nstates
          do k=1,nstates
            call idx(k,l,nstates,id)
            !calculate the multiplication of ckl
            cc=ckl(:,k,i)*ckl(:,l,j)
            WIJ(:,i,j)=WIJ(:,i,j)+cc*dpdh(:,ico,id)
            !construct the Hamiltonian contribution of VIJ 
            do igd=1,nvibs
              VIJ(:,igd,i,j)=VIJ(:,igd,i,j)+cc*dpdhdx(:,ico,igd,id)
            end do!igd
          end do!K
        end do!L

        !construct DIJ from WIJ. non zero only on the off diagonal
        if(j.ne.i) then
          do ipt=1,npoints
            if(dispgeoms(ipt)%ndeggrpHd .gt. 0) cycle
            DIJ(ipt,i,j)=WIJ(ipt,i,j)/(fitE(ipt,j,j)-fitE(ipt,i,i))
          end do
          !fill the lower triangle of D and W
          WIJ(:,j,i)   = WIJ(:,i,j)
          DIJ(:,j,i)   = -DIJ(:,i,j)
          VIJ(:,:,j,i) = VIJ(:,:,i,j)
        end if !(J.ne.I)

      end do!J
    end do!I

    !construct DIJ for points that have degeneracies
    do ipt=1,npoints
      if(dispgeoms(ipt)%ndeggrpHd .eq. 0) cycle
      call getDegDij(ipt,WIJ(ipt,:,:),VIJ(ipt,:,:,:),DIJ(ipt,:,:))
    end do !ipt

    !update the wave function dependent part of Wij and Vij
    do i=1,nstates
      do j=i,nstates
        do k=1,nstates
          if(k.ne.i) WIJ(:,i,j)=WIJ(:,i,j)+DIJ(:,k,i)*fitE(:,k,j)
          if(k.ne.j) WIJ(:,i,j)=WIJ(:,i,j)+DIJ(:,k,j)*fitE(:,k,i)
        end do
        if(i.ne.j) WIJ(:,j,i)=WIJ(:,i,j)
        !calculate the second piece of VIJ with DIJ and hIJ(=fitG)
        do igd = 1,nvibs
          do k = 1,nstates
            if(k.ne.i) VIJ(:,igd,i,j) = VIJ(:,igd,i,j) + DIJ(:,k,i)*fitG(:,igd,k,j)
            if(k.ne.j) VIJ(:,igd,i,j) = VIJ(:,igd,i,j) + DIJ(:,k,j)*fitG(:,igd,k,i)
          end do ! K = 1,nstates
          !complete lower triangle
          if(i.ne.j)  VIJ(:,igd,j,i) = VIJ(:,igd,i,j)
        end do!igd=1,nvibs 
      end do !J=I,nstates
    end do !I=1,nstates

    !calculate derivative vector with respect to the current variable
    do ieq = 1, neqs
      pt=eqmap(ieq,1)
      s1=eqmap(ieq,2)
      s2=eqmap(ieq,3)
      g =eqmap(ieq,4)
      if(g==0) then !is an energy
        if(s2>0) then !its adiabatic energy or adiabatic off diagonal terms(0)
          Jac(ico,ieq) = WIJ(pt,s1,s2)
        else !(s2>0)
          stop 's2>0 in clac_jac0!'
        end if !(s2>0)
      else !(g==0) is gradient
        Jac(ico,ieq) = VIJ(pt,g,s1,s2)
      end if!(g==0)
      Jac(ico,ieq)=Jac(ico,ieq)*weight(ieq)
    end do !ieq

  end do!ico

  return
end
!===============================================================================

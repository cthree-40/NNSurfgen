!==============================================================================
subroutine makesurf
  use progdata, only:printlvl
  use hddata
  use makesurfdata
  implicit none
  character(len=99) :: outfl
  real*8 :: P,P0
  integer :: ep,i,j,k,info,nt
  integer :: ifits
  
  call initmakesurf
  write(*,'("Training Neural Network with Levenberg-Marquardt Method")')
  i=ANN%M;j=ANN%L
  write(*,'("NN structure is ",<i+1>i6," with ",i6," parameters")') ANN%S(0:i),j

  nt=nstates*(nstates+1)/2
  nse=0
  do k=1,npoints
    !energies
    nse=nse+nstates
    !Hd
    nse=nse+nt
  end do
  if(allocated(Jac)) deallocate(Jac)
  if(allocated(err)) deallocate(err)
  allocate(Jac(npara,nse),err(nse))
  if(allocated(rmse)) deallocate(rmse)
  if(allocated(rmsh)) deallocate(rmsh)
  allocate(rmse(nstates),rmsh(nt))
  print *, "trainings = ", trainings
  do i=1,nc
    write(*,'(" Fitting Count = ",i3.1)') i
    !output
    write(outfl,"('WB',i2.2,'.txt')") i

    if (.not. restart) then
      call init_wb
    else
      write(*,'(" Restarting from WB-restart.txt.")')
      call readnet('WB-restart.txt')
    end if
    call calc_perf(P0)
    mu=1.d0

    if (trainings(i) .eq. 0) cycle
    
    do ep=1,epmax
      call copy_wb
      call calc_jac
      call calc_jtj
      call calc_jte
100   call calc_jtjmui
      call rec_wb
      call update_wb(info)
      if(info .ne. 0) then
        call rec_wb
        print*, 'Cholesky Decomposition Failed!'
        exit
      end if
      call calc_perf(P)

      if(P .lt. P0) then
        P0=P
        mu=mu*0.5d0
      else
        mu=mu*2.d0
        if(mu .gt. 1.d10) then
          write(*,"('Large mu=',1x,e10.2,1x,'exit loop.' )") mu
          exit
        end if
        goto 100
      end if

      write(*,"('Epoch: ',i4,'/',i4,2x,'P=',es9.2,2x,'rmse=',<nstates>es9.2,2x,&
                'rmsh=',<nt>es9.2,2x,'mu=',es7.0)") ep,epmax,P,rmse(:),rmsh(:),mu
      call savenet(outfl)
    end do
    write(*,"('Fitting No. ',i2.2,2x,'rmse=',<nstates>es9.2,2x,&
              'rmsh=',<nt>es9.2,2x)") i,rmse(:),rmsh(:)
    call savenet(outfl)
    ! Print error
    print *, "Printing error..."
    call printerrorE()

  end do

  return
end
!==============================================================================

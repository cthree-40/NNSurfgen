!==============================================================================
subroutine makesurf
  use progdata, only:printlvl
  use hddata
  use makesurfdata
  implicit none
  character(len=99) :: outfl
  real*8 :: P,P0
  integer :: ep,i,j,k,info
  logical :: here

  call initmakesurf
  write(*,'("Training Neural Network with Levenberg-Marquardt Method")')
  i=ANN%M;j=ANN%L
  write(*,'("NN structure is ",<i+1>i6," with ",i6," parameters")') ANN%S(0:i),j
  stop

  !first stage, without partial diagonalization
  inquire(file='NN.txt',exist=here)
  if(here) then
    call readcoef('NN.txt')
  else
    call init_coef
  end if

  !NN parameter output
  write(outfl,"('WB.txt')")

  mu=1.d0
  call calc_perf0(P0)
  do ep=1,epmax
    call copy_coef
    call calc_jac0
    call calc_jtj
    call calc_jte
100 call calc_jtjmui
    call rec_coef
    call update(info)
    if(info .ne. 0) then
      call rec_coef
      print*, 'Cholesky Decomposition Failed!'
      exit
    end if
    call calc_perf0(P)

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

    call savecoef(outfl)
    write(*,"('Epoch: ',i4,'/',i4,2x,'P=',es9.2,2x,'RMSEE=',es9.2,2x,&
              'MUEE=',es9.2,2x,'RMSEG=',es9.2,2x,'RMSEC=',es9.2,2x,'mu=',es7.0)") &
              ep,epmax,P,rmsee,mue,rmseg,rmsec,mu
  end do

  return
end
!==============================================================================

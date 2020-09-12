!  _   _  _____  _   _  ____      _     _     
! | \ | || ____|| | | ||  _ \    / \   | |    
! |  \| ||  _|  | | | || |_) |  / _ \  | |    
! | |\  || |___ | |_| ||  _ <  / ___ \ | |___ 
! |_| \_||_____| \___/ |_| \_\/_/   \_\|_____|
!  _   _  _____  _____ __        __ ___   ____   _  __ ____  
! | \ | || ____||_   _|\ \      / // _ \ |  _ \ | |/ // ___| 
! |  \| ||  _|    | |   \ \ /\ / /| | | || |_) || ' / \___ \ 
! | |\  || |___   | |    \ V  V / | |_| ||  _ < | . \  ___) |
! |_| \_||_____|  |_|     \_/\_/   \___/ |_| \_\|_|\_\|____/ 
!                         
! Feedforward Neural Network
! Fortran Module
! Structure: R-S1-S2-...-Sm
! See Hagan et al Neural Network Design and PhD thesis of Yafu Guan or
! Please cite the PCCP paper DOI:10.1039/c8cp06598e
! dydp: compute partial derivatives of output with respect to parameters
! dadx: compute gradients with respect to inputs
! dpdydx: compute partial derivatives of gradients with respect to parameters
! Yafu Guan, May 2nd, 2018
!===================================================================================
module nn_class
  implicit none
  type, public :: nn
        integer :: RX, M, L, maxS, RY
        integer, allocatable :: S(:)
        !weights and biases
        real*8, allocatable :: w(:,:,:)
        real*8, allocatable :: b(:,:)
        !weights and biases copy
        real*8, allocatable :: wc(:,:,:)
        real*8, allocatable :: bc(:,:)
        !input & output
        real*8, allocatable :: y(:), dydx(:,:)
        real*8, allocatable :: a(:,:), n(:,:)
        !back propagation
        real*8, allocatable :: dydn(:,:,:), dydp(:,:)
        !gradient
        real*8, allocatable :: dndx(:,:,:), dadx(:,:,:), dxdydn(:,:,:,:)
        real*8, allocatable :: dpdydx(:,:,:)
        !transfer
        character(len=99), dimension(:),  allocatable ::  activation
        character(len=99) :: structfl

        contains
        procedure, public :: init => nn_init
        procedure, public :: output => nn_output
        procedure, public :: wb_init => init_wb
        procedure, public :: wb_copy => copy_wb
        procedure, public :: wb_rec => rec_wb
        procedure, public :: wb => wb_inout
        procedure, public :: wb_update => update_wb
        procedure, public :: cal_dydn => calc_dydn
        procedure, public :: cal_dydp => calc_dydp
        procedure, public :: cal_dadx => calc_dadx
        procedure, public :: cal_dxdydn => calc_dxdydn
        procedure, public :: cal_dpdydx => calc_dpdydx
        procedure, public :: savenn => savepara
        procedure, public :: clean => clean_nn
  end type nn

  private :: nn_init, nn_output, calc_dydn, calc_dydp, init_wb, copy_wb, rec_wb,&
             update_wb, calc_dadx, calc_dxdydn, calc_dpdydx, wb_inout, savepara,&
             clean_nn

  contains
  !----------------------------------------------------------------------------
  subroutine nn_init(this)
        implicit none
        class(nn) :: this
        integer :: R,M,maxS,SM,L
        integer :: i,ios,fid

        call FLUnit(fid)
        open(unit=fid,file=trim(adjustl(this%structfl)),access='sequential',&
        form='formatted',STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
        if(ios.ne.0) stop 'Error opening NN structure file!'
        read(fid,*,iostat=ios) M
        if(ios.ne.0) stop 'Error reading structure file!'
        this%M=M
        allocate(this%S(0:M))
        allocate(this%activation(M))
        read(fid,*,iostat=ios) this%S(0:M)
        if(ios.ne.0) stop 'Error reading structure file!'
        this%RX=this%S(0)
        R=this%S(0)
        L=0
        do i=1,M
              L=L+this%S(i)*this%S(i-1)+this%S(i)
              read(fid,*,iostat=ios) this%activation(i)
              if(ios.ne.0) stop 'Error reading structure file!'
              call check_act(this%activation(i))
        end do
        this%L=L
        close(fid)

        maxS=maxval(this%S(0:M))
        this%maxS=maxS
        allocate(this%w(maxS,maxS,M),this%b(maxS,M))
        allocate(this%wc(maxS,maxS,M),this%bc(maxS,M))
        allocate(this%a(maxS,0:M),this%n(maxS,M))
        SM=this%S(M)
        this%RY=SM
        allocate(this%dydn(maxS,SM,M))
        allocate(this%dydp(L,SM))
        allocate(this%dndx(R,maxS,M),this%dydx(R,SM))
        allocate(this%dadx(R,maxS,0:M))
        allocate(this%dxdydn(R,maxS,SM,M))
        allocate(this%dpdydx(L,R,SM))
        allocate(this%y(SM))

        return
  end subroutine nn_init
  !-----------------------------------------------------------------------------
  subroutine clean_nn(this)
        implicit none
        class(nn) :: this
        if(allocated(this%S)) deallocate(this%S)
        if(allocated(this%w)) deallocate(this%w)
        if(allocated(this%b)) deallocate(this%b)
        if(allocated(this%wc)) deallocate(this%wc)
        if(allocated(this%bc)) deallocate(this%bc)
        if(allocated(this%y)) deallocate(this%y)
        if(allocated(this%dydx)) deallocate(this%dydx)
        if(allocated(this%a)) deallocate(this%a)
        if(allocated(this%n)) deallocate(this%n)
        if(allocated(this%dydn)) deallocate(this%dydn)
        if(allocated(this%dydp)) deallocate(this%dydp)
        if(allocated(this%dndx)) deallocate(this%dndx)
        if(allocated(this%dadx)) deallocate(this%dadx)
        if(allocated(this%dxdydn)) deallocate(this%dxdydn)
        if(allocated(this%dpdydx)) deallocate(this%dpdydx)
        if(allocated(this%activation)) deallocate(this%activation)
        return
  end subroutine clean_nn
  !-----------------------------------------------------------------------------
  subroutine init_wb(this,x1,x2)
        ![x1,x2]
        implicit none
        class(nn) :: this
        real*8, intent(in) :: x1,x2
        call random_number(this%w)
        call random_number(this%b)
        this%w=this%w*(x2-x1)+x1
        this%b=this%b*(x2-x1)+x1
        return
  end subroutine init_wb
  !----------------------------------------------------------------------------
  subroutine copy_wb(this)
        implicit none
        class(nn) :: this
        this%wc=this%w
        this%bc=this%b
        return
  end subroutine copy_wb
  !---------------------------------------------------------------------------
  subroutine rec_wb(this)
        implicit none
        class(nn) :: this
        this%w=this%wc
        this%b=this%bc
        return
  end subroutine rec_wb
  !---------------------------------------------------------------------------
  subroutine nn_output(this,x)
        implicit none
        class(nn) :: this
        real*8, intent(in) :: x(this%RX)
        integer :: i,R,M,SM

        R=this%RX
        M=this%M
        SM=this%RY
        !input
        this%a(1:R,0)=x(1:R)
        do i=1,M
              this%n(1:this%S(i),i)=this%b(1:this%S(i),i)
              call dgemv('N',this%S(i),this%S(i-1),1.d0,&
                         this%w(1:this%S(i),1:this%S(i-1),i),this%S(i),&
                         this%a(1:this%S(i-1),i-1),1,1.d0,this%n(1:this%S(i),i),1)
              call trans(this%S(i),this%n(1:this%S(i),i),this%a(1:this%S(i),i),&
                         this%activation(i))
        end do
        !output
        this%y(1:SM)=this%a(1:SM,M)

        return
  end subroutine nn_output
  !----------------------------------------------------------------------------
  subroutine calc_dydn(this)
        !depend on output
        implicit none
        class(nn) :: this
        integer :: i,j,k,l,M,SM
        real*8, allocatable :: wk1(:)

        allocate(wk1(this%maxS))

        !the last layer
        M=this%M
        SM=this%RY
        this%dydn(1:SM,1:SM,M)=0.d0
        do i=1,SM
              call dtrans(1,this%n(i,M),this%dydn(i,i,M),this%activation(M))
        end do
        !iteration
        do i=M,2,-1
              k=this%S(i-1)
              l=this%S(i)
              call dgemm('T','N',k,SM,l,1.d0,this%w(1:l,1:k,i),l,&
                         this%dydn(1:l,1:SM,i),l,0.d0,this%dydn(1:k,1:SM,i-1),k)
              call dtrans(k,this%n(1:k,i-1),wk1(1:k),this%activation(i-1))
              do j=1,SM
                    this%dydn(1:k,j,i-1)=this%dydn(1:k,j,i-1)*wk1(1:k)
              end do
        end do

        !deallocate(wk1)
        return
  end subroutine calc_dydn
  !-----------------------------------------------------------------------------
  subroutine calc_dydp(this)
        !depend on output & calc_dydn
        implicit none
        class(nn) :: this
        integer :: i,k,l,p,q,M,SM

        M=this%M
        SM=this%RY
        
        do i=1,SM
              p=0
              do k=1,M
                    q=this%S(k)
                    !dadw
                    do l=1,this%S(k-1)
                          this%dydp(p+1:p+q,i)=this%dydn(1:q,i,k)*this%a(l,k-1)
                          p=p+q
                    end do
                    !dadb
                    this%dydp(p+1:p+q,i)=this%dydn(1:q,i,k)
                    p=p+q
              end do
              if(p.ne.this%L) stop 'p.ne.this%L'
        end do

        return
  end subroutine calc_dydp
  !----------------------------------------------------------------------------
  subroutine calc_dadx(this)
        !depend on output
        implicit none
        class(nn) :: this
        integer :: i,j,k,l,R,M,SM
        real*8, allocatable :: wk1(:)

        allocate(wk1(this%maxS))

        R=this%RX
        M=this%M
        SM=this%RY

        !the zeroth layer
        this%dadx(1:R,1:R,0)=0.d0
        do i=1,R
              this%dadx(i,i,0)=1.d0
        end do
        !iteration
        do i=1,M
              !dndx(m)
              l=this%S(i)
              k=this%S(i-1)
              call dgemm('N','T',R,l,k,1.d0,this%dadx(1:R,1:k,i-1),R,&
                         this%w(1:l,1:k,i),l,0.d0,this%dndx(1:R,1:l,i),R)
              !dadx(m)
              call dtrans(l,this%n(1:l,i),wk1(1:l),this%activation(i))
              do j=1,R
                    this%dadx(j,1:l,i)=this%dndx(j,1:l,i)*wk1(1:l)
              end do
        end do

        this%dydx(1:R,1:SM)=this%dadx(1:R,1:SM,M)

        !deallocate(wk1)
        return
  end subroutine calc_dadx
  !-----------------------------------------------------------------------------
  subroutine calc_dxdydn(this)
        !depend on output, calc_dadx, calc_dydn
        implicit none
        class(nn) :: this
        integer :: i,j,k,l,o,p,R,M,SM
        real*8, allocatable :: wk(:,:), wk1(:), wk2(:,:)

        allocate(wk(this%maxS,this%maxS),wk1(this%maxS))
        allocate(wk2(this%maxS,this%maxS))

        R=this%RX
        M=this%M
        SM=this%RY
        !the last layer
        this%dxdydn(1:R,1:SM,1:SM,M)=0.d0
        call d2trans(SM,this%n(1:SM,M),wk1(1:SM),this%activation(M))
        do i=1,SM
              this%dxdydn(1:R,i,i,M)=wk1(i)*this%dndx(1:R,i,M)
        end do
        !iteration
        do l=1,R
              do i=M,2,-1
                    o=this%S(i-1)
                    call dtrans(o,this%n(1:o,i-1),wk1(1:o),this%activation(i-1))
                    p=this%S(i)
                    call dgemm('T','N',o,SM,p,1.d0,this%w(1:p,1:o,i),p,&
                                this%dxdydn(l,1:p,1:SM,i),p,0.d0,wk2(1:o,1:SM),o)
                    do k=1,SM
                          wk2(1:o,k)=wk2(1:o,k)*wk1(1:o)
                    end do
                    call d2trans(o,this%n(1:o,i-1),wk1(1:o),this%activation(i-1))
                    wk1(1:o)=wk1(1:o)*this%dndx(l,1:o,i-1)
                    call dgemm('T','N',o,SM,p,1.d0,this%w(1:p,1:o,i),p,&
                                this%dydn(1:p,1:SM,i),p,0.d0,wk(1:o,1:SM),o)
                    do k=1,SM
                          wk(1:o,k)=wk(1:o,k)*wk1(1:o)
                    end do
                    this%dxdydn(l,1:o,1:SM,i-1)=wk(1:o,1:SM)+wk2(1:o,1:SM)
              end do
        end do

        !deallocate(wk,wk1,wk2)
        return
  end subroutine calc_dxdydn
  !----------------------------------------------------------------------------
  subroutine calc_dpdydx(this)
        implicit none
        class(nn) :: this
        integer :: i,j,m,p,q,R,M0,SM,maxS
        real*8, allocatable :: wk1(:,:)

        R=this%RX
        M0=this%M
        SM=this%RY
        maxS=this%maxS

        allocate(wk1(maxS,maxS))

        do i=1,R
              p=0
              do m=1,M0
                    q=this%S(m)
                    !dw
                    wk1(1:q,1:SM)=this%dxdydn(i,1:q,1:SM,m)
                    do j=1,this%S(m-1)
                      this%dpdydx(p+1:p+q,i,1:SM)=wk1(1:q,1:SM)*this%a(j,m-1)+&
                      this%dydn(1:q,1:SM,m)*this%dadx(i,j,m-1)
                      p=p+q
                    end do
                    !db
                    this%dpdydx(p+1:p+q,i,1:SM)=wk1(1:q,1:SM)
                    p=p+q
              end do
              if(p.ne.this%L) stop 'p.ne.this%L'
        end do

        !deallocate(wk1)
        return
  end subroutine calc_dpdydx
  !----------------------------------------------------------------------------
  subroutine update_wb(this,dx)
        implicit none
        class(nn) :: this
        real*8, intent(in) :: dx(this%L)
        integer :: k,l,p,q,M

        M=this%M

        p=0
        do k=1,M
              q=this%S(k)
              !w
              do l=1,this%S(k-1)
                    this%w(1:q,l,k)=this%w(1:q,l,k)-dx(p+1:p+q)
                    p=p+q
              end do
              !b
              this%b(1:q,k)=this%b(1:q,k)-dx(p+1:p+q)
              p=p+q
        end do
        if(p.ne.this%L) stop 'p.ne.this%L'

        return
  end subroutine update_wb
  !--------------------------------------------------------------------------
  subroutine wb_inout(this,id,dx)
        implicit none
        class(nn) :: this
        integer, intent(in) :: id
        real*8 :: dx(this%L)
        integer :: k,l,p,q,M

        M=this%M

        p=0
        do k=1,M
              q=this%S(k)
              !w
              do l=1,this%S(k-1)
                    if(id.eq.0) then !in
                          this%w(1:q,l,k)=dx(p+1:p+q)
                    else !out
                          dx(p+1:p+q)=this%w(1:q,l,k)
                    end if
                    p=p+q
              end do
              !b
              if(id.eq.0) then !in
                    this%b(1:q,k)=dx(p+1:p+q)
              else !out
                    dx(p+1:p+q)=this%b(1:q,k)
              end if
              p=p+q
        end do
        if(p.ne.this%L) stop 'p.ne.this%L'

        return
  end subroutine wb_inout
  !---------------------------------------------------------------------------
  subroutine savepara(this,filename,io)
    implicit none
    class(nn) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in) :: io

    integer :: i,fid,ios
    real*8, allocatable :: para(:)

    allocate(para(this%L))
    if(io.eq.0) then
      !read parameters
      call FLUnit(fid)
      open(unit=fid,file=trim(adjustl(filename)),access='sequential',&
      form='formatted',STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
      if(ios.ne.0) stop 'Error opening NN parameter file!'
      do i=1,this%L
        read(fid,*,iostat=ios) para(i)
        if(ios.ne.0) stop 'Error reading NN parameter file!'
      end do
      call this%wb(0,para)
    else
      !save parameters
      call FLUnit(fid)
      open(unit=fid,file=trim(adjustl(filename)),iostat=ios)
      if(ios.ne.0) stop 'Error writing NN parameter file!'
      call this%wb(1,para)
      do i=1,this%L
        write(fid,*) para(i)
      end do
    end if

    close(fid)
    deallocate(para)
    return
  end subroutine savepara
  !---------------------------------------------------------------------------
end
!=================================================================================
subroutine check_act(activate)
      implicit none
      character(len=*), intent(in) :: activate

      if(trim(adjustl(activate)).ne.'tansig' .and. &
         trim(adjustl(activate)).ne.'algebraic' .and. &
         trim(adjustl(activate)).ne.'purelin' .and. &
         trim(adjustl(activate)).ne.'logsig' .and. &
         trim(adjustl(activate)).ne.'arctan') then
            write(*,*) trim(adjustl(activate))
            stop 'Wrong Activation Fuction!'
      end if

      return
end subroutine check_act
!=================================================================================
subroutine trans(M,n,a,activate)
      implicit none
      real*8, parameter :: fac=2.d0/acos(-1.d0)
      real*8, parameter :: fac1=0.5d0*acos(-1.d0)
      integer, intent(in) :: M
      real*8, intent(in) :: n(M)
      real*8, intent(out) :: a(M)
      character(len=*), intent(in) :: activate

      if(trim(adjustl(activate)).eq.'tansig') then
            a=tanh(n)
      else if(trim(adjustl(activate)).eq.'algebraic') then
            a=n/sqrt(abs(1.d0+n*n))
      else if(trim(adjustl(activate)).eq.'purelin') then
            a=n
      else if(trim(adjustl(activate)).eq.'logsig') then
            a=exp(n)
            a=a/(1.d0+a)
      else if(trim(adjustl(activate)).eq.'arctan') then
            a=fac*atan(fac1*n)
      else
            stop 'Wrong Activation Fuction!'
      end if

      return
end subroutine trans
!=================================================================================
subroutine dtrans(M,n,dn,activate)
      implicit none
      real*8, parameter :: fac1=0.5d0*acos(-1.d0)
      integer, intent(in) :: M
      real*8, intent(in) :: n(M)
      real*8, intent(out) :: dn(M)
      character(len=*), intent(in) :: activate

      if(trim(adjustl(activate)).eq.'tansig') then
            dn=1.d0/cosh(n)**2
      else if(trim(adjustl(activate)).eq.'algebraic') then
            dn=sqrt(abs(1.d0/(1.d0+n*n)))
            dn=dn**3
      else if(trim(adjustl(activate)).eq.'purelin') then
            dn=1.d0
      else if(trim(adjustl(activate)).eq.'logsig') then
            dn=exp(n)
            dn=dn/(1.d0+dn)**2
      else if(trim(adjustl(activate)).eq.'arctan') then
            dn=1.d0/(1.d0+(fac1*n)**2)
      else
            stop 'Wrong Activation Fuction!'
      end if

      return
end subroutine dtrans
!==================================================================================
subroutine d2trans(M,n,d2n,activate)
      implicit none
      real*8, parameter :: fac1=0.5d0*acos(-1.d0)
      integer, intent(in) :: M
      real*8, intent(in) :: n(M)
      real*8, intent(out) :: d2n(M)
      character(len=*), intent(in) :: activate

      if(trim(adjustl(activate)).eq.'tansig') then
            d2n=tanh(n)
            d2n=2.d0*d2n*(d2n**2-1.d0)
      else if(trim(adjustl(activate)).eq.'algebraic') then
            d2n=sqrt(abs(1.d0/(1.d0+n*n)))
            d2n=-3.d0*n*d2n**5
      else if(trim(adjustl(activate)).eq.'purelin') then
            d2n=0.d0
      else if(trim(adjustl(activate)).eq.'logsig') then
            d2n=exp(n)
            d2n=d2n*(1.d0-d2n)/(1.d0+d2n)**3
      else if(trim(adjustl(activate)).eq.'arctan') then
            d2n=-2.d0*fac1*n/(1.d0+(fac1*n)**2)**2
      else
            stop 'Wrong Activation Fuction!'
      end if

      return
end subroutine d2trans
!=================================================================================
subroutine FLUnit(fid)
      !get an available UNIT index for input/output
      implicit none
      integer, intent(out) :: fid
      integer :: i
      logical :: unitex,unitop
      fid=0
      do i=15,99999
            inquire(UNIT=i,EXIST=unitex,OPENED=unitop)
            if(unitex .and. .not. unitop) then
                  fid=i
                  exit
            end if
      end do
      if(fid.eq.0) stop "FLUnit: failed to find an available unit."
      return
end subroutine FLUnit
!==================================================================================
subroutine countline(filename,n)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: n
      integer :: fid,ios

      call FLUnit(fid)
      open(fid,file=trim(adjustl(filename)),iostat=ios)
      if(ios.ne.0) then
        print*, "Error opening ",trim(adjustl(filename))
        n=0
        return
      end if

      n=0
      do
            read (fid,*,end=10)
            n=n+1
      end do
10    close(fid)

      return
end subroutine countline
!==================================================================================

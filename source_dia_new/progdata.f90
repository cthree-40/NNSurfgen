! Shared Program Information
module progdata
! TYPE DEFINITIONS
! Type abpoint     
! [Purpose] 
!   Stores all information of an ab initio data point
! [Fields]
! cgeom      DOUBLE PRECISION,dimension(3*natoms)
!            Cartesian geometry of all the atoms
! igeom      DOUBLE PRECISION,dimension(ncoord)
!            Scaled global internal geometry
! energy     DOUBLE PRECISION,dimension(nstates)
!            Ab initio energies of all the states
! grads      DOUBLE PRECISIOIN,dimension(3*natoms,nstates,nstates)
!            Ab initio energy gradients and derivative couplings
! bmat       DOUBLE PRECISION,dimension(ncoords,3*natoms)
!            Wilson's B matrix from Cartesian to scaled global internal coords,
!            or from local linear internal coordinates to scaled global internal
!            coordinates if useIntGrad==.true. 
! cbmat      DOUBLE PRECISION,dimension(ncoords,3*natoms)
!            B matrix in the original cartesian coordinate
! eval       DOUBLE PRECISION,dimension(3*natoms)
!            eigenvalues of B^t.B
! lmat       DOUBLE PRECISION,dimension(3*natoms,3*natoms)
! scale      DOUBLE PRECISION,dimension(3*natoms)
!            Scaling factor for each of the internal motions
! nvibs      INTEGER
!            lmat Transformation to local internal coordinates, defined as the
!            eigenvectors of B^T.B
!            nvibs Number of internal degree of freedom at this geometry
!            Those eigenvectors with eigenvalue larger than threshold intGradT
!            (stored in MODULE makesurfdata) are considered internal.  
!            The non-internal modes are also stored in the last fields of lmat
! ndeggrp    INTEGER
!            Number of groups of degenerate states. Non-degenerate state counts
!            as one group by itself.
! deg_groups INTEGER,dimension(ndeggrp,2)
!            Stores the index of the first state and number of states for each
!            degenerate group            
! lb,ub      INTEGER
!            Specifies the index range between which ab initio data is available
!      
TYPE abpoint
  INTEGER                                       :: id
  DOUBLE PRECISION,dimension(:),allocatable     :: cgeom
  DOUBLE PRECISION,dimension(:),allocatable     :: igeom
  DOUBLE PRECISION,dimension(:,:),allocatable   :: energy
  DOUBLE PRECISION,dimension(:,:),allocatable   :: hd
  DOUBLE PRECISION,dimension(:),allocatable     :: eval
  DOUBLE PRECISION,dimension(:,:,:),allocatable :: grads
  INTEGER                                       :: lb,ub

  INTEGER                                       :: nvibs
  DOUBLE PRECISION,dimension(:,:),allocatable   :: bmat,cbmat
  DOUBLE PRECISION,dimension(:,:),allocatable   :: lmat
  DOUBLE PRECISION,dimension(:),allocatable     :: scale

  INTEGER,dimension(:,:),allocatable            :: deg_groups
  INTEGER                                       :: ndeggrp
end TYPE abpoint

! CONSTANTS
  real*8, parameter                             :: au2ev  = 27.21138602d0
  real*8, parameter                             :: au2cm  = 219474d0
  real*8, parameter                             :: ev2cm  = 8065.51d0
  real*8, parameter                             :: pi = dacos(-1.d0)

! GLOBAL VARIABLES
  INTEGER                                       :: printlvl
  INTEGER                                       :: natoms

END MODULE progdata
!===================================================================================

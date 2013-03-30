!> @file ex20.f90
!> @brief Test the FORTRAN interface for PEXSI.
!> @author Lin Lin
!> @date 2013-01-31
program ex20
implicit none
include 'mpif.h'

integer :: nrows, nnz, nnzLocal, numColLocal
integer, allocatable, dimension(:) ::  colptrLocal, rowindLocal
double precision, allocatable, dimension(:) :: &
	HnzvalLocal, SnzvalLocal, DMnzvalLocal, EDMnzvalLocal, &
	FDMnzvalLocal, muList, numElectronList, numElectronDrvList
integer :: numPole
double precision :: temperature, numElectronExact, numElectron,&
	gap, deltaE
double precision :: mu, muZeroT, muMin, muMax
integer:: muMaxIter, muIter
integer:: ordering
integer:: isInertiaCount
double precision :: poleTolerance, numElectronTolerance
integer:: npPerPole, nprow, npcol
integer :: mpirank, mpisize, ierr
double precision:: timeSta, timeEnd
character*32 :: Hfile, Sfile
! Communicator for reading the matrix, with size npPerPole
integer:: readComm
integer:: isProcRead
integer:: i

call mpi_init( ierr )
call mpi_comm_rank( MPI_COMM_WORLD, mpirank, ierr )
call mpi_comm_size( MPI_COMM_WORLD, mpisize, ierr )


! Data is for the DNA matrix.
temperature      = 3000.0d0
numElectronExact = 2442.0d0
numPole          = 20
gap              = 0.0d0
! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = 3.0d0
! Initial guess of chemical potential, also updated after pexsi.
mu               = -0.60d0
! Lower/Upper bound for the chemical potential.
muMin            = -1.0d0
muMax            =  0.0d0
! muMaxIter should be 1 or 2 later when combined with SCF.
muMaxIter        = 10
! Do not compute a pole if the corresponding weight is < poleTolerance.
poleTolerance    = 1d-8
! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
numElectronTolerance = 1d-1
! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
npPerPole        = 16
Hfile            = "H.matrix"
Sfile            = "S.matrix"
! Ordering 
!   0   : PARMETIS
!   1   : METIS_AT_PLUS_A
!   2   : MMD_AT_PLUS_A
ordering         = 1
isInertiaCount   = 0


! Read and compute the size/local size of the arrays 
! The conversion of the string to fit the C format is important.

! Split the processors to read matrix
if( mpirank < npPerPole ) then
	isProcRead = 1
else
	isProcRead = 0
endif

call mpi_comm_split( MPI_COMM_WORLD, isProcRead, mpirank, readComm, ierr )


allocate( muList( muMaxIter ) )
allocate( numElectronList( muMaxIter ) )
allocate( numElectronDrvList( muMaxIter ) )

if( isProcRead == 1 ) then
	write(*,*) "Proc ", mpirank, " is reading file..."
	call f_read_distsparsematrix_formatted_head( &
		trim(Hfile)//char(0),&
		nrows,&
		nnz,&
		nnzLocal,&
		numColLocal,&
		readComm )

	if( mpirank .eq. 0 ) then
		write(*,*) "Matrix size (local data on proc 0):" 
		write(*,*) "size = ", nrows
		write(*,*) "nnz  = ", nnz
		write(*,*) "nnzLocal = ", nnzLocal
		write(*,*) "numColLocal = ", numColLocal
	endif

	! Allocate memory
	allocate( colptrLocal( numColLocal + 1 ) )
	allocate( rowindLocal( nnzLocal ) )
	allocate( HnzvalLocal( nnzLocal ) )
	allocate( SnzvalLocal( nnzLocal ) ) 
	allocate( DMnzvalLocal( nnzLocal ) ) 
	allocate( EDMnzvalLocal( nnzLocal ) ) 
	allocate( FDMnzvalLocal( nnzLocal ) ) 

	timeSta = mpi_wtime()

	call f_read_distsparsematrix_formatted (&
		trim(Hfile)//char(0),&
		nrows,&
		nnz,&
		nnzLocal,&
		numColLocal,&
		colptrLocal,&
		rowindLocal,&
		HnzvalLocal,&
		readComm )

	call f_read_distsparsematrix_formatted (&
		trim(Sfile)//char(0),&
		nrows,&
		nnz,&
		nnzLocal,&
		numColLocal,&
		colptrLocal,&
		rowindLocal,&
		SnzvalLocal,&
	  readComm )	

	timeEnd = mpi_wtime()

endif

call mpi_barrier( MPI_COMM_WORLD, ierr )
call mpi_comm_free( readComm, ierr )

if( mpirank == 0 ) then
	write(*,*) "Time for reading H/S matrices is ", &
		timeEnd - timeSta, " [s]"
endif


call f_ppexsi_interface( &
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
  rowindLocal,&
	HnzvalLocal,&
	SnzvalLocal,&
	DMnzvalLocal,&
	EDMnzvalLocal,&
	FDMnzvalLocal,&
	numPole,&
	temperature,&
	numElectronExact,&
	numElectron,&
	gap,&
	deltaE,&
	mu,&
	muMin,&
  muMax,&
	muMaxIter,&
	ordering,&
	isInertiaCount,&
	muIter,&
	muList,&
	numElectronList,&
	numElectronDrvList,&
	muZeroT,&
	poleTolerance,&
	numElectronTolerance,&
	MPI_COMM_WORLD,&
	npPerPole )

if( mpirank == 0 ) then
	write(*, *) "mu          = ", mu
	write(*, *) "guess of mu(T=0) = ", muZeroT
	write(*, *) "numElectron = ", numElectron
	write(*, *) "number of mu iterations = ", muIter
	write(*, *) "mu history"
	do i = 1, muIter
		write(*,*) muList(i)
	enddo
	write(*, *) "numElectron history"
	do i = 1, muIter
		write(*,*) numElectronList(i)
	enddo
	write(*, *) "dN/dMu history"
	do i = 1, muIter
		write(*,*) numElectronDrvList(i)
	enddo
endif

deallocate( muList )
deallocate( numElectronList )
deallocate( numElectronDrvList )

if( isProcRead == 1 ) then
	deallocate( colptrLocal )
	deallocate( rowindLocal )
	deallocate( HnzvalLocal )
	deallocate( SnzvalLocal )
	deallocate( DMnzvalLocal )
	deallocate( EDMnzvalLocal )
	deallocate( FDMnzvalLocal )
endif


call mpi_finalize( ierr )

end program ex20


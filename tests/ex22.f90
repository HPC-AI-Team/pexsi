!> @file ex22.f90
!> @brief Test the new FORTRAN interface for PPEXSI.
!> @author Lin Lin
!> @date 2013-04-10
program ex22
implicit none
include 'mpif.h'

integer :: nrows, nnz, nnzLocal, numColLocal
integer, allocatable, dimension(:) ::  colptrLocal, rowindLocal
double precision, allocatable, dimension(:) :: &
	HnzvalLocal, SnzvalLocal, DMnzvalLocal, EDMnzvalLocal, &
	FDMnzvalLocal, muList, numElectronList, numElectronDrvList,&
	shiftList, inertiaList
integer :: numPole
double precision :: K2au, temperature, numElectronExact, numElectron,&
	gap, deltaE
double precision ::   muMin0, muMax0, muInertia, muMinInertia, muMaxInertia,&
	muLowerEdge, muUpperEdge, muPEXSI, muMinPEXSI, muMaxPEXSI
integer:: inertiaMaxIter, inertiaIter, muMaxIter, muIter
integer:: ordering
integer:: isInertiaCount
double precision :: poleTolerance, PEXSINumElectronTolerance, &
	inertiaNumElectronTolerance
integer:: npPerPole, nprow, npcol
integer :: mpirank, mpisize, ierr
double precision:: timeSta, timeEnd
character*32 :: Hfile, Sfile
integer:: isSIdentity
! Communicator for reading the matrix, with size npPerPole
integer:: readComm
integer:: isProcRead
integer:: i

call mpi_init( ierr )
call mpi_comm_rank( MPI_COMM_WORLD, mpirank, ierr )
call mpi_comm_size( MPI_COMM_WORLD, mpisize, ierr )


! Data is for the g20 matrix.
K2au             = 3.1668152d-6
! Temperature should be in the same unit as the H matrix. Here it is Hartree.
temperature      = 1000.0d0 * K2au

numElectronExact = 13.0d0
numPole          = 60
gap              = 0.0d0
! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = 20.0d0
! Initial searching interval for the chemical potential for inertia count.
muMin0           =  0.0d0
muMax0           =  2.0d0
! Maximum number of iterations for computing the inertia
inertiaMaxIter   = 3
! Maximum number of iterations for PEXSI iteration
muMaxIter        = 3
! Do not compute a pole if the corresponding weight is < poleTolerance.
poleTolerance    = 1d-8
! Stop inertia count if Ne(muMax) - Ne(muMin) < inertiaNumElectronTolerance
inertiaNumElectronTolerance = 10
! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
PEXSINumElectronTolerance = 1d-2
! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
npPerPole        = 1
Hfile            = "g20.matrix"
! Empty Sfile means the overlap matrix is identity
Sfile            = ''
isSIdentity      = 1

! Ordering 
!   0   : PARMETIS
!   1   : METIS_AT_PLUS_A
!   2   : MMD_AT_PLUS_A
ordering         = 2


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

allocate( shiftList( numPole ) )
allocate( inertiaList( numPole ) )


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

	if( isSIdentity == 0 ) then
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
	endif

	timeEnd = mpi_wtime()

endif

call mpi_barrier( MPI_COMM_WORLD, ierr )
call mpi_comm_free( readComm, ierr )

if( mpirank == 0 ) then
	write(*,*) "Time for reading H/S matrices is ", &
		timeEnd - timeSta, " [s]"
endif


call f_ppexsi_inertiacount_interface(&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
	rowindLocal,&
	HnzvalLocal,&
	isSIdentity,&
	SnzvalLocal,&
	temperature,&
	numElectronExact,&
	muMin0,&
	muMax0,&
	numPole,&
	inertiaMaxIter,&
	inertiaNumElectronTolerance,&
	ordering,&
	npPerPole,&
	MPI_COMM_WORLD,&
	muMinInertia,&
	muMaxInertia,&
	muLowerEdge,&
	muUpperEdge,&
	muInertia,&
	inertiaIter,&
	shiftList,&
	inertiaList)

call f_ppexsi_solve_interface(&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
  rowindLocal,&
	HnzvalLocal,&
	isSIdentity,&
	SnzvalLocal,&
	temperature,&
	numElectronExact,&
	muInertia,&
	muMinInertia,&
	muMaxInertia,&
	gap,&
	deltaE,&
	numPole,&
	muMaxIter,&
	PEXSINumElectronTolerance,&
	poleTolerance,&
	ordering,&
	npPerPole,&
	MPI_COMM_WORLD,&
	DMnzvalLocal,&
	EDMnzvalLocal,&
	FDMnzvalLocal,&
	muPEXSI,&
	numElectron,&
	muMinPEXSI,&
	muMaxPEXSI,&
	muIter,&
	muList,&
	numElectronList,&
	numElectronDrvList)

if( mpirank == 0 ) then
	write(*, *) "After inertia count,"
	write(*, *) "muMinInertia  = ", muMinInertia
	write(*, *) "muMaxInertia  = ", muMaxInertia
	write(*, *) "muLowerEdge   = ", muLowerEdge
	write(*, *) "muUpperEdge   = ", muUpperEdge
	write(*, *) "inertiaIter   = ", inertiaIter
	! write(*, *) "shift ,           inertia count"
	! do i = 1, numPole 
		! write(*,*) shiftList(i), inertiaList(i)
	! enddo
	write(*, *) 
	write(*, *) "After PEXSI,"
	write(*, *) "muPEXSI       = ", muPEXSI
	write(*, *) "numElectron   = ", numElectron
	write(*, *) "muMinPEXSI    = ", muMinPEXSI
	write(*, *) "muMaxPEXSI    = ", muMaxPEXSI
	write(*, *) "muIter        = ", muIter
endif

deallocate( muList )
deallocate( numElectronList )
deallocate( numElectronDrvList )

deallocate( shiftList )
deallocate( inertiaList )

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

end program ex22


!> @file ex19.f90
!> @brief Test the interface between FORTRAN and PEXSI.
!> @author Lin Lin
!> @date 2013-01-31
program ex19
implicit none
include 'mpif.h'

integer :: a
integer :: ierr

call mpi_init( ierr )

a = 2

call f_dummy_interface( MPI_COMM_WORLD, a )

call mpi_finalize( ierr )

end program ex19


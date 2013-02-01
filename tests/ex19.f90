!> @file ex19.f90
!> @brief Test the interface between FORTRAN and PEXSI.
!> @author Lin Lin
!> @date 2013-01-31
program ex19
implicit none

integer :: a

a = 2

call f_dummy_interface( a )


end program ex19


!	 Copyright (c) 2012 The Regents of the University of California,
!	 through Lawrence Berkeley National Laboratory.  
!
!  Author: Lin Lin
!	 
!  This file is part of PEXSI. All rights reserved.
!
!	 Redistribution and use in source and binary forms, with or without
!	 modification, are permitted provided that the following conditions are met:
!
!	 (1) Redistributions of source code must retain the above copyright notice, this
!	 list of conditions and the following disclaimer.
!	 (2) Redistributions in binary form must reproduce the above copyright notice,
!	 this list of conditions and the following disclaimer in the documentation
!	 and/or other materials provided with the distribution.
!	 (3) Neither the name of the University of California, Lawrence Berkeley
!	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
!	 be used to endorse or promote products derived from this software without
!	 specific prior written permission.
!
!	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
!	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!	 You are under no obligation whatsoever to provide any bug fixes, patches, or
!	 upgrades to the features, functionality or performance of the source code
!	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
!	 available either publicly, or directly to Lawrence Berkeley National
!	 Laboratory, without imposing a separate written license agreement for such
!	 Enhancements, then you hereby grant the following license: a non-exclusive,
!	 royalty-free perpetual license to install, use, modify, prepare derivative
!	 works, incorporate into other computer software, distribute, and sublicense
!	 such enhancements or derivative works thereof, in binary and source code form.
!
!> @file f_interface.f90
!> @brief FORTRAN interface for PEXSI library using ISO_C_BINDING.
!>
!> The ISO_C_BINDING feature is included in the FORTRAN 2003 standard, and is
!> implemented in most modern compilers.
!> 
!>
!> @note (From Alberto Garcia, 2013-04-10) Array arguments are *required* by the
!> standard to be explicit shape (e.g. a(3)), or assumed size (e.g. a(*)). This
!> avoids problems with the assumed shape specification (e.g. a(:)), which would
!> need to pass extra information in general. This permits the use of pointers
!> and array sections as actual arguments. The compiler would presumably make
!> on-the-fly copies in the case of non-contiguous data (a situation that should
!> be avoided on performance grounds).
!>
!> @see  c_pexsi_interface.h
!> @date 2014-04-01
module f_pexsi_interface

interface  
  subroutine f_read_distsparsematrix_formatted_head (&
      filename,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      fcomm ) &
      bind(C, Name="f_read_distsparsematrix_formatted_head")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(out)        :: nrows, nnz, nnzLocal, numColLocal
    integer, intent(in)                :: fcomm
  end subroutine 
end interface

interface  
  subroutine f_read_distsparsematrix_formatted(&
      filename,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      nzvalLocal,&
      fcomm ) &
      bind(C, Name="f_read_distsparsematrix_formatted")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(in), value  :: nrows, nnz, nnzLocal, numColLocal
    integer(c_int), intent(out)        :: colptrLocal(*), rowindLocal(*)
    real(c_double), intent(out)        :: nzvalLocal(*)
    integer, intent(in)                :: fcomm
  end subroutine 
end interface

interface  
  function f_ppexsi_plan_initialize(&
      fcomm,&
      numProcRow,&
      numProcCol,&
      outputFileIndex,&
      info) &
      bind(C, Name="f_ppexsi_plan_initialize")
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in)                :: fcomm
    integer(c_int), intent(in), value  :: numProcRow, numProcCol,outputFileIndex
    integer(c_int), intent(out)        :: info
    integer(c_intptr_t)                :: f_ppexsi_plan_initialize
  end function
end interface

interface  
  subroutine f_ppexsi_plan_finalize(&
      plan,&
      info) &
      bind(C, Name="PPEXSIPlanFinalize")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    integer(c_int), intent(out)            :: info
  end subroutine 
end interface


end module f_pexsi_interface

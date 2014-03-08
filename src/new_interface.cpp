/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Author: Lin Lin

   This file is part of PEXSI. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   (1) Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
   (2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   You are under no obligation whatsoever to provide any bug fixes, patches, or
   upgrades to the features, functionality or performance of the source code
   ("Enhancements") to anyone; however, if you choose to make your Enhancements
   available either publicly, or directly to Lawrence Berkeley National
   Laboratory, without imposing a separate written license agreement for such
   Enhancements, then you hereby grant the following license: a non-exclusive,
   royalty-free perpetual license to install, use, modify, prepare derivative
   works, incorporate into other computer software, distribute, and sublicense
   such enhancements or derivative works thereof, in binary and source code form.
*/
/// @file new_interface.cpp
/// @brief New interface subroutines of PPEXSI that can be called by both C and FORTRAN.
///
/// This file will eventually merge with interface.cpp.
/// @date 2014-03-07
#include "c_pexsi_interface.h"
#include "c_pexsi_new_interface.h"
#include "ppexsi.hpp"
#include "blas.hpp"

// FIXME
// Error handling used in the C interface that is different from the
// throw/catch system.
#define iC(fun)  { int ierr=fun; if(ierr!=0) exit(1); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<std::endl; std::cerr.flush(); exit(1); } }

using namespace PEXSI;

extern "C"
PPEXSIPlan PPEXSIPlanInitialize(
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       HnzvalLocal,                  
    int           isSIdentity,                  
    double*       SnzvalLocal,
    int           nprow,
    int           npcol,
    MPI_Comm      comm,
    int           outputFileIndex ){

	Int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );

  Int npPerPole = nprow * npcol;
  
  PPEXSIData *ptrData;
  GridType gridPole( comm, mpisize / npPerPole, npPerPole );
  ptrData = new PPEXSIData( &gridPole, nprow, npcol );
        
  return reinterpret_cast<PPEXSIPlan>(ptrData);
}   // -----  end of function PPEXSIPlanInitialize  ----- 

extern "C"
void PPEXSIPlanFinalize( PPEXSIPlan plan ){
  delete reinterpret_cast<PPEXSIData*>(plan);
}   // -----  end of function PPEXSIPlanFinalize  ----- 

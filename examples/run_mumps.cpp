/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Authors: Lin Lin and Mathias Jacquelin

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
/// @file run_pselinv.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @date 2013-04-15
#include  "ppexsi.hpp"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654


#define INFO(I) info[(I)-1]
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */


#include "timer.h"

//#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#include "zmumps_c.h"
#define MUMPS(a) zmumps_c(a)
#define MUMPS_STRUC_C ZMUMPS_STRUC_C
#else
#define MYSCALAR Real
#include "dmumps_c.h"
#define MUMPS(a) dmumps_c(a)
#define MUMPS_STRUC_C DMUMPS_STRUC_C
#endif


using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout << "Usage" << std::endl << "run_pselinv -T [isText] -F [doFacto -E [doTriSolve] -Sinv [doSelInv]]  -H <Hfile> -S [Sfile] -colperm [colperm] -r [nprow] -c [npcol] -npsymbfact [npsymbfact] -P [maxpipelinedepth] -SinvBcast [doSelInvBcast] -SinvPipeline [doSelInvPipeline] -SinvHybrid [doSelInvHybrid] -rshift [real shift] -ishift [imaginary shift] -ToDist [doToDist] -Diag [doDiag]" << std::endl;
}

int main(int argc, char **argv) 
{
  MUMPS_STRUC_C id;

  if( argc < 3 ) {
    Usage();
    return 0;
  }

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif

  MPI_Init( &argc, &argv );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


  try{
    MPI_Comm world_comm;

    // *********************************************************************
    // Input parameter
    // *********************************************************************
    std::map<std::string,std::string> options;

    OptionsCreate(argc, argv, options);

    // Default processor number
    Int nprow = 1;
    Int npcol = mpisize;

    if( options.find("-r") != options.end() ){
      if( options.find("-c") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -c option, -r also needs to be provided." );
      }
    }

    //Create a communicator with npcol*nprow processors
    MPI_Comm_split(MPI_COMM_WORLD, mpirank<nprow*npcol, mpirank, &world_comm);

    if (mpirank<nprow*npcol){

      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );
#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
      TAU_PROFILE_SET_CONTEXT(world_comm);
#endif

      stringstream  ss;
      ss << "logTest" << mpirank;
      statusOFS.open( ss.str().c_str() );

      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  throw std::runtime_error( "nprow == npcol is assumed in this test routine." );
      //}

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

      std::string Hfile, Sfile;
      int isCSC = true;
      if( options.find("-T") != options.end() ){ 
        isCSC= ! atoi(options["-T"].c_str());
      }


      /* load A in coordinates version */
      /* a[k] is as row irn[k] and col jcn[k] */

      IntNumVec irn, jcn;
      //  int * rhs_irn, * rhs_jcn;
      //  struct sparse_matrix_t* A;
      //
      //  int nz_rhs;
      //
      //  double * rhs;
      //  int * irhs_sparse;
      //  int * irhs_ptr;




      int icntl27 = 8;
      if(options.find("-B") != options.end() ){
        icntl27 = atoi(options["-B"].c_str());
      }

      int checkAccuracy = false;
      if( options.find("-E") != options.end() ){ 
        checkAccuracy= atoi(options["-E"].c_str());
      }

      int doFacto = true;
      if( options.find("-F") != options.end() ){ 
        doFacto= atoi(options["-F"].c_str());
      }

      int doSelInv = 1;
      if( options.find("-Sinv") != options.end() ){ 
        doSelInv= atoi(options["-Sinv"].c_str());
      }

      int doSymbfact = true;
      if( options.find("-Symb") != options.end() ){ 
        doSymbfact= atoi(options["-Symb"].c_str());
      }

      int doToDist = true;
      if( options.find("-ToDist") != options.end() ){ 
        doToDist= atoi(options["-ToDist"].c_str());
      }

      int doDiag = false;
      if( options.find("-Diag") != options.end() ){ 
        doDiag = atoi(options["-Diag"].c_str());
      }

      doFacto = doFacto && doSymbfact;

      if( options.find("-H") != options.end() ){ 
        Hfile = options["-H"];
      }
      else{
        throw std::logic_error("Hfile must be provided.");
      }

      if( options.find("-S") != options.end() ){ 
        Sfile = options["-S"];
      }
      else{
        statusOFS << "-S option is not given. " 
          << "Treat the overlap matrix as an identity matrix." 
          << std::endl << std::endl;
      }

      Int numProcSymbFact;
      if( options.find("-npsymbfact") != options.end() ){ 
        numProcSymbFact = atoi( options["-npsymbfact"].c_str() );
      }
      else{
        statusOFS << "-npsymbfact option is not given. " 
          << "Use default value (maximum number of procs)." 
          << std::endl << std::endl;
        numProcSymbFact = 0;
      }

      Real rshift = 0.0, ishift = 0.0;
      if( options.find("-rshift") != options.end() ){ 
        rshift = atof(options["-rshift"].c_str());
      }
      if( options.find("-ishift") != options.end() ){ 
        ishift = atof(options["-ishift"].c_str());
      }

      /*
         std::string ColPerm;
         if( options.find("-colperm") != options.end() ){ 
         ColPerm = options["-colperm"];
         }
         else{
         statusOFS << "-colperm option is not given. " 
         << "Use MMD_AT_PLUS_A." 
         << std::endl << std::endl;
         ColPerm = "MMD_AT_PLUS_A";
         }
       */


      // *********************************************************************
      // Read input matrix
      // *********************************************************************

      // Setup grid.
      SuperLUGrid<MYSCALAR> g( world_comm, nprow, npcol );
      //      SuperLUGrid<Complex> g1( world_comm, nprow, npcol );

      int      m, n;
      DistSparseMatrix<MYSCALAR>  AMat;

      DistSparseMatrix<Real> HMat;
      DistSparseMatrix<Real> SMat;
      Real timeSta, timeEnd;
      GetTime( timeSta );
      if(isCSC)
        ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
      else{
        ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
        ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
      }

      if( Sfile.empty() ){
        // Set the size to be zero.  This will tell PPEXSI.Solve to treat
        // the overlap matrix as an identity matrix implicitly.
        SMat.size = 0;  
      }
      else{
        if(isCSC)
          ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
        else{
          ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, world_comm ); 
          ParaWriteDistSparseMatrix( "S.csc", SMat, world_comm ); 
        }
      }

      GetTime( timeEnd );
      LongInt nnzH = HMat.Nnz();
      if( mpirank == 0 ){
        cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
        cout << "H.size = " << HMat.size << endl;
        cout << "H.nnz  = " << nnzH  << endl;
      }

      // Get the diagonal indices for H and save it n diagIdxLocal_

      std::vector<Int>  diagIdxLocal;
      { 
        Int numColLocal      = HMat.colptrLocal.m() - 1;
        Int numColLocalFirst = HMat.size / mpisize;
        Int firstCol         = mpirank * numColLocalFirst;

        diagIdxLocal.clear();

        for( Int j = 0; j < numColLocal; j++ ){
          Int jcol = firstCol + j + 1;
          for( Int i = HMat.colptrLocal(j)-1; 
              i < HMat.colptrLocal(j+1)-1; i++ ){
            Int irow = HMat.rowindLocal(i);
            if( irow == jcol ){
              diagIdxLocal.push_back( i );
            }
          }
        } // for (j)
      }


      GetTime( timeSta );

      AMat.size          = HMat.size;
      AMat.nnz           = HMat.nnz;
      AMat.nnzLocal      = HMat.nnzLocal;
      AMat.colptrLocal   = HMat.colptrLocal;
      AMat.rowindLocal   = HMat.rowindLocal;
      AMat.nzvalLocal.Resize( HMat.nnzLocal );
      AMat.comm = world_comm;

      MYSCALAR *ptr0 = AMat.nzvalLocal.Data();
      Real *ptr1 = HMat.nzvalLocal.Data();
      Real *ptr2 = SMat.nzvalLocal.Data();

#ifdef _MYCOMPLEX_
      Complex zshift = Complex(rshift, ishift);
#else
      Real zshift = Real(rshift);
#endif

      if( SMat.size != 0 ){
        // S is not an identity matrix
        for( Int i = 0; i < HMat.nnzLocal; i++ ){
          AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
        }
      }
      else{
        // S is an identity matrix
        for( Int i = 0; i < HMat.nnzLocal; i++ ){
          AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
        }

        for( Int i = 0; i < diagIdxLocal.size(); i++ ){
          AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
        }
      } // if (SMat.size != 0 )

      LongInt nnzA = AMat.Nnz();
      if( mpirank == 0 ){
        cout << "nonzero in A (DistSparseMatrix format) = " << nnzA << endl;
      }


      GetTime( timeEnd );

      // Compute the number of columns on each processor
      Int numColFirst;
      numColFirst = AMat.size / mpisize;


      if( mpirank == 0 )
        cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;


      IntNumVec irhs_ptr, irhs_sparse;
      NumVec<MYSCALAR> rhs_sparse;

      //MUMPS conversion
#define ICNTL18_3



#ifdef ICNTL18_3
      irn.Resize(AMat.nnzLocal);
      jcn.Resize(AMat.nnzLocal);
      //extract the coordinates
      Int index = 0;
      for(Int j =0; j<AMat.colptrLocal.m()-1; ++j){
        Int col = j +1 + (mpirank)*numColFirst; 
        for(Int i = AMat.colptrLocal[j]; i<AMat.colptrLocal[j+1]; ++i){
          Int row = AMat.rowindLocal[i-1];
          irn[index] = row;
          jcn[index] = col;
          index++;
        } 
      }


      //Get the structure of the matrix to compute the sparse rhs on the host

      IntNumVec numColLocalVec(mpisize);
      IntNumVec displs(mpisize);

      // Compute the number of columns on each processor
      numColLocalVec.Resize(mpisize);
      Int numColLocal;
      numColFirst = AMat.size / mpisize;
      SetValue( numColLocalVec, numColFirst );
      numColLocalVec[mpisize-1] = AMat.size - numColFirst * (mpisize-1);  // Modify the last entry	
      numColLocal = numColLocalVec[mpirank];

      displs[0]=0;
      for(int p = 1;p<mpisize;++p){
        displs[p]= displs[p-1]+ numColLocalVec[p-1];
      }

      for(int p = 0;p<mpisize;++p){
        numColLocalVec[p]*=sizeof(Int);
        displs[p]*=sizeof(Int);
      }


      if(mpirank==0){
        irhs_ptr.Resize(AMat.size+1);
        irhs_sparse.Resize(nnzA);
        rhs_sparse.Resize(nnzA);
      }

      //gatherv the colptrs
      if(mpirank==0){
        MPI_Gatherv(&AMat.colptrLocal[0],numColLocalVec[mpirank],MPI_BYTE,&irhs_ptr[0],&numColLocalVec[0],&displs[0],MPI_BYTE,0,world_comm);
      }
      else{
        MPI_Gatherv(&AMat.colptrLocal[0],numColLocalVec[mpirank],MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,0,world_comm);
      }

      if(mpirank==0){
        Int index = 0;

        for(int p=0;p<mpisize;++p){
          numColLocalVec[p] /= sizeof(Int); 
          displs[p] /= sizeof(Int); 
        }

        for(int p=1;p<mpisize;++p){
          Int * colptr = &irhs_ptr[0] + displs[p];
          Int * prevColptr = &irhs_ptr[0] + displs[p-1];
          Int offset = prevColptr[ numColLocalVec[p-1] -1 ]; 
          if(numColLocalVec[p]>0){
            for(int j = numColLocalVec[p]-1; j>=0;--j){
              colptr[j] += offset - colptr[0]; 
            }
          }
        }
        irhs_ptr[AMat.size]=nnzA+1;
      }

      //do the same for rowind
      //algather nzvalLocal
      IntNumVec nnzLocalVec(mpisize);
      MPI_Gather(&AMat.nnzLocal,sizeof(Int),MPI_BYTE,&nnzLocalVec[0],sizeof(Int),MPI_BYTE,0,world_comm);
      //compute displacements 
      displs[0]=0;
      for(int p = 1;p<mpisize;++p){
        displs[p]= displs[p-1]+ nnzLocalVec[p-1];
      }

      for(int p = 0;p<mpisize;++p){
        nnzLocalVec[p]*=sizeof(Int);
        displs[p]*=sizeof(Int);
      }

      if(mpirank==0){
        MPI_Gatherv(&AMat.rowindLocal[0],AMat.nnzLocal*sizeof(Int),MPI_BYTE,&irhs_sparse[0],&nnzLocalVec[0],&displs[0],MPI_BYTE,0,world_comm);
      }
      else{
        MPI_Gatherv(&AMat.rowindLocal[0],AMat.nnzLocal*sizeof(Int),MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,0,world_comm);
      }



      MPI_Barrier(world_comm);



      /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
      id.job=JOB_INIT; 
      id.par=1; 
      id.sym=2; /* symmetric matrix*/
      id.comm_fortran=world_comm;
      MUMPS(&id);

      id.n = AMat.size; 
      id.nz = nnzA; 
      id.nz_loc = AMat.nnzLocal;
      id.irn_loc = &irn[0];
      id.jcn_loc = &jcn[0];
      id.a_loc = &AMat.nzvalLocal[0];


      id.ICNTL(18)=3;

      GetTime( timeSta );
      id.job=1;
      /* call the mumps package for symbfact*/
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Symbolic factorization done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }

      GetTime( timeSta );
      id.job=2;
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Factorization done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }




#endif

#ifdef ICNTL18_1
      //convert from CSC to coordinates format
      if(mpirank==0){
        irn.Resize(nnzA);
        jcn.Resize(nnzA);

        irhs_ptr.Resize(AMat.size+1);
        irhs_sparse.Resize(nnzA);
        rhs_sparse.Resize(nnzA);

        //index for irn and jcn
        Int index = 0;
        Int index_ptr = 0;
        for(int p=0;p<mpisize;++p){
          std::stringstream sstr;
          DistSparseMatrix<MYSCALAR> tmpA;
          IntNumVec * pcolptr;
          IntNumVec * prowind;
          if(p>0){
            mpi::Recv(sstr,p,0,1,world_comm);
            deserialize(tmpA,sstr,NO_MASK);
            pcolptr = &tmpA.colptrLocal;
            prowind = &tmpA.rowindLocal;
          }
          else{
            pcolptr = &AMat.colptrLocal;
            prowind = &AMat.rowindLocal;
          }
          //extract the coordinates
          for(Int j =0; j<pcolptr->m()-1; ++j){
            Int col = j +1 + (p)*numColFirst; 
            irhs_ptr[index_ptr++] = index + 1;
            for(Int i = (*pcolptr)[j]; i<(*pcolptr)[j+1]; ++i){
              Int row = (*prowind)[i-1];
              irhs_sparse[index] = row;
              irn[index] = row;
              jcn[index] = col;
              index++;
            } 
          } 
        }
        irhs_ptr[index_ptr] = nnzA+1;
      }
      else{
        //send my AMat.colptrLocal / rowindLocal to P0
        std::stringstream sstr;
        serialize(AMat,sstr,NO_MASK);
        mpi::Send(sstr,0,0,1,world_comm);
      }


      /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
      id.job=JOB_INIT; 
      id.par=1; 
      id.sym=2; /* symmetric matrix*/
      id.comm_fortran=world_comm;
      MUMPS(&id);

      if(mpirank==0){
        id.n = AMat.size; 
        id.nz = nnzA; 
        id.irn=&irn[0]; 
        id.jcn=&jcn[0];
      }
      id.ICNTL(18)=1;

      GetTime( timeSta );
      id.job=1;
      /* call the mumps package for symbfact*/
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Symbolic factorization done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }


      IntNumVec mapping(nnzA);
      //Broadcast the mapping
      MPI_Bcast(&mapping[0],mapping.ByteSize(),MPI_BYTE,0,world_comm);


      //distribute data accordingly
      std::vector<Int> newIRN, newJCN;
      std::vector<MYSCALAR> newVal;

      Int index = 0;
      for(Int pSource = 0; pSource<mpisize;pSource++){
        std::stringstream sstm;
        if(mpirank==pSource){
          serialize(AMat,sstm,NO_MASK);
        }

        Int sstmsize = Size(sstm);
        MPI_Bcast(&sstmsize,sizeof(sstmsize),MPI_BYTE,pSource,world_comm);

        std::vector<char> sstr(sstmsize);
        if(mpirank==pSource){ 
          sstm.read(&sstr[0],sstmsize);
        }
        MPI_Bcast(&sstr[0],sstmsize,MPI_BYTE,pSource,world_comm);
        if(mpirank!=pSource){ 
          sstm.write(&sstr[0],sstmsize);
        }

        //rewind the stringstream
        if(mpirank==pSource){ 
          sstm.seekg (0, std::ios::beg);
        }

        DistSparseMatrix<MYSCALAR> tmpA;
        deserialize(tmpA,sstm,NO_MASK);

        for(Int j =0; j<tmpA.colptrLocal.m()-1; ++j){
          Int col = j +1 + (pSource)*numColFirst; 
          for(Int i = tmpA.colptrLocal[j]; i<tmpA.colptrLocal[j+1]; ++i){
            Int row = tmpA.rowindLocal[i-1];
            if(mapping[index]==mpirank){
              newIRN.push_back(row);
              newJCN.push_back(col);
              newVal.push_back(tmpA.nzvalLocal[i-1]);
            }
            index++;
          } 
        } 
      }

      id.nz_loc = newVal.size();
      id.irn_loc = &newIRN[0];
      id.jcn_loc = &newJCN[0];
      id.a_loc = &newVal[0];
      id.ICNTL(18)=1;

      GetTime( timeSta );
      id.job=2;
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Factorization done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }



#endif




      //proceed with inversion
      if(mpirank==0){
        id.nrhs =AMat.size; 
        id.rhs_sparse = &rhs_sparse[0];
        id.nz_rhs = nnzA;
        id.irhs_sparse = &irhs_sparse[0];
        id.irhs_ptr = &irhs_ptr[0];
      }



      GetTime( timeSta );
      id.ICNTL(11)=1;
      /* ask for inversion*/
      id.ICNTL(30)=1;
      //block size for the solve, ICNTL(27) is the parameter to be tweaked
      id.ICNTL(27) = icntl27; 
      /* Call the MUMPS package. */
      id.job=3;
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Inversion done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }




      id.job=JOB_END; MUMPS(&id); /* Terminate instance */





      /////
      /////      // *********************************************************************
      /////      // Symbolic factorization 
      /////      // *********************************************************************
      /////
      /////      GetTime( timeSta );
      /////      SuperLUOptions luOpt;
      /////      luOpt.ColPerm = ColPerm;
      /////      luOpt.maxPipelineDepth = maxPipelineDepth;
      /////      luOpt.numProcSymbFact = numProcSymbFact;
      /////
      /////
      /////      SuperLUMatrix<MYSCALAR> luMat(g, luOpt );
      /////
      /////
      /////      luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
      /////      GetTime( timeEnd );
      /////      if( mpirank == 0 )
      /////        cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;
      /////
      /////
      /////      if(doSymbfact){
      /////        GetTime( timeSta );
      /////        luMat.SymbolicFactorize();
      /////        luMat.DestroyAOnly();
      /////        GetTime( timeEnd );
      /////
      /////        if( mpirank == 0 )
      /////          cout << "Time for performing the symbolic factorization is " << timeEnd - timeSta << endl;
      /////      }
      /////
      /////      // *********************************************************************
      /////      // Numerical factorization only 
      /////      // *********************************************************************
      /////
      /////      if(doFacto){
      /////        Real timeTotalFactorizationSta, timeTotalFactorizationEnd; 
      /////
      /////
      /////        // Important: the distribution in pzsymbfact is going to mess up the
      /////        // A matrix.  Recompute the matrix A here.
      /////        luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
      /////
      /////        GetTime( timeTotalFactorizationSta );
      /////
      /////        GetTime( timeSta );
      /////        luMat.Distribute();
      /////        GetTime( timeEnd );
      /////        if( mpirank == 0 )
      /////          cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl; 
      /////
      /////
      /////
      /////        GetTime( timeSta );
      /////        luMat.NumericalFactorize();
      /////        GetTime( timeEnd );
      /////
      /////        if( mpirank == 0 )
      /////          cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 
      /////
      /////        GetTime( timeTotalFactorizationEnd );
      /////        if( mpirank == 0 )
      /////          cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 
      /////
      /////
      /////        // *********************************************************************
      /////        // Test the accuracy of factorization by solve
      /////        // *********************************************************************
      /////
      /////        if( checkAccuracy ) {
      /////          SuperLUMatrix<MYSCALAR> A1( g, luOpt );
      /////          SuperLUMatrix<MYSCALAR> GA( g, luOpt );
      /////          A1.DistSparseMatrixToSuperMatrixNRloc( AMat );
      /////          A1.ConvertNRlocToNC( GA );
      /////
      /////          int n = A1.n();
      /////          int nrhs = 5;
      /////          NumMat<MYSCALAR> xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
      /////          NumMat<MYSCALAR> xTrueLocal, bLocal;
      /////          DblNumVec berr;
      /////          UniformRandom( xTrueGlobal );
      /////
      /////          GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );
      /////
      /////          A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
      /////          A1.DistributeGlobalMultiVector( bGlobal,     bLocal );
      /////
      /////          luMat.SolveDistMultiVector( bLocal, berr );
      /////          luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );
      /////        }
      /////
      /////        // *********************************************************************
      /////        // Selected inversion
      /////        // *********************************************************************
      /////
      /////        for(int i=1; i<= doSelInv; ++i )
      /////        {
      /////
      /////          Real timeTotalSelInvSta, timeTotalSelInvEnd;
      /////
      //////*
      ///////          NumVec<MYSCALAR> diagBcast;
      ///////          PMatrix<MYSCALAR> * PMlocBcastPtr;
      ///////          SuperNodeType * superBcastPtr;
      ///////          GridType * g2Ptr;
      ///////
      ///////          if(doSinv_Bcast)
      ///////          {
      ///////            GetTime( timeTotalSelInvSta );
      ///////
      ///////            g2Ptr = new GridType( world_comm, nprow, npcol );
      ///////            GridType &g2 = *g2Ptr;
      ///////
      ///////            superBcastPtr = new SuperNodeType();
      ///////            SuperNodeType & superBcast = *superBcastPtr;
      ///////
      ///////            GetTime( timeSta );
      ///////            luMat.SymbolicToSuperNode( superBcast );
      ///////
      ///////            
      ///////
      ///////            PMlocBcastPtr = new PMatrix( &g2, &superBcast, &luOpt  );
      ///////            PMatrix & PMlocBcast = *PMlocBcastPtr;
      ///////
      ///////            luMat.LUstructToPMatrix( PMlocBcast );
      ///////            GetTime( timeEnd );
      ///////
      ///////            LongInt nnzLU = PMlocBcast.Nnz();
      ///////            if( mpirank == 0 ){
      ///////              cout << "nonzero in L+U  (PMatrix format) = " << nnzLU << endl;
      ///////            }
      ///////
      ///////
      ///////
      ///////            if( mpirank == 0 )
      ///////              cout << "Time for converting LUstruct to PMatrix (Collectives) is " << timeEnd  - timeSta << endl;
      ///////
      ///////            if(doConstructPattern){
      ///////              GetTime( timeSta );
      ///////              PMlocBcast.ConstructCommunicationPattern_Collectives();
      ///////              GetTime( timeEnd );
      ///////              if( mpirank == 0 )
      ///////                cout << "Time for constructing the communication pattern (Collectives) is " << timeEnd  - timeSta << endl;
      ///////
      ///////              if(doPreSelinv){
      ///////                GetTime( timeSta );
      ///////                PMlocBcast.PreSelInv();
      ///////                GetTime( timeEnd );
      ///////                if( mpirank == 0 )
      ///////                  cout << "Time for pre selected inversion (Collectives) is " << timeEnd  - timeSta << endl;
      ///////
      ///////                if(doSelinv){
      ///////                  GetTime( timeSta );
      ///////                  PMlocBcast.SelInv_Collectives();
      ///////                  GetTime( timeEnd );
      ///////                  if( mpirank == 0 )
      ///////                    cout << "Time for numerical selected inversion (Collectives) is " << timeEnd  - timeSta << endl;
      ///////
      ///////                  GetTime( timeTotalSelInvEnd );
      ///////                  if( mpirank == 0 )
      ///////                    cout << "Time for total selected inversion (Collectives) is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;
      ///////
      ///////
      ///////                  // Output the diagonal elements
      ///////                  if( doDiag ){
      ///////                    GetTime( timeSta );
      ///////                    PMlocBcast.GetDiagonal( diagBcast );
      ///////                    GetTime( timeEnd );
      ///////
      ///////
      ///////                    if( mpirank == 0 ){
      ///////                      statusOFS << std::endl << "Diagonal (Collectives) of inverse in natural order: " << std::endl << diagBcast << std::endl;
      ///////                      ofstream ofs("diag_bcast");
      ///////                      if( !ofs.good() ) 
      ///////                        throw std::runtime_error("file cannot be opened.");
      ///////                      serialize( diagBcast, ofs, NO_MASK );
      ///////                      ofs.close();
      ///////                    }
      ///////                  }
      ///////                }
      ///////
      ///////              }
      ///////
      ///////
      ///////
      /////////              PMlocBcast.DestructCommunicators_Collectives( );
      ///////            }
      ///////
      ///////          }
      ///////
      /////*/
      /////
      /////
      /////
      //////*
      ///////          NumVec<MYSCALAR> diagHybrid;
      ///////          PMatrix<MYSCALAR> * PMlocHybridPtr;
      ///////          SuperNodeType * superHybridPtr;
      ///////          GridType * gHybridPtr;
      ///////
      ///////          if(doSinv_Hybrid)
      ///////          {
      ///////            GetTime( timeTotalSelInvSta );
      ///////
      ///////            gHybridPtr = new GridType( world_comm, nprow, npcol );
      ///////            GridType &gHybrid = *gHybridPtr;
      ///////
      ///////            superHybridPtr = new SuperNodeType();
      ///////            SuperNodeType & superHybrid = *superHybridPtr;
      ///////
      ///////            GetTime( timeSta );
      ///////            luMat.SymbolicToSuperNode( superHybrid );
      ///////
      ///////            
      ///////
      ///////            PMlocHybridPtr = new PMatrix( &gHybrid, &superHybrid, &luOpt  );
      ///////            PMatrix & PMlocHybrid = *PMlocHybridPtr;
      ///////
      ///////            luMat.LUstructToPMatrix( PMlocHybrid );
      ///////            GetTime( timeEnd );
      ///////
      ///////            LongInt nnzLU = PMlocHybrid.Nnz();
      ///////            if( mpirank == 0 ){
      ///////              cout << "nonzero in L+U  (PMatrix format) = " << nnzLU << endl;
      ///////            }
      ///////
      ///////
      ///////
      ///////            if( mpirank == 0 )
      ///////              cout << "Time for converting LUstruct to PMatrix (Hybrid) is " << timeEnd  - timeSta << endl;
      ///////
      ///////            if(doConstructPattern){
      ///////              GetTime( timeSta );
      ///////              PMlocHybrid.ConstructCommunicationPattern_Collectives();
      ///////              GetTime( timeEnd );
      ///////              if( mpirank == 0 )
      ///////                cout << "Time for constructing the communication pattern (Hybrid) is " << timeEnd  - timeSta << endl;
      ///////
      ///////              if(doPreSelinv){
      ///////                GetTime( timeSta );
      ///////                PMlocHybrid.PreSelInv();
      ///////                GetTime( timeEnd );
      ///////                if( mpirank == 0 )
      ///////                  cout << "Time for pre selected inversion (Hybrid) is " << timeEnd  - timeSta << endl;
      ///////
      ///////                if(doSelinv){
      ///////                  GetTime( timeSta );
      ///////                  PMlocHybrid.SelInv_Hybrid(doSinv_Hybrid);
      ///////                  GetTime( timeEnd );
      ///////                  if( mpirank == 0 )
      ///////                    cout << "Time for numerical selected inversion (Hybrid) is " << timeEnd  - timeSta << endl;
      ///////
      ///////                  GetTime( timeTotalSelInvEnd );
      ///////                  if( mpirank == 0 )
      ///////                    cout << "Time for total selected inversion (Hybrid) is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;
      ///////
      ///////
      ///////                  // Output the diagonal elements
      ///////                  if( doDiag ){
      ///////                    GetTime( timeSta );
      ///////                    PMlocHybrid.GetDiagonal( diagHybrid );
      ///////                    GetTime( timeEnd );
      ///////
      ///////
      ///////                    if( mpirank == 0 ){
      ///////                      statusOFS << std::endl << "Diagonal (Hybrid) of inverse in natural order: " << std::endl << diagHybrid << std::endl;
      ///////                      ofstream ofs("diag_hybrid");
      ///////                      if( !ofs.good() ) 
      ///////                        throw std::runtime_error("file cannot be opened.");
      ///////                      serialize( diagHybrid, ofs, NO_MASK );
      ///////                      ofs.close();
      ///////                    }
      ///////                  }
      ///////                }
      ///////
      ///////              }
      ///////            }
      ///////
      ///////          }
      /////*/
      /////
      /////
      /////          NumVec<MYSCALAR> diag;
      /////          PMatrix<MYSCALAR> * PMlocPtr;
      /////          SuperNodeType * superPtr;
      /////          GridType * g1Ptr;
      /////
      /////          if(doSinvPipeline)
      /////          {
      /////            GetTime( timeTotalSelInvSta );
      /////
      /////            g1Ptr = new GridType( world_comm, nprow, npcol );
      /////            GridType &g1 = *g1Ptr;
      /////
      /////            superPtr = new SuperNodeType();
      /////            SuperNodeType & super = *superPtr;
      /////
      /////            GetTime( timeSta );
      /////            luMat.SymbolicToSuperNode( super );
      /////
      /////
      /////
      /////            PMlocPtr = new PMatrix<MYSCALAR>( &g1, &super, &luOpt  );
      /////
      /////
      /////
      /////
      /////            PMatrix<MYSCALAR> & PMloc = *PMlocPtr;
      /////
      /////            luMat.LUstructToPMatrix( PMloc );
      /////            GetTime( timeEnd );
      /////
      /////            LongInt nnzLU = PMloc.Nnz();
      /////            if( mpirank == 0 ){
      /////              cout << "nonzero in L+U  (PMatrix format) = " << nnzLU << endl;
      /////            }
      /////
      /////
      /////            if( mpirank == 0 )
      /////              cout << "Time for converting LUstruct to PMatrix is " << timeEnd  - timeSta << endl;
      /////
      /////            //        statusOFS << "perm: " << endl << super.perm << endl;
      /////            //        statusOFS << "permInv: " << endl << super.permInv << endl;
      /////            //        statusOFS << "superIdx:" << endl << super.superIdx << endl;
      /////            //        statusOFS << "superPtr:" << endl << super.superPtr << endl; 
      /////
      /////
      /////            // Preparation for the selected inversion
      /////            GetTime( timeSta );
      /////            PMloc.ConstructCommunicationPattern();
      ///////            PMloc.ConstructCommunicationPattern_Collectives();
      /////            GetTime( timeEnd );
      /////
      /////            if( mpirank == 0 )
      /////              cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;
      /////
      /////
      /////            GetTime( timeSta );
      /////            PMloc.PreSelInv();
      /////            GetTime( timeEnd );
      /////
      /////            if( mpirank == 0 )
      /////              cout << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;
      /////
      /////            // Main subroutine for selected inversion
      /////            GetTime( timeSta );
      /////            PMloc.SelInv();
      /////            GetTime( timeEnd );
      /////            if( mpirank == 0 )
      /////              cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;
      /////
      /////
      /////            GetTime( timeTotalSelInvEnd );
      /////            if( mpirank == 0 )
      /////              cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;
      /////
      /////
      ///// 
      /////             if(doToDist){
      /////               // Convert to DistSparseMatrix and get the diagonal
      /////               GetTime( timeSta );
      /////               DistSparseMatrix<MYSCALAR> Ainv;
      /////               PMloc.PMatrixToDistSparseMatrix( Ainv );
      /////               GetTime( timeEnd );
      ///// 
      /////               if( mpirank == 0 )
      /////                 cout << "Time for converting PMatrix to DistSparseMatrix is " << timeEnd  - timeSta << endl;
      ///// 
      /////               NumVec<MYSCALAR> diagDistSparse;
      /////               GetTime( timeSta );
      /////               GetDiagonal( Ainv, diagDistSparse );
      /////               GetTime( timeEnd );
      /////               if( mpirank == 0 )
      /////                 cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
      ///// 
      /////               if( mpirank == 0 ){
      /////                 statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format : " << std::endl << diagDistSparse << std::endl;
      /////                 Real diffNorm = 0.0;;
      /////                 for( Int i = 0; i < diag.m(); i++ ){
      /////                   diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
      /////                 }
      /////                 diffNorm = std::sqrt( diffNorm );
      /////                 statusOFS << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
      /////               }
      ///// 
      /////               // Convert to DistSparseMatrix in the 2nd format and get the diagonal
      /////               GetTime( timeSta );
      /////               DistSparseMatrix<MYSCALAR> Ainv2;
      /////               PMloc.PMatrixToDistSparseMatrix2( AMat, Ainv2 );
      /////               GetTime( timeEnd );
      ///// 
      /////               if( mpirank == 0 )
      /////                 cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;
      ///// 
      /////               NumVec<MYSCALAR> diagDistSparse2;
      /////               GetTime( timeSta );
      /////               GetDiagonal( Ainv2, diagDistSparse2 );
      /////               GetTime( timeEnd );
      /////               if( mpirank == 0 )
      /////                 cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
      ///// 
      /////               if( mpirank == 0 ){
      /////                 statusOFS << std::endl << "Diagonal of inverse from the 2nd conversion into DistSparseMatrix format : " << std::endl << diagDistSparse2 << std::endl;
      /////                 Real diffNorm = 0.0;;
      /////                 for( Int i = 0; i < diag.m(); i++ ){
      /////                   diffNorm += pow( std::abs( diag(i) - diagDistSparse2(i) ), 2.0 );
      /////                 }
      /////                 diffNorm = std::sqrt( diffNorm );
      /////                 statusOFS << std::endl << "||diag - diagDistSparse2||_2 = " << diffNorm << std::endl;
      /////               }
      ///// 
      /////               Complex traceLocal = blas::Dotu( AMat.nnzLocal, AMat.nzvalLocal.Data(), 1,
      /////                   Ainv2.nzvalLocal.Data(), 1 );
      /////               Complex trace = Z_ZERO;
      /////               mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, world_comm );
      ///// 
      /////               if( mpirank == 0 ){
      ///// 
      /////                 cout << "H.size = "  << HMat.size << endl;
      /////                 cout << std::endl << "Tr[Ainv2 * AMat] = " <<  trace << std::endl;
      /////                 statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;
      ///// 
      /////                 cout << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( Complex(HMat.size, 0.0) - trace ) << std::endl;
      /////                 statusOFS << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( Complex(HMat.size, 0.0) - trace ) << std::endl;
      ///// 
      /////               }
      /////             }
      /////
      /////
      ///////            GetTime( timeSta );
      ///////            luMat.LUstructToPMatrix( PMloc );
      ///////            GetTime( timeEnd );
      ///////
      ///////            nnzLU = PMloc.Nnz();
      ///////            if( mpirank == 0 ){
      ///////              cout << "nonzero in L+U  (PMatrix format) = " << nnzLU << endl;
      ///////            }
      ///////
      ///////
      ///////            if( mpirank == 0 )
      ///////              cout << "Time for converting LUstruct to PMatrix is " << timeEnd  - timeSta << endl;
      ///////
      ///////
      ///////
      ///////             GetTime( timeSta );
      ///////             PMloc.PreSelInv();
      ///////             GetTime( timeEnd );
      /////// 
      ///////             if( mpirank == 0 )
      ///////               cout << "Time for pre-selected inversion 2 is " << timeEnd  - timeSta << endl;
      /////// 
      ///////             // Main subroutine for selected inversion
      ///////             GetTime( timeSta );
      ///////             PMloc.SelInv();
      ///////             GetTime( timeEnd );
      ///////             if( mpirank == 0 )
      ///////               cout << "Time for numerical selected inversion 2 is " << timeEnd  - timeSta << endl;
      /////// 
      /////// 
      ///////             GetTime( timeTotalSelInvEnd );
      ///////             if( mpirank == 0 )
      ///////               cout << "Time for total selected inversion 2 is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;
      /////// 
      /////// 
      /////// 
      /////// 
      /////// 
      ///////             if(doToDist){
      ///////               // Convert to DistSparseMatrix and get the diagonal
      ///////               GetTime( timeSta );
      ///////               DistSparseMatrix<MYSCALAR> Ainv;
      ///////               PMloc.PMatrixToDistSparseMatrix( Ainv );
      ///////               GetTime( timeEnd );
      /////// 
      ///////               if( mpirank == 0 )
      ///////                 cout << "Time for converting PMatrix to DistSparseMatrix is " << timeEnd  - timeSta << endl;
      /////// 
      ///////               NumVec<MYSCALAR> diagDistSparse;
      ///////               GetTime( timeSta );
      ///////               GetDiagonal( Ainv, diagDistSparse );
      ///////               GetTime( timeEnd );
      ///////               if( mpirank == 0 )
      ///////                 cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
      /////// 
      ///////               if( mpirank == 0 ){
      ///////                 statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format : " << std::endl << diagDistSparse << std::endl;
      ///////                 Real diffNorm = 0.0;;
      ///////                 for( Int i = 0; i < diag.m(); i++ ){
      ///////                   diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
      ///////                 }
      ///////                 diffNorm = std::sqrt( diffNorm );
      ///////                 statusOFS << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
      ///////               }
      /////// 
      ///////               // Convert to DistSparseMatrix in the 2nd format and get the diagonal
      ///////               GetTime( timeSta );
      ///////               DistSparseMatrix<MYSCALAR> Ainv2;
      ///////               PMloc.PMatrixToDistSparseMatrix2( AMat, Ainv2 );
      ///////               GetTime( timeEnd );
      /////// 
      ///////               if( mpirank == 0 )
      ///////                 cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;
      /////// 
      ///////               NumVec<MYSCALAR> diagDistSparse2;
      ///////               GetTime( timeSta );
      ///////               GetDiagonal( Ainv2, diagDistSparse2 );
      ///////               GetTime( timeEnd );
      ///////               if( mpirank == 0 )
      ///////                 cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
      /////// 
      ///////               if( mpirank == 0 ){
      ///////                 statusOFS << std::endl << "Diagonal of inverse from the 2nd conversion into DistSparseMatrix format : " << std::endl << diagDistSparse2 << std::endl;
      ///////                 Real diffNorm = 0.0;;
      ///////                 for( Int i = 0; i < diag.m(); i++ ){
      ///////                   diffNorm += pow( std::abs( diag(i) - diagDistSparse2(i) ), 2.0 );
      ///////                 }
      ///////                 diffNorm = std::sqrt( diffNorm );
      ///////                 statusOFS << std::endl << "||diag - diagDistSparse2||_2 = " << diffNorm << std::endl;
      ///////               }
      ///////
      ///////
      ///////                //dump the last supernode
      ///////                 statusOFS << "Ainv2 = "  << Ainv2.nzvalLocal << endl;
      ///////
      ///////               Complex traceLocal = blas::Dotu( AMat.nnzLocal, AMat.nzvalLocal.Data(), 1,
      ///////                   Ainv2.nzvalLocal.Data(), 1 );
      ///////               Complex trace = Z_ZERO;
      ///////               mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, world_comm );
      /////// 
      ///////               if( mpirank == 0 ){
      /////// 
      ///////                 cout << "H.size = "  << HMat.size << endl;
      ///////                 cout << std::endl << "Tr[Ainv2 * AMat] = " <<  trace << std::endl;
      ///////                 statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;
      /////// 
      ///////                 cout << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( Complex(HMat.size, 0.0) - trace ) << std::endl;
      ///////                 statusOFS << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( Complex(HMat.size, 0.0) - trace ) << std::endl;
      /////// 
      ///////               }
      ///////             }
      ///// 
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////            // Output the diagonal elements
      /////            if( doDiag ){
      /////              NumVec<MYSCALAR> diag;
      /////
      /////              GetTime( timeSta );
      /////              PMloc.GetDiagonal( diag );
      /////              GetTime( timeEnd );
      /////
      /////
      /////              if( mpirank == 0 )
      /////                cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;
      /////
      /////
      /////              if( mpirank == 0 ){
      /////                statusOFS << std::endl << "Diagonal (pipeline) of inverse in natural order: " << std::endl << diag << std::endl;
      /////                ofstream ofs("diag");
      /////                if( !ofs.good() ) 
      /////                  throw std::runtime_error("file cannot be opened.");
      /////                serialize( diag, ofs, NO_MASK );
      /////                ofs.close();
      /////              }
      /////            }
      /////
      /////
      /////
      /////          }
      /////
      /////
      /////
      /////          if(doSinvPipeline || doSinv_Bcast || doSinv_Hybrid){
      /////
      ///////            PMatrix<MYSCALAR> * PMloc = doSinvPipeline?PMlocPtr:(doSinv_Bcast?PMlocBcastPtr:PMlocHybridPtr);
      /////            PMatrix<MYSCALAR> * PMloc = PMlocPtr;
      /////
      /////            if(doToDist){
      /////              // Convert to DistSparseMatrix and get the diagonal
      /////              GetTime( timeSta );
      /////              DistSparseMatrix<MYSCALAR> Ainv;
      /////              PMloc->PMatrixToDistSparseMatrix( Ainv );
      /////              GetTime( timeEnd );
      /////
      /////              if( mpirank == 0 )
      /////                cout << "Time for converting PMatrix to DistSparseMatrix is " << timeEnd  - timeSta << endl;
      /////
      /////              NumVec<MYSCALAR> diagDistSparse;
      /////              GetTime( timeSta );
      /////              GetDiagonal( Ainv, diagDistSparse );
      /////              GetTime( timeEnd );
      /////              if( mpirank == 0 )
      /////                cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
      /////
      /////              if( mpirank == 0 ){
      /////                statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format : " << std::endl << diagDistSparse << std::endl;
      /////                Real diffNorm = 0.0;;
      /////                for( Int i = 0; i < diag.m(); i++ ){
      /////                  diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
      /////                }
      /////                diffNorm = std::sqrt( diffNorm );
      /////                statusOFS << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
      /////              }
      /////
      /////              // Convert to DistSparseMatrix in the 2nd format and get the diagonal
      /////              GetTime( timeSta );
      /////              DistSparseMatrix<MYSCALAR> Ainv2;
      /////              PMloc->PMatrixToDistSparseMatrix2( AMat, Ainv2 );
      /////              GetTime( timeEnd );
      /////
      /////              if( mpirank == 0 )
      /////                cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;
      /////
      /////              NumVec<MYSCALAR> diagDistSparse2;
      /////              GetTime( timeSta );
      /////              GetDiagonal( Ainv2, diagDistSparse2 );
      /////              GetTime( timeEnd );
      /////              if( mpirank == 0 )
      /////                cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
      /////
      /////              if( mpirank == 0 ){
      /////                statusOFS << std::endl << "Diagonal of inverse from the 2nd conversion into DistSparseMatrix format : " << std::endl << diagDistSparse2 << std::endl;
      /////                Real diffNorm = 0.0;;
      /////                for( Int i = 0; i < diag.m(); i++ ){
      /////                  diffNorm += pow( std::abs( diag(i) - diagDistSparse2(i) ), 2.0 );
      /////                }
      /////                diffNorm = std::sqrt( diffNorm );
      /////                statusOFS << std::endl << "||diag - diagDistSparse2||_2 = " << diffNorm << std::endl;
      /////              }
      /////
      /////              MYSCALAR traceLocal = blas::Dotu( AMat.nnzLocal, AMat.nzvalLocal.Data(), 1, 
      /////                  Ainv2.nzvalLocal.Data(), 1 );
      /////              MYSCALAR trace = ZERO<MYSCALAR>();
      /////              mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, world_comm );
      /////
      /////              if( mpirank == 0 ){
      /////
      /////                cout << "H.size = "  << HMat.size << endl;
      /////                cout << std::endl << "Tr[Ainv2 * AMat] = " <<  trace << std::endl;
      /////                statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;
      /////
      /////                cout << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( static_cast<MYSCALAR>(HMat.size) - trace ) << std::endl;
      /////                statusOFS << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( static_cast<MYSCALAR>(HMat.size) - trace ) << std::endl;
      /////
      /////              }
      /////            }
      /////
      /////          }
      /////
      /////
      /////
      /////
      //////*
      ///////          if(doSinv_Bcast){
      ///////            delete PMlocBcastPtr;
      ///////            delete superBcastPtr;
      ///////            delete g2Ptr;
      ///////          }
      /////*/
      //////*
      ///////          if(doSinv_Hybrid){
      ///////            delete PMlocHybridPtr;
      ///////            delete superHybridPtr;
      ///////            delete gHybridPtr;
      ///////          }
      /////*/
      /////
      /////
      /////          if(doSinvPipeline){
      /////            delete PMlocPtr;
      /////            delete superPtr;
      /////            delete g1Ptr;
      /////          }
      /////
      /////        }
      /////
      /////
      /////      }


      statusOFS.close();
    }
  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
#ifndef _RELEASE_
    DumpCallStack();
#endif
  }

  MPI_Finalize();

  return 0;
}

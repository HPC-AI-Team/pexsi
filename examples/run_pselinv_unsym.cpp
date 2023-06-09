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

#include "pexsi/timer.h"

#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#else
#define MYSCALAR Real
#endif


using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout << "Usage" << std::endl << "run_pselinv -T [isText] -F [doFacto -E [doTriSolve] -Sinv [doSelInv]]  -H <Hfile> -S [Sfile] -colperm [colperm] -r [nprow] -c [npcol] -npsymbfact [npsymbfact] -P [maxpipelinedepth] -SinvBcast [doSelInvBcast] -SinvPipeline [doSelInvPipeline] -SinvHybrid [doSelInvHybrid] -rshift [real shift] -ishift [imaginary shift] -ToDist [doToDist] -Diag [doDiag] -Real [isReal] -Symm [isSymm] -rowperm [LargeDiag|NOROWPERM] -equil [YES|NO] -transpose [0|1]" << std::endl;
}

int main(int argc, char **argv) 
{

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
          ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        ErrorHandling( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        ErrorHandling( "When using -c option, -r also needs to be provided." );
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

#ifdef GEMM_PROFILE
      stringstream  ss2;
      ss2 << "gemm_stat" << mpirank;
      statOFS.open( ss2.str().c_str());
#endif

#ifdef COMM_PROFILE
      stringstream  ss3;
      ss3 << "comm_stat" << mpirank;
      commOFS.open( ss3.str().c_str());
#endif


      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  ErrorHandling( "nprow == npcol is assumed in this test routine." );
      //}

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

      std::string Hfile, Sfile;
      int isCSC = true;
      if( options.find("-T") != options.end() ){ 
        isCSC= ! atoi(options["-T"].c_str());
      }


      int isSym = 0;
      if( options.find("-Symm") != options.end() ){ 
        isSym = atoi(options["-Symm"].c_str());
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
        ErrorHandling("Hfile must be provided.");
      }

      if( options.find("-S") != options.end() ){ 
        Sfile = options["-S"];
      }
      else{
        statusOFS << "-S option is not given. " 
          << "Treat the overlap matrix as an identity matrix." 
          << std::endl << std::endl;
      }

      Int maxPipelineDepth = -1;
      if( options.find("-P") != options.end() ){ 
        maxPipelineDepth = atoi(options["-P"].c_str());
      }
      else{
        statusOFS << "-P option is not given. " 
          << "Do not limit SelInv pipelining depth." 
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

      Int isReal = 0;
      if( options.find("-Real") != options.end() ){ 
        isReal = atoi(options["-Real"].c_str());
      }

      Int doSelinv = 1;
      if( options.find("-Selinv") != options.end() ){ 
        doSelinv = atoi(options["-Selinv"].c_str());
      }

      Real rshift = 0.0, ishift = 0.0;
      if( options.find("-rshift") != options.end() ){ 
        rshift = atof(options["-rshift"].c_str());
      }
      if( options.find("-ishift") != options.end() ){ 
        ishift = atof(options["-ishift"].c_str());
      }


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



      std::string RowPerm;
      if( options.find("-rowperm") != options.end() ){ 
        RowPerm = options["-rowperm"];
      }
      else{
        statusOFS << "-rowperm option is not given. " 
          << "Use NOROWPERM." 
          << std::endl << std::endl;
        RowPerm = "NOROWPERM";
      }

      std::string Equil;
      if( options.find("-equil") != options.end() ){ 
        Equil = options["-equil"];
      }
      else{
        statusOFS << "-equil option is not given. " 
          << "Use NO." 
          << std::endl << std::endl;
        Equil = "NO";
      }


      Int transpose;
      if( options.find("-transpose") != options.end() ){ 
        transpose = atoi( options["-transpose"].c_str() );
      }
      else{
        statusOFS << "-transpose option is not given. " 
          << "Use default value (0)." 
          << std::endl << std::endl;
        transpose = 0;
      }

      // *********************************************************************
      // Read input matrix
      // *********************************************************************

      // Setup grid.
      SuperLUGrid<MYSCALAR> g( world_comm, nprow, npcol );

      int      m, n;
      DistSparseMatrix<MYSCALAR>  AMat;
      Real timeSta, timeEnd;


#ifdef _MYCOMPLEX_
      if(!isReal){
        DistSparseMatrix<Complex> HMat;
        DistSparseMatrix<Complex> SMat;
        GetTime( timeSta );
        if(isCSC){
          ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
        }
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
          if(isCSC){
            ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
          }
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
        Complex *ptr1 = HMat.nzvalLocal.Data();
        Complex *ptr2 = SMat.nzvalLocal.Data();

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
      }
      else
#endif
      {
        DistSparseMatrix<Real> HMat;
        DistSparseMatrix<Real> SMat;
        Real timeSta, timeEnd;
        GetTime( timeSta );
        if(isCSC){
          ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
        }
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
          if(isCSC){
            ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
          }
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
      }

      LongInt nnzA = AMat.Nnz();
      if( mpirank == 0 ){
        cout << "nonzero in A (DistSparseMatrix format) = " << nnzA << endl;
      }


      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;



      GetTime( timeSta );
      SuperLUOptions luOpt;
      luOpt.ColPerm = ColPerm;
      luOpt.RowPerm = RowPerm;
      luOpt.Equil = Equil;
      luOpt.numProcSymbFact = numProcSymbFact;
      luOpt.Symmetric = isSym;
      luOpt.Transpose = transpose;

      FactorizationOptions factOpt;
      factOpt.ColPerm = ColPerm;
      factOpt.RowPerm = RowPerm;
      factOpt.Symmetric = isSym;


      //Initialize SuperLU data structures
      //SuperLUGrid<MYSCALAR> * pLuGrid;
      //SuperLUMatrix<MYSCALAR> * pLuMat;
      //SuperNodeType * pSuper;
      //PEXSICreator<MYSCALAR>::CreateSuperLUMatrix(world_comm, nprow, npcol, luOpt, pLuMat, pLuGrid, pSuper);

      SuperLUGrid<MYSCALAR> * pLuGrid = new SuperLUGrid<MYSCALAR>(world_comm,nprow,npcol);
      SuperLUMatrix<MYSCALAR> * pLuMat = new SuperLUMatrix<MYSCALAR>(*pLuGrid, luOpt);
      SuperNodeType * pSuper = new SuperNodeType();

      // *********************************************************************
      // Symbolic factorization 
      // *********************************************************************

      pLuMat->DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );
      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;


      if(doSymbfact){
        GetTime( timeSta );
        pLuMat->SymbolicFactorize();
        pLuMat->DestroyAOnly();
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for performing the symbolic factorization is " << timeEnd - timeSta << endl;

      }

      // *********************************************************************
      // Numerical factorization only 
      // *********************************************************************

      if(doFacto){
        Real timeTotalFactorizationSta, timeTotalFactorizationEnd; 


        // Important: the distribution in pzsymbfact is going to mess up the
        // A matrix.  Recompute the matrix A here.
        pLuMat->DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );

        GetTime( timeTotalFactorizationSta );

        GetTime( timeSta );
        pLuMat->Distribute();
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl; 

        //        {
        //          SuperNodeType * pSuper = new SuperNodeType();
        //          pLuMat->SymbolicToSuperNode( *pSuper );
        //          GridType * pGrid = new GridType(world_comm,nprow,npcol);
        //          PMatrix<MYSCALAR> * pMat = PMatrix<MYSCALAR>::Create(pGrid,pSuper, &luOpt);
        //          pLuMat->LUstructToPMatrix( *pMat );
        //          pMat->DumpLU();
        //          delete pMat;
        //          delete pGrid;
        //          delete pSuper;
        //        }

        GetTime( timeSta );
        pLuMat->NumericalFactorize();
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 

        GetTime( timeTotalFactorizationEnd );
        if( mpirank == 0 )
          cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 


        // *********************************************************************
        // Test the accuracy of factorization by solve
        // *********************************************************************

        if( checkAccuracy ) {
          SuperLUMatrix<MYSCALAR> A1( *pLuGrid, luOpt );
          SuperLUMatrix<MYSCALAR> GA( *pLuGrid, luOpt );

          A1.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );
          A1.ConvertNRlocToNC( GA );

          int n = A1.n();
          int nrhs = 1;
          NumMat<MYSCALAR> xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
          NumMat<MYSCALAR> xTrueLocal, bLocal;
          DblNumVec berr;
          UniformRandom( xTrueGlobal );

          GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );

          A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
          A1.DistributeGlobalMultiVector( bGlobal,     bLocal );

          if(mpirank==0){std::cout<<"Starting solve"<<std::endl;}
          pLuMat->SolveDistMultiVector( bLocal, berr );
          pLuMat->CheckErrorDistMultiVector( bLocal, xTrueLocal );
        }

        // *********************************************************************
        // Selected inversion
        // *********************************************************************

        for(int i=1; i<= doSelInv; ++i )
        {

          Real timeTotalSelInvSta, timeTotalSelInvEnd;
          NumVec<MYSCALAR> diag;


          GetTime( timeSta );
          pLuMat->SymbolicToSuperNode( *pSuper );

          if(1){
            statusOFS<<"superIdx: "<<pSuper->superIdx<<std::endl;
            statusOFS<<"superPtr: "<<pSuper->superPtr<<std::endl;
            statusOFS<<"colperm: "<<pSuper->perm<<std::endl;
            statusOFS<<"rowperm: "<<pSuper->perm_r<<std::endl;
          }

          GetTime( timeTotalSelInvSta );


          //Initialize PEXSI/PSelInv data structures
          GridType * pGrid = new GridType(world_comm,nprow,npcol);

          PSelInvOptions selInvOpt;
          selInvOpt.maxPipelineDepth = maxPipelineDepth;

          PMatrix<MYSCALAR> * pMat = PMatrix<MYSCALAR>::Create(pGrid,pSuper, &selInvOpt, &factOpt);

          pLuMat->LUstructToPMatrix( *pMat );
          GetTime( timeEnd );

          LongInt nnzLU = pMat->Nnz();
          if( mpirank == 0 ){
            cout << "nonzero in L+U  (PMatrix format) = " << nnzLU << endl;
          }


          if( mpirank == 0 )
            cout << "Time for converting LUstruct to PMatrix is " << timeEnd  - timeSta << endl;

          if(0){
            statusOFS.close();
            stringstream  sslu;
            sslu << "LUfactDump_" << mpirank<<".m";
            statusOFS.open( sslu.str().c_str() );
            pMat->DumpLU();
            statusOFS.close();
            statusOFS.open( ss.str().c_str(),std::ofstream::out | std::ofstream::app );
          }

          // Preparation for the selected inversion
          GetTime( timeSta );
          pMat->ConstructCommunicationPattern();
          GetTime( timeEnd );

          if( mpirank == 0 )
            cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;


          GetTime( timeSta );
          pMat->PreSelInv();
          GetTime( timeEnd );

          if( mpirank == 0 )
            cout << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;


          if(0){
            statusOFS.close();
            stringstream  sslu;
            sslu << "PreLUDump_" << mpirank<<".m";
            statusOFS.open( sslu.str().c_str() );
            pMat->DumpLU();
            statusOFS.close();
            statusOFS.open( ss.str().c_str(),std::ofstream::out | std::ofstream::app );
          }


          // Main subroutine for selected inversion
          GetTime( timeSta );
          pMat->SelInv();
          GetTime( timeEnd );
          if( mpirank == 0 )
            cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;


          GetTime( timeTotalSelInvEnd );
          if( mpirank == 0 )
            cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;

#ifdef _PRINT_STATS_
          double flops = pMat->GetTotalFlops();
          if( mpirank == 0 )
            cout << "Total FLOPs for selected inversion is " << flops << endl;
#endif

          if(0){
            statusOFS.close();
            stringstream  sslu;
            sslu << "LUDump_" << mpirank<<".m";
            statusOFS.open( sslu.str().c_str() );
            pMat->DumpLU();
            statusOFS.close();
            statusOFS.open( ss.str().c_str(),std::ofstream::out | std::ofstream::app );
          }

          if(1){
            DistSparseMatrix<MYSCALAR> Ainv;
            //pMat->PMatrixToDistSparseMatrix(Ainv );
            //pMat->PMatrixToDistSparseMatrix( AMat, Ainv );
            //if( mpirank == 0 )
            //WriteDistSparseMatrixMatlab("AinvDump",Ainv,world_comm);

            //WriteDistSparseMatrixMatlab("ADump",AMat,world_comm);

            pMat->PMatrixToDistSparseMatrix( AMat, Ainv );
            WriteDistSparseMatrixMatlab("AinvDump",Ainv,world_comm);
          }

          if(doToDist){
            // Convert to DistSparseMatrix in the 2nd format and get the diagonal
            DistSparseMatrix<MYSCALAR> Ainv;
            MYSCALAR traceLocal;


            DistSparseMatrix<MYSCALAR> * Aptr;
#ifndef OLD_SELINV
            Aptr = &AMat;
#else
            if(factOpt.Symmetric==0 && 0){
              Aptr = new DistSparseMatrix<MYSCALAR>();
              //compute the transpose
              CSCToCSR(AMat,*Aptr);
            }
            else{
              Aptr = &AMat;
            }
#endif
            GetTime( timeSta );
            pMat->PMatrixToDistSparseMatrix( *Aptr, Ainv );
            // Note: At this point, the matrix Ainv actually stores the
            // selected elements of A^{-T}. 
            // According to the rules of selected inversion, the set of selected
            // elements being computed are
            // A^{-T}_{ij} such that A_{ij} \ne 0.
            // 
            // This is why the computation of the trace is computed as
            //
            // Tr[AA^{-1}] = \sum_{ij} A_{ij} A^{-1}_{ji}
            //             = \sum_{ij} A_{ij} A^{-T}{ij},
            // 
            // and is implemented using dotu below
            GetTime( timeEnd );

#ifndef OLD_SELINV
            traceLocal = ZERO<MYSCALAR>();
            traceLocal = blas::Dotu( Aptr->nnzLocal, Aptr->nzvalLocal.Data(), 1,
                Ainv.nzvalLocal.Data(), 1 );
#else
            traceLocal = ZERO<MYSCALAR>();
            traceLocal = blas::Dotu( Aptr->nnzLocal, Ainv.nzvalLocal.Data(), 1,
                Aptr->nzvalLocal.Data(), 1 );
#endif

#ifdef OLD_SELINV
            if(factOpt.Symmetric==0 && 0){//factOpt.Transpose==0){
              delete Aptr;
            }
#endif


            if( mpirank == 0 )
              cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;

            MYSCALAR trace = ZERO<MYSCALAR>();
            mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, world_comm );

            if( mpirank == 0 ){

              cout << "H.size = "  << AMat.size << endl;
              cout << std::endl << "Tr[Ainv * AMat] = " <<  trace << std::endl;
              statusOFS << std::endl << "Tr[Ainv * AMat] = " << std::endl << trace << std::endl;
#ifdef _MYCOMPLEX_ 
              cout << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( Complex(AMat.size, 0.0) - trace ) << std::endl;
              statusOFS << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( Complex(AMat.size, 0.0) - trace ) << std::endl;
#else
              cout << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( AMat.size - trace ) << std::endl;
              statusOFS << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( AMat.size - trace ) << std::endl;
#endif
            }

            if( doDiag ){
              NumVec<MYSCALAR> diag;

              GetTime( timeSta );
              pMat->GetDiagonal( diag );
              GetTime( timeEnd );


              if( mpirank == 0 )
                cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;


              NumVec<MYSCALAR> diagDistSparse;
              GetTime( timeSta );
#ifndef OLD_SELINV
              GetDiagonal( Ainv, diagDistSparse );
#else
              GetDiagonal( Ainv, diagDistSparse );
#endif
              GetTime( timeEnd );
              if( mpirank == 0 )
                cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;

              if( mpirank == 0 ){
                statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format: " << std::endl << diagDistSparse << std::endl;
                Real diffNorm = 0.0;;
                for( Int i = 0; i < diag.m(); i++ ){
                  diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
                }
                diffNorm = std::sqrt( diffNorm );
                cout << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
              }


              if( mpirank == 0 ){
                statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;
                ofstream ofs("diag");
                if( !ofs.good() ) 
                  ErrorHandling("file cannot be opened.");
                serialize( diag, ofs, NO_MASK );
                ofs.close();
              }
            }


            }
            else if( doDiag ){
              NumVec<MYSCALAR> diag;

              GetTime( timeSta );
              pMat->GetDiagonal( diag );
              GetTime( timeEnd );


              if( mpirank == 0 )
                cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;

              if( mpirank == 0 ){
                statusOFS << std::endl << "Diagonal (pipeline) of inverse in natural order: " << std::endl << diag << std::endl;
                ofstream ofs("diag");
                if( !ofs.good() ) 
                  ErrorHandling("file cannot be opened.");
                serialize( diag, ofs, NO_MASK );
                ofs.close();
              }
            }


            delete pMat;
            delete pGrid;

            }


          }

          delete pSuper;
          delete pLuMat;
          delete pLuGrid;

          statusOFS.close();
#ifdef GEMM_PROFILE
          statOFS.close();
#endif

#ifdef COMM_PROFILE
          commOFS.close();
#endif


        }
      }
      catch( std::exception& e )
      {
        std::cerr << "Processor " << mpirank << " caught exception with message: "
          << e.what() << std::endl;
      }

      MPI_Finalize();

      return 0;
    }

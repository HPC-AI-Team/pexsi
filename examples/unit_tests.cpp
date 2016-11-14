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
/// @file unit_tests.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @date 2013-04-15
#include  "ppexsi.hpp"
#include "pexsi/timer.h"

//#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#else
#define MYSCALAR Real
#endif


using namespace PEXSI;
using namespace std;



void check_broadcasts(MPI_Comm world_comm,
      const std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > > & bcastTrees, 
            const NumMat<MYSCALAR> & sbuf, const NumMat<MYSCALAR> &rbuf){
  //validate the output
  int mpisize = 0;
  int mpirank = 0;
      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );

  NumMat<MYSCALAR> checkBuf(sbuf.m(),sbuf.n());
  SetValue(checkBuf,MYSCALAR(0.0));

  for(Int t = 0; t<bcastTrees.size(); t++){
    Int root = t%2;
    if(t == bcastTrees.size()-1){
      root = 0;
    }

    MPI_Bcast(mpirank==root?sbuf.Data():checkBuf.Data(),sbuf.ByteSize(),MPI_BYTE,root,world_comm);

    auto & bcastTree =  bcastTrees[t];
    if(bcastTree!=nullptr){
      if(!bcastTree->IsRoot()){
        for(Int i = 0; i<checkBuf.m(); i++){
          for(Int j = 0; j<checkBuf.n(); j++){
            assert(checkBuf(i,j)==rbuf(i,j));
          }
        }
      }
    }
  }
  statusOFS<<"Broadcast trees test passed"<<std::endl;
}


void check_reduces(MPI_Comm world_comm,
      const std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > > & redTrees, 
      const NumMat<MYSCALAR> & sbuf, const NumMat<MYSCALAR> &rbuf){
  //validate the output
  int mpisize = 0;
  int mpirank = 0;
  MPI_Comm_rank(world_comm, &mpirank );
  MPI_Comm_size(world_comm, &mpisize );

  NumMat<MYSCALAR> checkBuf(sbuf.m(),sbuf.n());
  NumMat<MYSCALAR> tmpBuf(sbuf.m(),sbuf.n());
  SetValue(checkBuf,MYSCALAR(0.0));

  for(Int t = 0; t<redTrees.size(); t++){
    std::vector<Int> ranks;
    if(t==redTrees.size()-1){
      ranks.push_back(0);
      for(Int i = 1; i<mpisize; i++){
        ranks.push_back(i);
      }
    }
    else{
      ranks.push_back(t%2);
      for(Int i = t%2+2; i<mpisize; i+=2){
        ranks.push_back(i);
      }
    }
    Int root = ranks[0];

    for(auto p : ranks){
      if(mpirank==root){
        if(mpirank!=p){
          MPI_Recv(tmpBuf.Data(),tmpBuf.ByteSize(),MPI_BYTE,MPI_ANY_SOURCE,t,world_comm,MPI_STATUS_IGNORE);
          //axpy
          blas::Axpy(tmpBuf.Size(), ONE<MYSCALAR>(), tmpBuf.Data(), 1, checkBuf.Data(), 1 );
        }
//        else{
//          blas::Axpy(rbuf.Size(), ONE<MYSCALAR>(), rbuf.Data(), 1, checkBuf.Data(), 1 );
//        }
      }
      else if(mpirank == p){
        MPI_Send(sbuf.Data(),sbuf.ByteSize(),MPI_BYTE,root,t,world_comm,MPI_STATUS_IGNORE);
      }
    }

    auto & redTree =  redTrees[t];
    if(redTree!=nullptr){
      if(redTree->IsRoot()){
        for(Int i = 0; i<checkBuf.m(); i++){
          for(Int j = 0; j<checkBuf.n(); j++){
            assert(checkBuf(i,j)==rbuf(i,j));
          }
        }
      }
    }

  }
  statusOFS<<"Reduce trees test passed"<<std::endl;
}





int main(int argc, char **argv) 
{

  if( argc < 3 ) {
    return 0;
  }


  MPI_Init( &argc, &argv );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  if(mpirank==mpisize-1){
        //gdb_lock();
    }

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

      stringstream  ss;
      ss << "logTest" << mpirank;
      statusOFS.open( ss.str().c_str() );

      NumMat<MYSCALAR> sbuf(2,5);
      NumMat<MYSCALAR> rbuf(2,5);
      SetValue(sbuf,MYSCALAR(0.0));
      SetValue(rbuf,MYSCALAR(0.0));

      //TEST BROADCAST TREES
      {
        Int numTrees =2;
        std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > >bcastTrees(numTrees);

        std::vector<int> treeIdx;
        for(Int t = 0; t<bcastTrees.size(); t++){
          std::vector<Int> ranks;
          if(t==bcastTrees.size()-1){
            ranks.push_back(0);
            for(Int i = 1; i<mpisize; i++){
              ranks.push_back(i);
            }
          }
          else{
            ranks.push_back(t%2);
            for(Int i = t%2+2; i<mpisize; i+=2){
              ranks.push_back(i);
            }
          }
          statusOFS<<"Ranks: "<<ranks<<std::endl;
          Int msgSize = sbuf.Size();
          Int seed = mpisize-1;
          std::sort(ranks.begin()+1,ranks.end());

          if(std::find(ranks.begin(),ranks.end(),mpirank)!=ranks.end()){
            double rseed = (double)seed / (double)mpisize;
            auto & bcastTree = bcastTrees[t];
            bcastTree.reset(TreeBcast_v2<MYSCALAR>::Create(world_comm,ranks.data(),ranks.size(),msgSize,rseed));
            bcastTree->SetTag(t);
          }
          treeIdx.push_back(t);
        }

        for(Int t = 0; t<bcastTrees.size(); t++){
          auto & bcastTree =  bcastTrees[t];
          if(bcastTree!=nullptr){
            if(bcastTree->IsRoot()){
              SetValue(sbuf,MYSCALAR((double)t+1));
              //statusOFS<<"sbuf: "<<sbuf<<std::endl;
              bcastTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              bcastTree->SetLocalBuffer(rbuf.Data());
            }
          }
        }

        double elapsed = -MPI_Wtime();

        //Now test the tree
        for(Int t = 0; t<bcastTrees.size(); t++){
          auto & bcastTree =  bcastTrees[t];
          if(bcastTree!=nullptr){
            if(bcastTree->IsRoot()){
              bcastTree->SetDataReady(true);
            }
            bcastTree->Progress(); 
          }
        }

        //Use the Waitsome function
        std::list<int> rdyIdx;
        std::vector<bool> bcastdone(numTrees,false);

        bool all_done = std::all_of(bcastdone.begin(), bcastdone.end(), [](bool v) { return v; });
        while(!all_done){
          TreeBcast_Waitsome( treeIdx, bcastTrees, rdyIdx, bcastdone);
          for(auto idx : rdyIdx){
            auto & bcastTree = bcastTrees[idx];
            bcastTree->cleanupBuffers();
          }
          all_done = std::all_of(bcastdone.begin(), bcastdone.end(), [](bool v) { return v; });
        }

        elapsed += MPI_Wtime();

        statusOFS<<"Elapsed time: "<<elapsed<<std::endl;

        check_broadcasts(world_comm,bcastTrees, sbuf, rbuf);

        //Reset the trees
        for(Int t = 0; t<bcastTrees.size(); t++){
          auto & bcastTree =  bcastTrees[t];
          if(bcastTree!=nullptr){
            bcastTree->Reset();
          }
        }

        for(Int t = 0; t<bcastTrees.size(); t++){
          auto & bcastTree =  bcastTrees[t];
          if(bcastTree!=nullptr){
            if(bcastTree->IsRoot()){
              SetValue(sbuf,MYSCALAR((double)t+42));
              bcastTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              bcastTree->SetLocalBuffer(rbuf.Data());
            }
          }
        }

        //Disable overlap
        elapsed = -MPI_Wtime();

        for(Int t = 0; t<bcastTrees.size(); t++){
          auto & bcastTree =  bcastTrees[t];
          if(bcastTree!=nullptr){
            if(bcastTree->IsRoot()){
              bcastTree->SetDataReady(true);
            }
            bcastTree->Wait();
            bcastTree->cleanupBuffers();
          }
        }

        elapsed += MPI_Wtime();

        statusOFS<<"Elapsed time: "<<elapsed<<std::endl;

        check_broadcasts(world_comm,bcastTrees, sbuf, rbuf);
      }

      //TEST REDUCE TREES
      {
        Int numTrees =2;
        std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > >redTrees(numTrees);

        std::vector<int> treeIdx;
        for(Int t = 0; t<redTrees.size(); t++){
          std::vector<Int> ranks;

          if(t==redTrees.size()-1){
            ranks.push_back(0);
            for(Int i = 1; i<mpisize; i++){
              ranks.push_back(i);
            }
          }
          else{
            ranks.push_back(t%2);
            for(Int i = t%2+2; i<mpisize; i+=2){
              ranks.push_back(i);
            }
          }
          statusOFS<<"Ranks: "<<ranks<<std::endl;
          Int msgSize = sbuf.Size();
          Int seed = mpisize-1;
          std::sort(ranks.begin()+1,ranks.end());

          if(std::find(ranks.begin(),ranks.end(),mpirank)!=ranks.end()){
            double rseed = (double)seed / (double)mpisize;
            auto & redTree = redTrees[t];
            redTree.reset(TreeReduce_v2<MYSCALAR>::Create(world_comm,ranks.data(),ranks.size(),msgSize,rseed));
            redTree->SetTag(t);
          }
          treeIdx.push_back(t);
        }

        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            if(!redTree->IsRoot()){
              SetValue(sbuf,MYSCALAR((double)t+1));
              //statusOFS<<"sbuf: "<<sbuf<<std::endl;
              redTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              SetValue(rbuf,MYSCALAR((double)0.0));
//              SetValue(rbuf,MYSCALAR((double)t+1));
              redTree->SetLocalBuffer(rbuf.Data());
            }
          }
        }

        double elapsed = -MPI_Wtime();


        //Now test the tree
        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            if(!redTree->IsRoot()){
              redTree->SetDataReady(true);
            }
            redTree->Progress(); 
          }
        }

        //Use the Waitsome function
        std::list<int> rdyIdx;
        std::vector<bool> reddone(numTrees,false);

        bool all_done = std::all_of(reddone.begin(), reddone.end(), [](bool v) { return v; });
        while(!all_done){
          TreeBcast_Waitsome( treeIdx, redTrees, rdyIdx, reddone);
          for(auto idx : rdyIdx){
            auto & redTree = redTrees[idx];
            redTree->cleanupBuffers();
          }
          all_done = std::all_of(reddone.begin(), reddone.end(), [](bool v) { return v; });
        }

        elapsed += MPI_Wtime();

        statusOFS<<"Elapsed time: "<<elapsed<<std::endl;


////        //validate the output
////
////        NumMat<MYSCALAR> checkBuf(sbuf.m(),sbuf.n());
////        NumMat<MYSCALAR> tmpBuf(sbuf.m(),sbuf.n());
////        SetValue(checkBuf,MYSCALAR(0.0));
////
////        for(Int t = 0; t<redTrees.size(); t++){
////          Int root = t%2;
////
////          
////          std::vector<Int> ranks;
////          for(Int i = t%2+2; i<mpisize; i+=2){
////            ranks.push_back(i);
////          }
////
////          for(auto p : ranks){
////            if(mpirank==root){
////              MPI_Recv(tmpBuf.Data(),tmpBuf.ByteSize(),MPI_BYTE,MPI_ANY_SOURCE,t,world_comm,MPI_STATUS_IGNORE);
////              //axpy
////              blas::Axpy(tmpBuf.Size(), ONE<MYSCALAR>(), tmpBuf.Data(), 1, checkBuf.Data(), 1 );
////            }
////            else{
////              MPI_Send(sbuf.Data(),sbuf.ByteSize(),MPI_BYTE,root,t,world_comm,MPI_STATUS_IGNORE);
////            }
////          }
////
////          auto & redTree =  redTrees[t];
////          if(redTree!=nullptr){
////            if(redTree->IsRoot()){
////              for(Int i = 0; i<checkBuf.m(); i++){
////                for(Int j = 0; j<checkBuf.n(); j++){
////                  assert(checkBuf(i,j)==rbuf(i,j));
////                }
////              }
////            }
////          }
////        }
////        statusOFS<<"Reduce trees test passed"<<std::endl;


        check_reduces(world_comm,redTrees, sbuf, rbuf);



        //Reset the trees
        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            redTree->Reset();
          }
        }

        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            if(!redTree->IsRoot()){
              SetValue(sbuf,MYSCALAR((double)t+42));
              redTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              redTree->SetLocalBuffer(rbuf.Data());
            }
          }
        }

        //Disable overlap
        elapsed = -MPI_Wtime();

        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            if(!redTree->IsRoot()){
              redTree->SetDataReady(true);
            }
            redTree->Wait();
            redTree->cleanupBuffers();
          }
        }

        elapsed += MPI_Wtime();

        statusOFS<<"Elapsed time: "<<elapsed<<std::endl;
        statusOFS<<"rbuf: "<<rbuf<<std::endl;
      }

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

      statusOFS.close();
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

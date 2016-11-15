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

#define FTREE_LIMIT 1

#include  "ppexsi.hpp"
#include "pexsi/timer.h"

//#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#else
#define MYSCALAR Real
#endif

#include <random>

using namespace PEXSI;
using namespace std;



bool check_broadcasts(MPI_Comm world_comm,
      const std::vector<std::vector<Int > > & arrRanks, 
      const std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > > & bcastTrees, 
            const std::vector<NumMat<MYSCALAR> > & sbufs, const std::vector<NumMat<MYSCALAR> > &rbufs){
  //validate the output
  int mpisize = 0;
  int mpirank = 0;
      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );


  bool all_valid = true;

  NumMat<MYSCALAR> checkBuf;
  for(Int t = 0; t<bcastTrees.size(); t++){

    auto & ranks = arrRanks[t];
    auto & sbuf = sbufs[t];
    auto & rbuf = rbufs[t];

    Int root = ranks[0];


    if(mpirank==root){
      checkBuf.Resize(sbuf.m(),sbuf.n());
      SetValue(checkBuf,MYSCALAR(0.0));
    }
    else{
      checkBuf.Resize(rbuf.m(),rbuf.n());
      SetValue(checkBuf,MYSCALAR(0.0));
    }


    for(auto p : ranks){
      if(mpirank==root){
        if(mpirank!=p){
          MPI_Send(sbuf.Data(),sbuf.ByteSize(),MPI_BYTE,p,t,world_comm);
        }
      }
      else if(mpirank == p){
        MPI_Recv(checkBuf.Data(),checkBuf.ByteSize(),MPI_BYTE,root,t,world_comm,MPI_STATUS_IGNORE);
      }
    }

    MPI_Barrier(world_comm);

    bool valid = true;
    auto & bcastTree =  bcastTrees[t];
    if(bcastTree!=nullptr){
      if(!bcastTree->IsRoot()){
        for(Int i = 0; i<checkBuf.m(); i++){
          for(Int j = 0; j<checkBuf.n(); j++){
            valid = valid && checkBuf(i,j)==rbuf(i,j);
          }
        }
        if(mpirank!=0){
          MPI_Send(&valid,sizeof(valid),MPI_BYTE,0,t,world_comm);
        }
      }
    }
    
    if(mpirank==0){
      for(int ip = 1; ip<ranks.size();ip++){
        auto p = ranks[ip];
        if(p!=mpirank){
          bool tvalid = false;
          MPI_Recv(&tvalid,sizeof(valid),MPI_BYTE,p,t,world_comm,MPI_STATUS_IGNORE);
          valid = valid && tvalid;
        }
      }
    }

    if(!valid){
      if(mpirank==0){
        std::cout<<"Bcast Tree "<<t<<" failed"<<std::endl;
      }
      if(mpirank!=root){
        statusOFS<<rbuf<<std::endl;
        statusOFS<<checkBuf<<std::endl;
      }
    }
    all_valid = all_valid && valid;
  }

  MPI_Barrier(world_comm);
  return all_valid;
}


bool check_reduces(MPI_Comm world_comm,
      const std::vector<std::vector<Int > > & arrRanks, 
      const std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > > & redTrees, 
      const std::vector<NumMat<MYSCALAR> > & sbufs, const std::vector<NumMat<MYSCALAR> > &rbufs){
  //validate the output
  int mpisize = 0;
  int mpirank = 0;
  MPI_Comm_rank(world_comm, &mpirank );
  MPI_Comm_size(world_comm, &mpisize );

  NumMat<MYSCALAR> checkBuf;
  NumMat<MYSCALAR> tmpBuf;

  bool all_valid = true;
  for(Int t = 0; t<redTrees.size(); t++){
    auto & ranks = arrRanks[t];
    auto & sbuf = sbufs[t];
    auto & rbuf = rbufs[t];

    Int root = ranks[0];


    if(mpirank==root){
      checkBuf.Resize(rbuf.m(),rbuf.n());
      tmpBuf.Resize(rbuf.m(),rbuf.n());
      SetValue(checkBuf,MYSCALAR(0.0));
    }


    for(auto p : ranks){
      if(mpirank==root){
        if(mpirank!=p){
          MPI_Status stat;
          MPI_Recv(tmpBuf.Data(),tmpBuf.ByteSize(),MPI_BYTE,MPI_ANY_SOURCE,t,world_comm,&stat);
          //axpy
          statusOFS<<"Tree "<<t<<" Recv from P"<<stat.MPI_SOURCE<<": "<<tmpBuf<<std::endl;
          blas::Axpy(tmpBuf.Size(), ONE<MYSCALAR>(), tmpBuf.Data(), 1, checkBuf.Data(), 1 );
        }
//        else{
//          blas::Axpy(rbuf.Size(), ONE<MYSCALAR>(), rbuf.Data(), 1, checkBuf.Data(), 1 );
//        }
      }
      else if(mpirank == p){
        statusOFS<<"Tree "<<t<<" Send to P"<<root<<": "<<sbuf<<std::endl;
        MPI_Send(sbuf.Data(),sbuf.ByteSize(),MPI_BYTE,root,t,world_comm);
      }
    }

    MPI_Barrier(world_comm);

    bool valid = true;
    auto & redTree =  redTrees[t];

    if(redTree!=nullptr){
      if(redTree->IsRoot()){
        for(Int i = 0; i<checkBuf.m(); i++){
          for(Int j = 0; j<checkBuf.n(); j++){
            valid = valid && checkBuf(i,j)==rbuf(i,j);
          }
        }
        if(mpirank!=0){
          MPI_Send(&valid,sizeof(valid),MPI_BYTE,0,t,world_comm);
        }
      }
    }
    else if(mpirank==0){
       MPI_Recv(&valid,sizeof(valid),MPI_BYTE,root,t,world_comm,MPI_STATUS_IGNORE);
    }

    if(!valid){
      if(mpirank==0){
        std::cout<<"Reduce Tree "<<t<<" failed"<<". Root is P"<<root<< std::endl;
      }
      if(mpirank==root){
        statusOFS<<rbuf<<std::endl;
        statusOFS<<checkBuf<<std::endl;
      }
    }

    all_valid = all_valid && valid;
  }

  MPI_Barrier(world_comm);
  return all_valid;
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
    std::random_device rd;
    std::mt19937 gen(rd());
             
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
        Int numTrees =10;
        std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > >bcastTrees(numTrees);
        std::vector<std::vector<Int> > arrRanks(numTrees);
        std::vector<NumMat<MYSCALAR> > sbufs(numTrees);
        std::vector<NumMat<MYSCALAR> > rbufs(numTrees);


        std::vector<int> treeIdx;
        for(Int t = 0; t<bcastTrees.size(); t++){
          auto & ranks = arrRanks[t];
          if(t==bcastTrees.size()-1){
            ranks.push_back(0);
            for(Int i = 1; i<mpisize; i++){
              ranks.push_back(i);
            }
          }
          else{
            ranks.resize(mpisize);
            
            if(mpirank==0){
              std::iota(ranks.begin(),ranks.end(),0);
              std::random_shuffle ( ranks.begin(), ranks.end() );
              std::uniform_int_distribution<> dis(1, ranks.size());
              int count = dis(gen);
              MPI_Bcast(&count,sizeof(count),MPI_BYTE,0,world_comm);
              ranks.resize(count);
              MPI_Bcast(ranks.data(),count*sizeof(Int),MPI_BYTE,0,world_comm);
            }
            else{
              int count = 0;
              MPI_Bcast(&count,sizeof(count),MPI_BYTE,0,world_comm);
              ranks.resize(count);
              MPI_Bcast(ranks.data(),count*sizeof(Int),MPI_BYTE,0,world_comm);
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
              auto & sbuf = sbufs[t];
              sbuf.Resize(rbuf.m(),rbuf.n());
              SetValue(sbuf,MYSCALAR((double)t+1));
              bcastTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              auto & rbuf = rbufs[t];
              rbuf.Resize(sbuf.m(),sbuf.n());
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

        if(check_broadcasts(world_comm,arrRanks,bcastTrees, sbufs, rbufs)){
          if(mpirank==0){
            std::cout<<"Bcast trees test passed"<<std::endl;
          }
        }
        else{
          if(mpirank==0){
            std::cout<<"Bcast trees test failed"<<std::endl;
            return -1;
          }
        }




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
              auto & sbuf = sbufs[t];
              sbuf.Resize(rbuf.m(),rbuf.n());
              SetValue(sbuf,MYSCALAR((double)t+42));
              bcastTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              auto & rbuf = rbufs[t];
              rbuf.Resize(sbuf.m(),sbuf.n());
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

        if(check_broadcasts(world_comm,arrRanks,bcastTrees, sbufs, rbufs)){
          if(mpirank==0){
            std::cout<<"Bcast trees test passed"<<std::endl;
          }
        }
        else{
          if(mpirank==0){
            std::cout<<"Bcast trees test failed"<<std::endl;
            return -1;
          }
        }


      }

      //TEST REDUCE TREES
      {
        Int numTrees =10;
        std::vector<std::unique_ptr<TreeBcast_v2<MYSCALAR> > >redTrees(numTrees);
        std::vector<std::vector<Int> > arrRanks(numTrees);
        std::vector<NumMat<MYSCALAR> > rbufs(numTrees);
        std::vector<NumMat<MYSCALAR> > sbufs(numTrees);

        std::vector<int> treeIdx;
        for(Int t = 0; t<redTrees.size(); t++){
          auto & ranks = arrRanks[t];

          


          if(t==redTrees.size()-1){
            ranks.push_back(0);
            for(Int i = 1; i<mpisize; i++){
              ranks.push_back(i);
            }
          }
          else{

            ranks.resize(mpisize);
            
            if(mpirank==0){
              std::iota(ranks.begin(),ranks.end(),0);
              std::random_shuffle ( ranks.begin(), ranks.end() );
              std::uniform_int_distribution<> dis(1, ranks.size());
              int count = dis(gen);
              MPI_Bcast(&count,sizeof(count),MPI_BYTE,0,world_comm);
              ranks.resize(count);
              MPI_Bcast(ranks.data(),count*sizeof(Int),MPI_BYTE,0,world_comm);
            }
            else{
              int count = 0;
              MPI_Bcast(&count,sizeof(count),MPI_BYTE,0,world_comm);
              ranks.resize(count);
              MPI_Bcast(ranks.data(),count*sizeof(Int),MPI_BYTE,0,world_comm);
            }
            

          }
          statusOFS<<"Reduce "<<t<<" Ranks: "<<ranks<<std::endl;
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
              auto & sbuf = sbufs[t];
              sbuf.Resize(rbuf.m(),rbuf.n());
              SetValue(sbuf,MYSCALAR((double)t+1));
              redTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              auto & rbuf = rbufs[t];
              rbuf.Resize(sbuf.m(),sbuf.n());
              SetValue(rbuf,MYSCALAR((double)0.0));
              redTree->SetLocalBuffer(rbuf.Data());
            }
          }
        }

        //sbufs_back = sbufs;

        double elapsed = -MPI_Wtime();

        //Now test the tree
        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            redTree->SetDataReady(true);
            //if(!redTree->IsRoot()){
              redTree->SetDataReady(true);
            //}
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

        if(check_reduces(world_comm,arrRanks,redTrees, sbufs, rbufs)){
          if(mpirank==0){
            std::cout<<"Reduce trees test passed"<<std::endl;
          }
        }
        else{
          if(mpirank==0){
            std::cout<<"Reduce trees test failed"<<std::endl;
            return -1;
          }
        }

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
              auto & sbuf = sbufs[t];
              sbuf.Resize(rbuf.m(),rbuf.n());
              SetValue(sbuf,MYSCALAR((double)t+42));
              redTree->SetLocalBuffer(sbuf.Data());
            }
            else{
              auto & rbuf = rbufs[t];
              rbuf.Resize(sbuf.m(),sbuf.n());
              SetValue(rbuf,MYSCALAR((double)0.0));
              redTree->SetLocalBuffer(rbuf.Data());
            }
          }
        }


        //Disable overlap
        elapsed = -MPI_Wtime();

        for(Int t = 0; t<redTrees.size(); t++){
          auto & redTree =  redTrees[t];
          if(redTree!=nullptr){
            //if(!redTree->IsRoot()){
              redTree->SetDataReady(true);
            //}
            redTree->Wait();
            redTree->cleanupBuffers();
          }
        }

        elapsed += MPI_Wtime();

        statusOFS<<"Elapsed time: "<<elapsed<<std::endl;

        if(check_reduces(world_comm,arrRanks,redTrees, sbufs, rbufs)){
          if(mpirank==0){
            std::cout<<"Reduce trees test passed"<<std::endl;
          }
        }
        else{
          if(mpirank==0){
            std::cout<<"Reduce trees test failed"<<std::endl;
            return -1;
          }
        }
      }

      if( mpirank == 0 )
        cout << "All tests were successfull"<< endl;

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

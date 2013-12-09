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
/// @file pselinv.cpp
/// @brief Implementation of the parallel SelInv.
/// @date 2013-08-05
#include "pselinv.hpp"

#define IDX_TO_TAG(lidx,tag) (SELINV_TAG_COUNT*(lidx)+(tag)) 
#define TAG_TO_IDX(tag,typetag) (((tag)-(typetag))/SELINV_TAG_COUNT) 

#define MPI_MAX_COMM 6000

#define HYBRID
#define BCAST_THRESHOLD 20


#ifdef USE_TAU
#include "TAU.h"
#elif defined (PROFILE) || defined(PMPI)
#define TAU
#include "timer.h"
#endif

#define VAL(str) #str
#define TOSTRING(str) VAL(str)


#ifdef USE_TAU 
#define TIMER_START(a) TAU_START(TOSTRING(a));
#define TIMER_STOP(a) TAU_STOP(TOSTRING(a));
#elif defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#else
#define TIMER_START(a)
#define TIMER_STOP(a)
#endif



#define MOD(a,b) \
  ( ((a)%(b)+(b))%(b))






















namespace PEXSI{




template <typename F> inline std::vector<char>::iterator serialize( F val, const std::vector<char>::iterator & bufpos)
{
  std::copy((char *)&val,(char*)&val+sizeof(F),bufpos);
  return bufpos+sizeof(F);
}




template <typename F> inline std::vector<char>::iterator deserialize( F & val, const std::vector<char>::iterator & bufpos)
{
  val = ((F*)&(*bufpos))[0];
  return bufpos+sizeof(F);
}


template <typename F> inline std::vector<char>::iterator serialize( NumMat<F> & mat, const std::vector<char>::iterator & bufpos)
{
  std::vector<char>::iterator curpos = bufpos;
  //serialize the size first
  curpos = serialize( mat.m(), curpos);
  curpos = serialize( mat.n(), curpos);
  //serialize the content
  std::copy((char*)mat.Data(),(char*)mat.Data()+mat.ByteSize(),curpos);
  return curpos+mat.ByteSize();
}



template <typename F> inline std::vector<char>::iterator deserialize( NumMat<F> & mat, const std::vector<char>::iterator & bufpos)
{
  std::vector<char>::iterator curpos = bufpos;
  //deserialize the size first
  Int mat_n,mat_m;
  curpos = deserialize( mat_m, curpos);
  curpos = deserialize( mat_n, curpos);
  //allocate the array
  mat.Resize(mat_m,mat_n);
  std::copy(curpos,curpos+mat_n*mat_m*sizeof(F),mat.Data());
  return curpos+mat_n*mat_m*sizeof(F);
}

template <typename F> inline std::vector<char>::iterator serialize( Int mat_m, Int mat_n, F * mat, const std::vector<char>::iterator & bufpos)
{
  std::vector<char>::iterator curpos = bufpos;
  //serialize the size first
  curpos = serialize( mat_m, curpos);
  curpos = serialize( mat_n, curpos);
  if(mat!=NULL){
    //serialize the content
    std::copy((char*)mat,(char*)mat+mat_m*mat_n*sizeof(F),curpos);
  }
  return curpos+mat_n*mat_m*sizeof(F);
}

template <typename F> inline std::vector<char>::iterator deserialize( Int & mat_m, Int & mat_n, F * & mat, const std::vector<char>::iterator & bufpos)
{
  std::vector<char>::iterator curpos = bufpos;
  //deserialize the size first
  curpos = deserialize( mat_m, curpos);
  curpos = deserialize( mat_n, curpos);
  mat = (F*)&(*curpos);  
  return curpos+mat_n*mat_m*sizeof(F);
}

template <typename F> inline std::vector<char>::iterator serialize( Int vec_size, F * vec, const std::vector<char>::iterator & bufpos)
{
  std::vector<char>::iterator curpos = bufpos;
  //serialize the size first
  curpos = serialize( vec_size, curpos);
  if(vec!=NULL){
    //serialize the content
    std::copy((char*)vec,(char*)vec+vec_size*sizeof(F),curpos);
  }
  return curpos+vec_size*sizeof(F);
}

template <typename F> inline std::vector<char>::iterator deserialize( Int & vec_size, F * & vec, const std::vector<char>::iterator & bufpos)
{
  std::vector<char>::iterator curpos = bufpos;
  //deserialize the size first
  curpos = deserialize( vec_size, curpos);
  vec = (F*)&(*curpos);  
  return curpos+vec_size*sizeof(F);
}

















  GridType::GridType	( MPI_Comm Bcomm, int nprow, int npcol )
  {
#ifndef _RELEASE_
    PushCallStack("GridType::GridType");
#endif
    Int info;
    MPI_Initialized( &info );
    if( !info ){
      throw std::logic_error( "MPI has not been initialized." );
    }
    MPI_Group  comm_group;
    MPI_Comm_group( Bcomm, &comm_group );
    MPI_Comm_create( Bcomm, comm_group, &comm );
    //		comm = Bcomm;

    MPI_Comm_rank( comm, &mpirank );
    MPI_Comm_size( comm, &mpisize );
    if( mpisize != nprow * npcol ){
      throw std::logic_error( "mpisize != nprow * npcol." ); 
    }

    numProcRow = nprow;
    numProcCol = npcol;

    Int myrow = mpirank / npcol;
    Int mycol = mpirank % npcol;

    MPI_Comm_split( comm, myrow, mycol, &rowComm );
    MPI_Comm_split( comm, mycol, myrow, &colComm );

    MPI_Group_free( &comm_group );

#ifndef _RELEASE_
    PopCallStack();
#endif

    return ;
  } 		// -----  end of method GridType::GridType  ----- 


  GridType::~GridType	(  )
  {
#ifndef _RELEASE_
    PushCallStack("GridType::~GridType");
#endif
    // Dot not free grid.comm which is not generated by GridType().

    MPI_Comm_free( &rowComm );
    MPI_Comm_free( &colComm ); 
    MPI_Comm_free( &comm );

#ifndef _RELEASE_
    PopCallStack();
#endif
    return ;
  } 		// -----  end of method GridType::~GridType  ----- 

} // namespace PEXSI


namespace PEXSI{


  PMatrix::PMatrix ( const GridType* g, const SuperNodeType* s, const PEXSI::SuperLUOptions * o ):grid_(g), super_(s), options_(o)
  {
#ifndef _RELEASE_
    PushCallStack("PMatrix::PMatrix");
#endif

    //    if( grid_->numProcRow != grid_->numProcCol ){
    //      throw std::runtime_error( "The current version of SelInv only works for square processor grids." ); }


    L_.resize( this->NumLocalBlockCol() );
    U_.resize( this->NumLocalBlockRow() );
    //workingSet_.resize(this->NumSuper());
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "PMatrix is constructed. The grid information: " << std::endl;
    statusOFS << "mpirank = " << MYPROC(grid_) << std::endl;
    statusOFS << "myrow   = " << MYROW(grid_) << std::endl; 
    statusOFS << "mycol   = " << MYCOL(grid_) << std::endl; 
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif
    return ;
  } 		// -----  end of method PMatrix::PMatrix  ----- 

  PMatrix::~PMatrix() {}	



  void PMatrix::GetEtree(std::vector<Int> & etree_supno )
  {

#ifndef _RELEASE_
    PushCallStack("PMatrix::GetEtree");
    double begin =  MPI_Wtime( );
#endif
    Int nsupers = this->NumSuper();

    if( options_->ColPerm != "PARMETIS" ) {
      /* Use the etree computed from serial symb. fact., and turn it
         into supernodal tree.  */
      const SuperNodeType * superNode = this->SuperNode();


      //translate from columns to supernodes etree using supIdx
      etree_supno.resize(this->NumSuper());
      for(Int i = 0; i < superNode->etree.m(); ++i){
        Int curSnode = superNode->superIdx[i];
        Int parentSnode = (superNode->etree[i]>= superNode->etree.m()) ?this->NumSuper():superNode->superIdx[superNode->etree[i]];
        if( curSnode != parentSnode){
          etree_supno[curSnode] = parentSnode;
        }
      }

    } else { /* ParSymbFACT==YES and SymPattern==YES  and RowPerm == NOROWPERM */
      /* Compute an "etree" based on struct(L), 
         assuming struct(U) = struct(L').   */

      /* find the first block in each supernodal-column of local L-factor */
      std::vector<Int> etree_supno_l( nsupers, nsupers  );
      for( Int ksup = 0; ksup < nsupers; ksup++ ){
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          // L part
          std::vector<LBlock>& Lcol = this->L( LBj(ksup, grid_) );
          if(Lcol.size()>0){
            Int firstBlk = 0;
            if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
              firstBlk=1;

            for( Int ib = firstBlk; ib < Lcol.size(); ib++ ){
              etree_supno_l[ksup] = std::min(etree_supno_l[ksup] , Lcol[ib].blockIdx);
            }
          }
        }
      }


#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << " Local supernodal elimination tree is " << etree_supno_l <<std::endl<<std::endl;

#endif
      /* form global e-tree */
      etree_supno.resize( nsupers );
      mpi::Allreduce( (Int*) &etree_supno_l[0],(Int *) &etree_supno[0], nsupers, MPI_MIN, grid_->comm );
      etree_supno[nsupers-1]=nsupers;
    }

#ifndef _RELEASE_
    double end =  MPI_Wtime( );
    statusOFS<<"Building the list took "<<end-begin<<"s"<<std::endl;
#endif
#ifndef _RELEASE_
    PopCallStack();
#endif
  } 		// -----  end of method PMatrix::GetEtree  ----- 


  void PMatrix::ConstructCommunicationPattern_P2p	(  )
  {
#ifndef _RELEASE_
    PushCallStack("PMatrix::ConstructCommunicationPattern_P2p");
#endif
    Int numSuper = this->NumSuper();
#ifndef _RELEASE_
    PushCallStack( "Initialize the communication pattern" );
#endif
    isSendToBelow_.Resize(grid_->numProcRow, numSuper);
    isSendToRight_.Resize(grid_->numProcCol, numSuper);
    isSendToDiagonal_.Resize( numSuper );
    SetValue( isSendToBelow_, false );
    SetValue( isSendToRight_, false );
    SetValue( isSendToDiagonal_, false );

    isSendToCrossDiagonal_.Resize(grid_->numProcCol+1, numSuper );
    SetValue( isSendToCrossDiagonal_, false );
    isRecvFromCrossDiagonal_.Resize(grid_->numProcRow+1, numSuper );
    SetValue( isRecvFromCrossDiagonal_, false );

    isRecvFromAbove_.Resize( numSuper );
    isRecvFromLeft_.Resize( numSuper );
    isRecvFromBelow_.Resize( grid_->numProcRow, numSuper );
    SetValue( isRecvFromAbove_, false );
    SetValue( isRecvFromBelow_, false );
    SetValue( isRecvFromLeft_, false );
#ifndef _RELEASE_
    PopCallStack();
#endif


      std::vector<Int> snodeEtree(this->NumSuper());
      GetEtree(snodeEtree);

#ifndef _RELEASE_
    PushCallStack( "Local column communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local column communication" << std::endl;
#endif
    // localColBlockRowIdx stores the nonzero block indices for each local block column.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockCol x numNonzeroBlock
    std::vector<std::set<Int> >   localColBlockRowIdx;

    localColBlockRowIdx.resize( this->NumLocalBlockCol() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        std::vector<Int>  tBlockRowIdx;
        tBlockRowIdx.clear();

        // L part
        std::vector<LBlock>& Lcol = this->L( LBj(ksup, grid_) );
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          tBlockRowIdx.push_back( Lcol[ib].blockIdx );
        }

        // U part
        for( Int ib = 0; ib < this->NumLocalBlockRow(); ib++ ){
          std::vector<UBlock>& Urow = this->U(ib);
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            if( Urow[jb].blockIdx == ksup ){
              tBlockRowIdx.push_back( GBi( ib, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockRowIdx;
        mpi::Allgatherv( tBlockRowIdx, tAllBlockRowIdx, grid_->colComm );

        localColBlockRowIdx[LBj( ksup, grid_ )].insert(
            tAllBlockRowIdx.begin(), tAllBlockRowIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Column block " << ksup 
          << " has the following nonzero block rows" << std::endl;
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for(ksup)


#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack( "Local row communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local row communication" << std::endl;
#endif
    // localRowBlockColIdx stores the nonzero block indices for each local block row.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockRow x numNonzeroBlock
    std::vector<std::set<Int> >   localRowBlockColIdx;

    localRowBlockColIdx.resize( this->NumLocalBlockRow() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        std::vector<Int>  tBlockColIdx;
        tBlockColIdx.clear();

        // U part
        std::vector<UBlock>& Urow = this->U( LBi(ksup, grid_) );
        for( Int jb = 0; jb < Urow.size(); jb++ ){
          tBlockColIdx.push_back( Urow[jb].blockIdx );
        }

        // L part
        for( Int jb = 0; jb < this->NumLocalBlockCol(); jb++ ){
          std::vector<LBlock>& Lcol = this->L(jb);
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            if( Lcol[ib].blockIdx == ksup ){
              tBlockColIdx.push_back( GBj( jb, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockColIdx;
        mpi::Allgatherv( tBlockColIdx, tAllBlockColIdx, grid_->rowComm );

        localRowBlockColIdx[LBi( ksup, grid_ )].insert(
            tAllBlockColIdx.begin(), tAllBlockColIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Row block " << ksup 
          << " has the following nonzero block columns" << std::endl;
        for( std::set<Int>::iterator si = localRowBlockColIdx[LBi( ksup, grid_ )].begin();
            si != localRowBlockColIdx[LBi( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
    } // for(ksup)

#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack("SendToBelow / RecvFromAbove");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes to the right of ksup


      Int jsup = snodeEtree[ksup];
      while(jsup<numSuper){
        Int jsupLocalBlockCol = LBj( jsup, grid_ );
        Int jsupProcCol = PCOL( jsup, grid_ );
        if( MYCOL( grid_ ) == jsupProcCol ){

          // SendToBelow / RecvFromAbove only if (ksup, jsup) is nonzero.
          if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ) {
            for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
                si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
              Int isup = *si;
              Int isupProcRow = PROW( isup, grid_ );
              if( isup > ksup ){
                if( MYROW( grid_ ) == isupProcRow ){
                  isRecvFromAbove_(ksup) = true;
                }
                if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
                  isSendToBelow_( isupProcRow, ksup ) = true;
                }
              } // if( isup > ksup )
            } // for (si)
          } // if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 )

        } // if( MYCOL( grid_ ) == PCOL( jsup, grid_ ) )
        jsup = snodeEtree[jsup];

      } // for(jsup)
    } // for(ksup)

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToBelow:" << std::endl;
    for(int j = 0;j< isSendToBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToBelow_.m();i++){
        statusOFS<< isSendToBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromAbove:" << std::endl;
    for(int j = 0;j< isRecvFromAbove_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromAbove_(j)<<std::endl;
    }
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif












#ifndef _RELEASE_
    PushCallStack("SendToRight / RecvFromLeft");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes below ksup

      Int isup = snodeEtree[ksup];
      while(isup<numSuper){
        Int isupLocalBlockRow = LBi( isup, grid_ );
        Int isupProcRow       = PROW( isup, grid_ );
        if( MYROW( grid_ ) == isupProcRow ){
          // SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.
          if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
            for( std::set<Int>::iterator si = localRowBlockColIdx[isupLocalBlockRow].begin();
                si != localRowBlockColIdx[isupLocalBlockRow].end(); si++ ){
              Int jsup = *si;
              Int jsupProcCol = PCOL( jsup, grid_ );
              if( jsup > ksup ){

                if( MYCOL( grid_ ) == jsupProcCol ){
                  isRecvFromLeft_(ksup) = true;
                }
                if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
                  isSendToRight_( jsupProcCol, ksup ) = true;
                }
              }
            } // for (si)
          } // if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 )
        } // if( MYROW( grid_ ) == isupProcRow )


        if( MYCOL( grid_ ) == PCOL(ksup, grid_) ){

          if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){ 
            isRecvFromBelow_(isupProcRow,ksup) = true;
          }    
          else if (MYROW(grid_) == isupProcRow){
            isSendToDiagonal_(ksup)=true;
          }    
        } // if( MYCOL( grid_ ) == PCOL(ksup, grid_) )
        isup = snodeEtree[isup];

      } // for (isup)
    }	 // for (ksup)


#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToRight:" << std::endl;
    for(int j = 0;j< isSendToRight_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToRight_.m();i++){
        statusOFS<< isSendToRight_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromLeft:" << std::endl;
    for(int j = 0;j< isRecvFromLeft_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromLeft_(j)<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromBelow:" << std::endl;
    for(int j = 0;j< isRecvFromBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isRecvFromBelow_.m();i++){
        statusOFS<< isRecvFromBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif








#ifndef _RELEASE_
    PushCallStack("SendToCrossDiagonal / RecvFromCrossDiagonal");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end(); si++ ){
          Int isup = *si;
          Int isupProcRow = PROW( isup, grid_ );
          Int isupProcCol = PCOL( isup, grid_ );
          if( isup > ksup && MYROW( grid_ ) == isupProcRow ){
            isSendToCrossDiagonal_(grid_->numProcCol, ksup ) = true;
            isSendToCrossDiagonal_(isupProcCol, ksup ) = true;
          }
        } // for (si)
      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for (ksup)

    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localRowBlockColIdx[ LBi(ksup, grid_) ].begin();
            si != localRowBlockColIdx[ LBi(ksup, grid_) ].end(); si++ ){
          Int jsup = *si;
          Int jsupProcCol = PCOL( jsup, grid_ );
          Int jsupProcRow = PROW( jsup, grid_ );
          if( jsup > ksup && MYCOL(grid_) == jsupProcCol ){
            isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) = true;
            isRecvFromCrossDiagonal_(jsupProcRow, ksup ) = true;
          }
        } // for (si)
      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
    } // for (ksup)
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToCrossDiagonal:" << std::endl;
    for(int j =0; j < isSendToCrossDiagonal_.n();j++){
      if(isSendToCrossDiagonal_(grid_->numProcCol,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isSendToCrossDiagonal_.m()-1;i++){
          if(isSendToCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(PROW(j,grid_),i,grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }

    statusOFS << std::endl << "isRecvFromCrossDiagonal:" << std::endl;
    for(int j =0; j < isRecvFromCrossDiagonal_.n();j++){
      if(isRecvFromCrossDiagonal_(grid_->numProcRow,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isRecvFromCrossDiagonal_.m()-1;i++){
          if(isRecvFromCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(i,PCOL(j,grid_),grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }


#endif

#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PopCallStack();
#endif

    //Build the list of supernodes based on the elimination tree from SuperLU
    GetWorkSet(snodeEtree,this->WorkingSet());
 
    return ;
  } 		// -----  end of method PMatrix::ConstructCommunicationPattern_P2p  ----- 

  void PMatrix::SelInv_P2p	(  )
  {
//#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
//    TAU_PROFILE_SET_CONTEXT(grid_->comm);
//#endif

    TIMER_START(SelInv_P2p);

#ifndef _RELEASE_
    PushCallStack("PMatrix::SelInv_P2p");
#endif


    Int numSuper = this->NumSuper(); 

    // Main loop
    std::vector<std::vector<Int> > & superList = this->WorkingSet();
    Int numSteps = superList.size();

    for (Int lidx=0; lidx<numSteps ; lidx++){
      Int stepSuper = superList[lidx].size(); 

      SelInvIntra_P2p(lidx);
    }

#ifndef _RELEASE_
    PopCallStack();
#endif

    TIMER_STOP(SelInv_P2p);

    return ;
  } 		// -----  end of method PMatrix::SelInv_P2p  ----- 




  void PMatrix::DestructCommunicators_Collectives	(  )
  {
#ifndef _RELEASE_
    PushCallStack("Destructing communicators");
#endif
    {
      for(int i = 0; i< commSendToBelow_.size(); ++i){
        if(commSendToBelow_[i]!=MPI_COMM_NULL){
          MPI_Comm_free(&commSendToBelow_[i]);
        }
      }

      for(int i = 0; i< commRecvFromBelow_.size(); ++i){
        if(commRecvFromBelow_[i]!=MPI_COMM_NULL){
          MPI_Comm_free(&commRecvFromBelow_[i]);
        }
      }

      for(int i = 0; i< commSendToRight_.size(); ++i){
        if(commSendToRight_[i]!=MPI_COMM_NULL){
          MPI_Comm_free(&commSendToRight_[i]);
        }
      }
    }
#ifndef _RELEASE_
    PopCallStack();
#endif


  } 		// -----  end of method PMatrix::DestructCommunicators_Collectives  ----- 



  inline void PMatrix::GetWorkSet(std::vector<Int> & snodeEtree, std::vector<std::vector<Int> > & WSet){
    TIMER_START(Compute_WorkSet);
    Int numSuper = this->NumSuper();


    if (options_->maxPipelineDepth!=1){





      //find roots in the supernode etree (it must be postordered)
      //initialize the parent we are looking at 
      Int rootParent = snodeEtree[numSuper-2];

      //compute the level of each supernode and the total number of levels
      IntNumVec level(numSuper);
      level(rootParent)=0;
      Int numLevel = 0; 
      for(Int i=rootParent-1; i>=0; i-- ){ level(i) = level(snodeEtree[i])+1; numLevel = std::max(numLevel, level(i)); }
      numLevel++;

      //Compute the number of supernodes at each level
      IntNumVec levelSize(numLevel);
      SetValue(levelSize,I_ZERO);
      for(Int i=rootParent-1; i>=0; i-- ){ levelSize(level(i))++; }

      //Allocate memory
      WSet.resize(numLevel,std::vector<Int>());
      for(Int i=0; i<numLevel; i++ ){WSet[i].reserve(levelSize(i));}

      //Fill the worklist based on the level of each supernode
      for(Int i=rootParent-1; i>=0; i-- ){
        WSet[level(i)].push_back(i);  
      }

      //Constrain the size of each list to be min(MPI_MAX_COMM,options_->maxPipelineDepth)
      Int limit = (options_->maxPipelineDepth>0)?std::min(MPI_MAX_COMM,options_->maxPipelineDepth):MPI_MAX_COMM;
      for (Int lidx=0; lidx<WSet.size() ; lidx++){
        if(WSet[lidx].size()>limit)
        {
          std::vector<std::vector<Int> >::iterator pos = WSet.begin()+lidx+1;               
          WSet.insert(pos,std::vector<Int>());
          WSet[lidx+1].insert(WSet[lidx+1].begin(),WSet[lidx].begin() + limit ,WSet[lidx].end());
          WSet[lidx].erase(WSet[lidx].begin()+limit,WSet[lidx].end());
        }
      }



    }
    else{
      for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
        WSet.push_back(std::vector<Int>());
        WSet.back().push_back(ksup);
      }

    }




    TIMER_STOP(Compute_WorkSet);
#if ( _DEBUGlevel_ >= 1 )
    for (Int lidx=0; lidx<WSet.size() ; lidx++){
      statusOFS << std::endl << "L"<< lidx << " is: {";
      for (Int supidx=0; supidx<WSet[lidx].size() ; supidx++){
        statusOFS << WSet[lidx][supidx] << " ["<<snodeEtree[WSet[lidx][supidx]]<<"] ";
      }
      statusOFS << " }"<< std::endl;
    }
#endif
  }










#ifdef WIP
  void PMatrix::ConstructCommunicationPattern_Collectives	(  )
  {
#ifndef _RELEASE_
    PushCallStack("PMatrix::ConstructCommunicationPattern_Collectives");
#endif
    Int numSuper = this->NumSuper();
#ifndef _RELEASE_
    PushCallStack( "Initialize the communication pattern" );
#endif
    isSendToBelow_.Resize(grid_->numProcRow, numSuper);
    isSendToRight_.Resize(grid_->numProcCol, numSuper);

    BolNumMat isSendToBelow2_(numSuper,grid_->numProcRow);
    BolNumMat isSendToRight2_(numSuper,grid_->numProcCol);
    BolNumMat isRecvFromBelow2_(numSuper,grid_->numProcRow);
    SetValue( isSendToBelow2_, false );
    SetValue( isSendToRight2_, false );
    SetValue( isRecvFromBelow2_, false );

    SetValue( isSendToBelow_, false );
    SetValue( isSendToRight_, false );


    NumVec<bool> isSendToAbove_( numSuper );
    SetValue( isSendToAbove_, false );

    isRecvFromAbove_.Resize( numSuper );
    isRecvFromLeft_.Resize( numSuper );
    SetValue( isRecvFromAbove_, false );
    SetValue( isRecvFromLeft_, false );

    isSendToCrossDiagonal_.Resize(grid_->numProcCol+1, numSuper );
    SetValue( isSendToCrossDiagonal_, false );
    isRecvFromCrossDiagonal_.Resize(grid_->numProcRow+1, numSuper );
    SetValue( isRecvFromCrossDiagonal_, false );

    countSendToBelow_.Resize(numSuper);
    countSendToRight_.Resize(numSuper);
    countRecvFromBelow_.Resize( numSuper );
    SetValue( countSendToBelow_, 0 );
    SetValue( countSendToRight_, 0 );
    SetValue( countRecvFromBelow_, 0 );

#ifndef _RELEASE_
    PopCallStack();
#endif



    std::vector<Int> snodeEtree(this->NumSuper());
    GetEtree(snodeEtree);

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << " Supernodal elimination tree is " << snodeEtree <<std::endl<<std::endl;
#endif




#ifndef _RELEASE_
    PushCallStack( "Local column communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local column communication" << std::endl;
#endif
    // localColBlockRowIdx stores the nonzero block indices for each local block column.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockCol x numNonzeroBlock
    std::vector<std::set<Int> >   localColBlockRowIdx;
    localColBlockRowIdx.resize( this->NumLocalBlockCol() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        std::vector<Int>  tBlockRowIdx;
        tBlockRowIdx.clear();

        // L part
        std::vector<LBlock>& Lcol = this->L( LBj(ksup, grid_) );
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          tBlockRowIdx.push_back( Lcol[ib].blockIdx );
        }

        // U part
        for( Int ib = 0; ib < this->NumLocalBlockRow(); ib++ ){
          std::vector<UBlock>& Urow = this->U(ib);
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            if( Urow[jb].blockIdx == ksup ){
              tBlockRowIdx.push_back( GBi( ib, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockRowIdx;
        mpi::Allgatherv( tBlockRowIdx, tAllBlockRowIdx, grid_->colComm );

        localColBlockRowIdx[LBj( ksup, grid_ )].insert(
            tAllBlockRowIdx.begin(), tAllBlockRowIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Column block " << ksup 
          << " has the following nonzero block rows" << std::endl;
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for(ksup)


#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack( "Local row communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local row communication" << std::endl;
#endif
    // localRowBlockColIdx stores the nonzero block indices for each local block row.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockRow x numNonzeroBlock
    std::vector<std::set<Int> >   localRowBlockColIdx;

    localRowBlockColIdx.resize( this->NumLocalBlockRow() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        std::vector<Int>  tBlockColIdx;
        tBlockColIdx.clear();

        // U part
        std::vector<UBlock>& Urow = this->U( LBi(ksup, grid_) );
        for( Int jb = 0; jb < Urow.size(); jb++ ){
          tBlockColIdx.push_back( Urow[jb].blockIdx );
        }

        // L part
        for( Int jb = 0; jb < this->NumLocalBlockCol(); jb++ ){
          std::vector<LBlock>& Lcol = this->L(jb);
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            if( Lcol[ib].blockIdx == ksup ){
              tBlockColIdx.push_back( GBj( jb, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockColIdx;
        mpi::Allgatherv( tBlockColIdx, tAllBlockColIdx, grid_->rowComm );

        localRowBlockColIdx[LBi( ksup, grid_ )].insert(
            tAllBlockColIdx.begin(), tAllBlockColIdx.end() );
#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Row block " << ksup 
          << " has the following nonzero block columns" << std::endl;
        for( std::set<Int>::iterator si = localRowBlockColIdx[LBi( ksup, grid_ )].begin();
            si != localRowBlockColIdx[LBi( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
    } // for(ksup)

#ifndef _RELEASE_
    PopCallStack();
#endif

#ifndef _RELEASE_
    PushCallStack("SendToBelow / RecvFromAbove");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes to the right of ksup


      Int jsup = snodeEtree[ksup];
      while(jsup<numSuper){
        Int jsupLocalBlockCol = LBj( jsup, grid_ );
        Int jsupProcCol = PCOL( jsup, grid_ );
        if( MYCOL( grid_ ) == jsupProcCol ){

          // SendToBelow / RecvFromAbove only if (ksup, jsup) is nonzero.
          if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ) {
            for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
                si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
              Int isup = *si;
              Int isupProcRow = PROW( isup, grid_ );
              if( isup > ksup ){
                if( MYROW( grid_ ) == isupProcRow ){
                  isRecvFromAbove_(ksup) = true;
                }
              } // if( isup > ksup )
            } // for (si)
          } // if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 )

        } // if( MYCOL( grid_ ) == PCOL( jsup, grid_ ) )
jsup = snodeEtree[jsup];

      } // for(jsup)
    } // for(ksup)




#ifndef _RELEASE_
    PushCallStack("SendToRight / RecvFromLeft");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes below ksup

      Int isup = snodeEtree[ksup];
      while(isup<numSuper){
        Int isupLocalBlockRow = LBi( isup, grid_ );
        Int isupProcRow       = PROW( isup, grid_ );
        if( MYROW( grid_ ) == isupProcRow ){
          // SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.

          if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
            for( std::set<Int>::iterator si = localRowBlockColIdx[isupLocalBlockRow].begin();
                si != localRowBlockColIdx[isupLocalBlockRow].end(); si++ ){

              Int jsup = *si;
              Int jsupProcCol = PCOL( jsup, grid_ );
              if( jsup > ksup ){

                if( MYCOL( grid_ ) == jsupProcCol ){
                  isRecvFromLeft_(ksup) = true;
                }
              }
            } // for (si)
          } // if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 )

        } // if( MYROW( grid_ ) == isupProcRow )

        isup = snodeEtree[isup];


      } // for (isup)
    }	 // for (ksup)



    //do a logical and between isRecvFromAbove and isRecvFromLeft to 
    //determine whether a processor takes part in the computations or not
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      isRecvFromLeft_(ksup) = (isRecvFromLeft_(ksup) &&  isRecvFromAbove_(ksup));// || ( MYCOL( grid_ ) == PCOL( ksup, grid_ ) );
      isRecvFromAbove_(ksup) = (isRecvFromAbove_(ksup) && isRecvFromLeft_(ksup));// || ( MYROW( grid_ ) == PROW( ksup, grid_ ) );
    }

    //reduce isRecvFromLeft and isRecvFromAbove to form isSendToRight and isSendToBelow
    MPI_Allgather(isRecvFromLeft_.Data(), numSuper*sizeof(bool),MPI_BYTE,isSendToRight2_.Data(), numSuper*sizeof(bool), MPI_BYTE, grid_->rowComm );
    MPI_Allgather(isRecvFromAbove_.Data(), numSuper*sizeof(bool),MPI_BYTE,isSendToBelow2_.Data(), numSuper*sizeof(bool), MPI_BYTE, grid_->colComm );

    Transpose(isSendToRight2_,isSendToRight_);
    Transpose(isSendToBelow2_,isSendToBelow_);

    for(int ksup=0; ksup<numSuper-1; ksup++){
      assert(isRecvFromLeft_(ksup)==isSendToRight2_(ksup,MYCOL(grid_)));
      assert(isRecvFromAbove_(ksup)==isSendToBelow2_(ksup,MYROW(grid_)));
    }

    for(int ksup=0; ksup<numSuper-1; ksup++){


#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<"["<<ksup<<"]      ";
      for(int i=0;i<grid_->numProcCol;i++){if(i==PCOL(ksup,grid_)){ statusOFS<<"["<<PNUM(MYROW(grid_),i,grid_)<<"] ";}else if(i==MYCOL(grid_)){ statusOFS<<"("<<PNUM(MYROW(grid_),i,grid_)<<") ";}else{statusOFS<<PNUM(MYROW(grid_),i,grid_)<<" ";}}
      statusOFS<<std::endl; 
      statusOFS<<"["<<ksup<<"] isSendToRight: ";
      for(int i=0;i<isSendToRight_.m();i++){statusOFS<<isSendToRight_(i,ksup)<<" ";}
      statusOFS<<std::endl;
      statusOFS<<"["<<ksup<<"] isRecvFromLeft: "<< isRecvFromLeft_(ksup)<<std::endl; 
      statusOFS<<"["<<ksup<<"]      ";


      for(int i=0;i<grid_->numProcRow;i++){if(i==PROW(ksup,grid_)){ statusOFS<<"["<<PNUM(i,MYCOL(grid_),grid_)<<"] ";}else if(i==MYROW(grid_)){ statusOFS<<"("<<PNUM(i,MYCOL(grid_),grid_)<<") ";}else{statusOFS<<PNUM(i,MYCOL(grid_),grid_)<<" ";}}
      statusOFS<<std::endl; 
      statusOFS<<"["<<ksup<<"] isSendToBelow: ";
      for(int i=0;i<isSendToBelow_.m();i++){statusOFS<<isSendToBelow_(i,ksup)<<" ";}
      statusOFS<<std::endl; 
      statusOFS<<"["<<ksup<<"] isRecvFromAbove: "<< isRecvFromAbove_(ksup)<<std::endl; 
#endif

      std::vector<bool> sTB(isSendToBelow_.VecData(ksup),isSendToBelow_.VecData(ksup+1));
      Int count= std::count(sTB.begin(), sTB.end(), true);
      Int color = 0;
      if(count>=1){
        if (!sTB[PROW(ksup,grid_)]){
          sTB[PROW(ksup,grid_)]=true;
          count++;
        }
        color = sTB[MYROW(grid_)];
        std::vector<Int> & snodeList = maskSendToBelow_[sTB];
        snodeList.push_back(ksup);
      }
      countSendToBelow_(ksup) = count * color;


      std::vector<bool> sTR(isSendToRight_.VecData(ksup),isSendToRight_.VecData(ksup+1));
      count= std::count(sTR.begin(), sTR.end(), true);

      color = 0;
      if(count>=1){
        if (!sTR[PCOL(ksup,grid_)]){
          sTR[PCOL(ksup,grid_)]=true;
          count++;
        }
        color = sTR[MYCOL(grid_)];
        std::vector<Int> & snodeList = maskSendToRight_[sTR];
        snodeList.push_back(ksup);
      }
      countSendToRight_(ksup) = count * color;

      isSendToAbove_(ksup)= (MYCOL(grid_)==PCOL(ksup,grid_)) && (countSendToRight_(ksup)>=1);


    }


    MPI_Allgather(isSendToAbove_.Data(), numSuper*sizeof(bool),MPI_BYTE,isRecvFromBelow2_.Data(), numSuper*sizeof(bool), MPI_BYTE, grid_->colComm );
    Transpose(isRecvFromBelow2_,isRecvFromBelow_);

    for(int ksup=0; ksup<numSuper-1; ksup++){


#if ( _DEBUGlevel_ >= 1 )
      for(int i=0;i<grid_->numProcRow;i++){if(i==PROW(ksup,grid_)){ statusOFS<<"["<<PNUM(i,MYCOL(grid_),grid_)<<"] ";}else if(i==MYROW(grid_)){ statusOFS<<"("<<PNUM(i,MYCOL(grid_),grid_)<<") ";}else{statusOFS<<PNUM(i,MYCOL(grid_),grid_)<<" ";}}
      statusOFS<<std::endl; 
      statusOFS<<"["<<ksup<<"] isRecvFromBelow: ";
      for(int i=0;i<isRecvFromBelow_.m();i++){statusOFS<<isRecvFromBelow_(i,ksup)<<" ";}
      statusOFS<<std::endl; 
#endif

      std::vector<bool> rFB(isRecvFromBelow_.VecData(ksup),isRecvFromBelow_.VecData(ksup+1));
      Int count= std::count(rFB.begin(), rFB.end(), true);
      Int color = 0;
      if(count>=1){
        if (!rFB[PROW(ksup,grid_)]){
          rFB[PROW(ksup,grid_)]=true;
          count++;
        }
        color = rFB[MYROW(grid_)];
        std::vector<Int> & snodeList = maskRecvFromBelow_[rFB];
        snodeList.push_back(ksup);
      }
      countRecvFromBelow_(ksup) = count * color;
    }

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToBelow:" << std::endl;
    for(int j = 0;j< isSendToBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToBelow_.m();i++){
        statusOFS<< isSendToBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromAbove:" << std::endl;
    for(int j = 0;j< isRecvFromAbove_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromAbove_(j)<<std::endl;
    }
#endif



#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToRight:" << std::endl;
    for(int j = 0;j< isSendToRight_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToRight_.m();i++){
        statusOFS<< isSendToRight_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromLeft:" << std::endl;
    for(int j = 0;j< isRecvFromLeft_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromLeft_(j)<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromBelow:" << std::endl;
    for(int j = 0;j< isRecvFromBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isRecvFromBelow_.m();i++){
        statusOFS<< isRecvFromBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }
#endif



#ifdef PRINT_COMMUNICATOR_STAT
    {
      statusOFS << std::endl << "countSendToBelow:" << countSendToBelow_ << std::endl;
      statusOFS << std::endl << "maskSendToBelow_:" << maskSendToBelow_.size() <<std::endl; 
      bitMaskSet::iterator it;
      for(it = maskSendToBelow_.begin(); it != maskSendToBelow_.end(); it++){
        //print the involved processors
        for(int curi = 0; curi < it->first.size(); curi++){
          statusOFS << it->first[curi] << " "; 
        }

        statusOFS<< "    ( ";
        //print the involved supernode indexes
        for(int curi = 0; curi < it->second.size(); curi++){
          statusOFS<< it->second[curi]<<" ";
        }

        statusOFS << ")"<< std::endl;
      }
    }
#endif
#ifndef _RELEASE_
    PopCallStack();
#endif

#ifndef _RELEASE_
    PushCallStack("Creating SendToBelow communicator");
#endif
    {
      MPI_Group colCommGroup;  
      MPI_Comm_group(grid_->colComm,&colCommGroup);

      commSendToBelow_.resize(maskSendToBelow_.size());
      commSendToBelowRoot_.resize(numSuper);
      commSendToBelowPtr_.resize(numSuper);

      Int commIdx = 0;
      bitMaskSet::iterator it;
      for(it = maskSendToBelow_.begin(); it != maskSendToBelow_.end(); it++){
        Int count= std::count(it->first.begin(), it->first.end(), true);

        if(count>1){
          MPI_Group commGroup;  

          Int color = it->first[MYROW(grid_)];

          MPI_Comm_split(grid_->colComm, color  ,MYROW(grid_) , &commSendToBelow_[commIdx]);
          MPI_Comm_group(commSendToBelow_[commIdx],&commGroup);
          //now for each supernode, we need to store the pointer to the communnicator and the rank of the root
          std::vector<Int> & snodeList = it->second;
          for(int curi = 0; curi < snodeList.size(); curi++){
            commSendToBelowPtr_[snodeList[curi]] = &commSendToBelow_[commIdx];
            Int ksup = snodeList[curi];

            Int curRoot = PROW(snodeList[curi],grid_);

            Int newRank = -1;
            if(color>0){
              MPI_Group_translate_ranks(colCommGroup, 1,&curRoot,commGroup, &newRank);
              if(newRank==MPI_UNDEFINED){
                statusOFS<<"["<<ksup<<"] Root ROW "<<curRoot<<" has no matching rank"<<std::endl;
              }
            }

            commSendToBelowRoot_[snodeList[curi]] = newRank;
          }

          commIdx++;

        }
      }

#ifdef CLEAR_MASKS
      //reduce memory usage
      maskSendToBelow_.clear();
#endif

    }
#ifndef _RELEASE_
    PopCallStack();
#endif







#ifdef PRINT_COMMUNICATOR_STAT
    {
      statusOFS << std::endl << "countSendToRight:" << countSendToRight_ << std::endl;
      statusOFS << std::endl << "maskSendToRight_:" << maskSendToRight_.size() <<std::endl; 
      bitMaskSet::iterator it;
      for(it = maskSendToRight_.begin(); it != maskSendToRight_.end(); it++){
        //print the involved processors
        for(int curi = 0; curi < it->first.size(); curi++){
          statusOFS << it->first[curi] << " "; 
        }

        statusOFS<< "    ( ";
        //print the involved supernode indexes
        std::vector<Int> & snodeList = it->second;
        for(int curi = 0; curi < snodeList.size(); curi++){
          statusOFS<<snodeList[curi]<<" ";
        }

        statusOFS << ")"<< std::endl;
      }
    }
#endif

#ifdef PRINT_COMMUNICATOR_STAT
    {
      statusOFS << std::endl << "countRecvFromBelow:" << countRecvFromBelow_ << std::endl;
      statusOFS << std::endl << "maskRecvFromBelow_:" << maskRecvFromBelow_.size() <<std::endl; 
      bitMaskSet::iterator it;
      for(it = maskRecvFromBelow_.begin(); it != maskRecvFromBelow_.end(); it++){
        //print the involved processors
        for(int curi = 0; curi < it->first.size(); curi++){
          statusOFS << it->first[curi] << " "; 
        }

        statusOFS<< "    ( ";
        //print the involved supernode indexes
        std::vector<Int> & snodeList = it->second;
        for(int curi = 0; curi < snodeList.size(); curi++){
          statusOFS<<snodeList[curi]<<" ";
        }

        statusOFS << ")"<< std::endl;
      }


    }
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif





#ifndef _RELEASE_
    PushCallStack("Creating SendToRight communicator");
#endif
    {
      MPI_Group rowCommGroup;  
      MPI_Comm_group(grid_->rowComm,&rowCommGroup);

      commSendToRight_.resize(maskSendToRight_.size());
      commSendToRightRoot_.resize(numSuper);
      commSendToRightPtr_.resize(numSuper);

      Int commIdx = 0;
      bitMaskSet::iterator it;
      for(it = maskSendToRight_.begin(); it != maskSendToRight_.end(); it++){
        Int count= std::count(it->first.begin(), it->first.end(), true);

        if(count>1){
          MPI_Group commGroup;  

          Int color = it->first[MYCOL(grid_)];

          MPI_Comm_split(grid_->rowComm, color  ,MYCOL(grid_) , &commSendToRight_[commIdx]);
          MPI_Comm_group(commSendToRight_[commIdx],&commGroup);
          //now for each supernode, we need to store the pointer to the communnicator and the rank of the root
          std::vector<Int> & snodeList = it->second;
          for(int curi = 0; curi < snodeList.size(); curi++){
            commSendToRightPtr_[snodeList[curi]] = &commSendToRight_[commIdx];
            Int ksup = snodeList[curi];
            Int curRoot = PCOL(snodeList[curi],grid_);
            Int newRank = -1;
            if(color>0){
              MPI_Group_translate_ranks(rowCommGroup, 1,&curRoot,commGroup, &newRank);

              if(newRank==MPI_UNDEFINED){
                statusOFS<<"["<<ksup<<"] Root COL "<<curRoot<<" has no matching rank"<<std::endl;
              }
            }

            commSendToRightRoot_[snodeList[curi]] = newRank;
          }

          commIdx++;
        }
      }

#ifdef CLEAR_MASKS
      //reduce memory usage
      maskSendToRight_.clear();
#endif
    }
#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack("Creating RecvFromBelow communicator");
#endif
    {
      MPI_Group colCommGroup;  
      MPI_Comm_group(grid_->colComm,&colCommGroup);

      commRecvFromBelow_.resize(maskRecvFromBelow_.size());
      commRecvFromBelowRoot_.resize(numSuper);
      commRecvFromBelowPtr_.resize(numSuper);

      Int commIdx = 0;
      bitMaskSet::iterator it;
      for(it = maskRecvFromBelow_.begin(); it != maskRecvFromBelow_.end(); it++){
        Int count= std::count(it->first.begin(), it->first.end(), true);

        if(count>1){
          MPI_Group commGroup;  

          Int color = it->first[MYROW(grid_)];

          MPI_Comm_split(grid_->colComm, color  ,MYROW(grid_) , &commRecvFromBelow_[commIdx]);
          MPI_Comm_group(commRecvFromBelow_[commIdx],&commGroup);
          //now for each supernode, we need to store the pointer to the communnicator and the rank of the root
          std::vector<Int> & snodeList = it->second;
          for(int curi = 0; curi < snodeList.size(); curi++){
            commRecvFromBelowPtr_[snodeList[curi]] = &commRecvFromBelow_[commIdx];
            Int ksup = snodeList[curi];
            Int curRoot = PROW(snodeList[curi],grid_);
            Int newRank = -1;
            if(color>0){
              MPI_Group_translate_ranks(colCommGroup, 1,&curRoot,commGroup, &newRank);

              if(newRank==MPI_UNDEFINED){
                statusOFS<<"["<<ksup<<"] Root ROW "<<curRoot<<" has no matching rank"<<std::endl;
              }
            }

            commRecvFromBelowRoot_[snodeList[curi]] = newRank;
          }

          commIdx++;
        }
      }


#ifdef CLEAR_MASKS
      //reduce memory usage
      maskRecvFromBelow_.clear();
#endif

    }
#ifndef _RELEASE_
    PopCallStack();
#endif







#ifndef _RELEASE_
    PushCallStack("SendToCrossDiagonal / RecvFromCrossDiagonal");
#endif


    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end(); si++ ){
          Int isup = *si;
          Int isupProcRow = PROW( isup, grid_ );
          Int isupProcCol = PCOL( isup, grid_ );
          if( isup > ksup ){
            if( MYROW( grid_ ) == isupProcRow){ 
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<"["<<ksup<<"] should send "<<isup<<" to "<<PNUM(PROW(ksup,grid_),isupProcCol,grid_)<<std::endl;
#endif
              isSendToCrossDiagonal_(grid_->numProcCol, ksup ) = true;
              isSendToCrossDiagonal_(isupProcCol, ksup ) = true;
            }
          }
        } // for (si)
      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )

      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localRowBlockColIdx[ LBi(ksup, grid_) ].begin();
            si != localRowBlockColIdx[ LBi(ksup, grid_) ].end(); si++ ){
          Int jsup = *si;
          Int jsupProcCol = PCOL( jsup, grid_ );
          Int jsupProcRow = PROW( jsup, grid_ );
          if( jsup > ksup){
            if( MYCOL(grid_) == jsupProcCol ){
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<"["<<ksup<<"] should receive "<<jsup<<" from "<<PNUM(jsupProcRow,PCOL(ksup,grid_),grid_)<<std::endl;
#endif
              isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) = true;
              isRecvFromCrossDiagonal_(jsupProcRow, ksup ) = true;
            }
          }
        } // for (si)
      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )

    } // for (ksup)




#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToCrossDiagonal:" << std::endl;
    for(int j =0; j < isSendToCrossDiagonal_.n();j++){
      if(isSendToCrossDiagonal_(grid_->numProcCol,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isSendToCrossDiagonal_.m()-1;i++){
          if(isSendToCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(PROW(j,grid_),i,grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }

    statusOFS << std::endl << "isRecvFromCrossDiagonal:" << std::endl;
    for(int j =0; j < isRecvFromCrossDiagonal_.n();j++){
      if(isRecvFromCrossDiagonal_(grid_->numProcRow,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isRecvFromCrossDiagonal_.m()-1;i++){
          if(isRecvFromCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(i,PCOL(j,grid_),grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }
#endif



#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PopCallStack();
#endif

    //Build the list of supernodes based on the elimination tree from SuperLU
    std::vector<std::vector<Int> > & WSet = this->WorkingSet();


    if (options_->maxPipelineDepth!=1){

      //find roots in the supernode etree (it must be postordered)
      //initialize the parent we are looking at 
      Int rootParent = snodeEtree[this->NumSuper()-2];

      //look for roots in the forest
      std::vector< Int>  initialRootList(1,rootParent);
      std::vector< Int>  mergeRootBuf;
      Int prevRootIdx = -1;
      std::vector< Int> & prevRoot = initialRootList;

      /* initialize the num of child for each node */
      Int nsupers = this->NumSuper();
      std::vector<Int> num_child;
      num_child.resize(nsupers,0);
      for(Int i=0; i<nsupers; i++ ) if( snodeEtree[i] != nsupers ) num_child[snodeEtree[i]] ++;

      while(prevRoot.size()>0){
        WSet.push_back(std::vector<Int>());
        Int totalChild =0;
        for(Int i = 0; i<prevRoot.size();++i){ totalChild += num_child[prevRoot[i]]; }
        WSet.back().reserve(totalChild);

        for(Int i = 0; i<prevRoot.size();++i){
          rootParent = prevRoot[i];
          std::vector<Int>::iterator parentIt = snodeEtree.begin()+rootParent;
          std::vector<Int>::iterator curRootIt = std::find (snodeEtree.begin() ,parentIt, rootParent);
          while(curRootIt != parentIt){
            Int curNode = curRootIt - snodeEtree.begin();
            WSet.back().push_back(curNode);
            //switch the sign to remove this root
            *curRootIt =-*curRootIt;
            //look for next root
            curRootIt = std::find (snodeEtree.begin() ,parentIt, rootParent);
          }
        }
        //No we have now several roots >> must maintain a vector of roots

          prevRootIdx++;
          prevRoot = WSet[prevRootIdx];
      }
      if(WSet.back().size()==0){
        WSet.pop_back();
      }


      for (Int lidx=0; lidx<WSet.size() ; lidx++){
        if(options_->maxPipelineDepth){
          if(WSet[lidx].size()>options_->maxPipelineDepth)
          {
            std::vector<std::vector<Int> >::iterator pos = WSet.begin()+lidx+1;               
            WSet.insert(pos,std::vector<Int>());
            WSet[lidx+1].insert(WSet[lidx+1].begin(),WSet[lidx].begin() +options_->maxPipelineDepth ,WSet[lidx].end());
            WSet[lidx].erase(WSet[lidx].begin()+options_->maxPipelineDepth,WSet[lidx].end());
          }
        }
      }

#if ( _DEBUGlevel_ >= 1 )
      for (Int lidx=0; lidx<WSet.size() ; lidx++){
        statusOFS << std::endl << "L"<< lidx << " is: {";
        for (Int supidx=0; supidx<WSet[lidx].size() ; supidx++){
          statusOFS << WSet[lidx][supidx] << " ["<<snodeEtree[WSet[lidx][supidx]]<<"] ";
        }
        statusOFS << " }"<< std::endl;
      }
#endif


    }
    else{
      for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
        WSet.push_back(std::vector<Int>());
        WSet.back().push_back(ksup);
      }
#if ( _DEBUGlevel_ >= 1 )
      for (Int lidx=0; lidx<WSet.size() ; lidx++){
        statusOFS << std::endl << "L"<< lidx << " is: {";
        for (Int supidx=0; supidx<WSet[lidx].size() ; supidx++){
          statusOFS << WSet[lidx][supidx] << " ";
        }
        statusOFS << " }"<< std::endl;
      }
#endif



    }


    return ;
  } 		// -----  end of method PMatrix::ConstructCommunicationPattern_Collectives  ----- 
#else
  void PMatrix::ConstructCommunicationPattern_Collectives	(  )
  {
#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
//    TAU_PROFILE_SET_CONTEXT(grid_->comm);
#endif


    TIMER_START(ConstructCommunicationPattern_Collectives);
#ifndef _RELEASE_
    PushCallStack("PMatrix::ConstructCommunicationPattern_Collectives");
#endif
    Int numSuper = this->NumSuper();
#ifndef _RELEASE_
    PushCallStack( "Initialize the communication pattern" );
#endif
    isSendToBelow_.Resize(grid_->numProcRow, numSuper);
    isSendToRight_.Resize(grid_->numProcCol, numSuper);
    SetValue( isSendToBelow_, false );
    SetValue( isSendToRight_, false );


#ifdef HYBRID
    isRecvFromBelow_.Resize(grid_->numProcRow, numSuper);
    SetValue( isRecvFromBelow_, false );
    isSendToDiagonal_.Resize( numSuper );
    SetValue( isSendToDiagonal_, false );
#endif

    isRecvFromAbove_.Resize( numSuper );
    isRecvFromLeft_.Resize( numSuper );
    SetValue( isRecvFromAbove_, false );
    SetValue( isRecvFromLeft_, false );

    isSendToCrossDiagonal_.Resize(grid_->numProcCol+1, numSuper );
    SetValue( isSendToCrossDiagonal_, false );
    isRecvFromCrossDiagonal_.Resize(grid_->numProcRow+1, numSuper );
    SetValue( isRecvFromCrossDiagonal_, false );

    countSendToBelow_.Resize(numSuper);
    countSendToRight_.Resize(numSuper);
    countRecvFromBelow_.Resize( numSuper );
    SetValue( countSendToBelow_, 0 );
    SetValue( countSendToRight_, 0 );
    SetValue( countRecvFromBelow_, 0 );


    
    commSendToRightMask_.resize(numSuper);
    commSendToBelowMask_.resize(numSuper);
    commRecvFromBelowMask_.resize(numSuper);


#ifndef _RELEASE_
    PopCallStack();
#endif



    TIMER_START(Get_ETREE);

      std::vector<Int> snodeEtree(this->NumSuper());
      GetEtree(snodeEtree);
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << " Supernodal elimination tree is " << snodeEtree <<std::endl<<std::endl;
#endif

    TIMER_STOP(Get_ETREE);








#ifndef _RELEASE_
    PushCallStack( "Local column communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local column communication" << std::endl;
#endif


    TIMER_START(Local_Column_Communication);

    // localColBlockRowIdx stores the nonzero block indices for each local block column.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockCol x numNonzeroBlock
    std::vector<std::set<Int> >   localColBlockRowIdx;
    localColBlockRowIdx.resize( this->NumLocalBlockCol() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        std::vector<Int>  tBlockRowIdx;
        tBlockRowIdx.clear();

        // L part
        std::vector<LBlock>& Lcol = this->L( LBj(ksup, grid_) );
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          tBlockRowIdx.push_back( Lcol[ib].blockIdx );
        }

        // U part
        for( Int ib = 0; ib < this->NumLocalBlockRow(); ib++ ){
          std::vector<UBlock>& Urow = this->U(ib);
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            if( Urow[jb].blockIdx == ksup ){
              tBlockRowIdx.push_back( GBi( ib, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockRowIdx;
        mpi::Allgatherv( tBlockRowIdx, tAllBlockRowIdx, grid_->colComm );

        localColBlockRowIdx[LBj( ksup, grid_ )].insert(
            tAllBlockRowIdx.begin(), tAllBlockRowIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Column block " << ksup 
          << " has the following nonzero block rows" << std::endl;
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for(ksup)

    TIMER_STOP(Local_Column_Communication);

#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack( "Local row communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local row communication" << std::endl;
#endif

    TIMER_START(Local_Row_Communication);

    // localRowBlockColIdx stores the nonzero block indices for each local block row.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockRow x numNonzeroBlock
    std::vector<std::set<Int> >   localRowBlockColIdx;

    localRowBlockColIdx.resize( this->NumLocalBlockRow() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        std::vector<Int>  tBlockColIdx;
        tBlockColIdx.clear();

        // U part
        std::vector<UBlock>& Urow = this->U( LBi(ksup, grid_) );
        for( Int jb = 0; jb < Urow.size(); jb++ ){
          tBlockColIdx.push_back( Urow[jb].blockIdx );
        }

        // L part
        for( Int jb = 0; jb < this->NumLocalBlockCol(); jb++ ){
          std::vector<LBlock>& Lcol = this->L(jb);
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            if( Lcol[ib].blockIdx == ksup ){
              tBlockColIdx.push_back( GBj( jb, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockColIdx;
        mpi::Allgatherv( tBlockColIdx, tAllBlockColIdx, grid_->rowComm );

        localRowBlockColIdx[LBi( ksup, grid_ )].insert(
            tAllBlockColIdx.begin(), tAllBlockColIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Row block " << ksup 
          << " has the following nonzero block columns" << std::endl;
        for( std::set<Int>::iterator si = localRowBlockColIdx[LBi( ksup, grid_ )].begin();
            si != localRowBlockColIdx[LBi( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
    } // for(ksup)

    TIMER_STOP(Local_Row_Communication);
#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack("SendToBelow / RecvFromAbove");
#endif

    TIMER_START(SendToBelow_RecvFromAbove);

    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes to the right of ksup

      std::vector<bool> sTB(grid_->numProcRow,false);

      //Explore the ancestors of ksup
      Int jsup = snodeEtree[ksup];
      while(jsup<numSuper)
      {
        Int jsupLocalBlockCol = LBj( jsup, grid_ );
        Int jsupProcCol = PCOL( jsup, grid_ );
        if( MYCOL( grid_ ) == jsupProcCol ){

          // SendToBelow / RecvFromAbove only if (ksup, jsup) is nonzero.
          if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ) {
            for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
                si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
              Int isup = *si;
              Int isupProcRow = PROW( isup, grid_ );
              if( isup > ksup ){
                if( MYROW( grid_ ) == isupProcRow ){
                  isRecvFromAbove_(ksup) = true;
                }
                if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
                  isSendToBelow_( isupProcRow, ksup ) = true;
                }
              } // if( isup > ksup )
            } // for (si)
          } // if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 )

          sTB[ PROW(ksup,grid_) ] = true;
          if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ) {
            for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
                si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
              Int isup = *si;
              Int isupProcRow = PROW( isup, grid_ );
              if( isup > ksup ){
                sTB[isupProcRow] = true;
              } // if( isup > ksup )
            } // for (si)
          }
        } // if( MYCOL( grid_ ) == PCOL( jsup, grid_ ) )

        jsup = snodeEtree[jsup];
      } // while(jsup)

#if ( _DEBUGlevel_ >= 1 )
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        statusOFS<<"["<<ksup<<"] sTB is: "; for(int curi=0; curi<sTB.size();curi++){statusOFS<<sTB[curi]<<" ";} statusOFS<<std::endl;
        statusOFS<<"["<<ksup<<"] isSendToBelow_ is: "; for(int curi=0; curi<sTB.size();curi++){statusOFS<<isSendToBelow_(curi,ksup)<<" ";} statusOFS<<std::endl;
      }
#endif
    TIMER_START(Creating_MaskSendToBelow);
      Int count= std::count(sTB.begin(), sTB.end(), true);
      Int color = sTB[MYROW(grid_)];
      if(count>1){
      commSendToBelowMask_[ksup] = sTB;
      }
      countSendToBelow_(ksup) = count * color;
    TIMER_STOP(Creating_MaskSendToBelow);
    } // for(ksup)

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToBelow:" << std::endl;
    for(int j = 0;j< isSendToBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToBelow_.m();i++){
        statusOFS<< isSendToBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromAbove:" << std::endl;
    for(int j = 0;j< isRecvFromAbove_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromAbove_(j)<<std::endl;
    }
#endif
#ifdef PRINT_COMMUNICATOR_STAT
    {
      statusOFS << std::endl << "countSendToBelow:" << countSendToBelow_ << std::endl;
    }
#endif
#ifndef _RELEASE_
    PopCallStack();
#endif



    TIMER_STOP(SendToBelow_RecvFromAbove);


#ifndef _RELEASE_
    PushCallStack("SendToRight / RecvFromLeft");
#endif


    TIMER_START(SendToRight_RecvFromLeft);


    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes below ksup
      std::vector<bool> sTR(grid_->numProcCol,false);
      std::vector<bool> rFB(grid_->numProcRow,false);

      Int isup = snodeEtree[ksup];
      while(isup<numSuper)
      {
        Int isupLocalBlockRow = LBi( isup, grid_ );
        Int isupProcRow       = PROW( isup, grid_ );
        if( MYROW( grid_ ) == isupProcRow ){
          // SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.
          if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
            for( std::set<Int>::iterator sj = localRowBlockColIdx[isupLocalBlockRow].begin();
                sj != localRowBlockColIdx[isupLocalBlockRow].end(); sj++ ){
              Int jsup = *sj;
              Int jsupProcCol = PCOL( jsup, grid_ );
              if( jsup > ksup ){

                if( MYCOL( grid_ ) == jsupProcCol ){
                  isRecvFromLeft_(ksup) = true;
                }
                if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
                  isSendToRight_( jsupProcCol, ksup ) = true;
                }
              }
            } // for (sj)
          } // if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 )


          sTR[ PCOL(ksup,grid_) ] = true;
          if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
            for( std::set<Int>::iterator sj = localRowBlockColIdx[isupLocalBlockRow].begin();
                sj != localRowBlockColIdx[isupLocalBlockRow].end(); sj++ ){
              Int jsup = *sj;
              Int jsupProcCol = PCOL( jsup, grid_ );
              if( jsup > ksup ){
                sTR[ jsupProcCol ] = true;
              } // if( jsup > ksup )
            } // for (sj)
          }
        } // if( MYROW( grid_ ) == isupProcRow )


        rFB[ PROW(ksup,grid_) ] = true;
        if( MYCOL( grid_ ) == PCOL(ksup, grid_) ){
          rFB[ isupProcRow ] = true;

#ifdef HYBRID
          if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){ 
            isRecvFromBelow_(isupProcRow,ksup) = true;
          }    
          else if (MYROW(grid_) == isupProcRow){
            isSendToDiagonal_(ksup)=true;
          }
#endif
        } // if( MYCOL( grid_ ) == PCOL(ksup, grid_) )

        isup = snodeEtree[isup];
      } // while(isup)
#if ( _DEBUGlevel_ >= 1 )
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        statusOFS<<"["<<ksup<<"] sTR is: "; for(int curi=0; curi<sTR.size();curi++){statusOFS<<sTR[curi]<<" ";} statusOFS<<std::endl;
        statusOFS<<"["<<ksup<<"] rFB is: "; for(int curi=0; curi<rFB.size();curi++){statusOFS<<rFB[curi]<<" ";} statusOFS<<std::endl;
      }
#endif
      //      std::vector<bool> mask( sTR.Data(), sTR.Data() + sTR.m() );
    TIMER_START(Creating_MaskSendToRight_RecvFromBelow);
      Int count= std::count(sTR.begin(), sTR.end(), true);
      Int color = sTR[MYCOL(grid_)];
      if(count>1){
        commSendToRightMask_[ksup] = sTR;
      }
      countSendToRight_(ksup) = count * color;



      count= std::count(rFB.begin(), rFB.end(), true);
      color = rFB[MYROW(grid_)];
      if(count>1){
        commRecvFromBelowMask_[ksup] = rFB;
      }
      countRecvFromBelow_(ksup) = count * color;
    TIMER_STOP(Creating_MaskSendToRight_RecvFromBelow);

    }	 // for (ksup)


#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToRight:" << std::endl;
    for(int j = 0;j< isSendToRight_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToRight_.m();i++){
        statusOFS<< isSendToRight_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromLeft:" << std::endl;
    for(int j = 0;j< isRecvFromLeft_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromLeft_(j)<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromBelow:" << std::endl;
    for(int j = 0;j< isRecvFromBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isRecvFromBelow_.m();i++){
        statusOFS<< isRecvFromBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }
#endif

#ifdef PRINT_COMMUNICATOR_STAT
    {
      statusOFS << std::endl << "countSendToRight:" << countSendToRight_ << std::endl;


      statusOFS << std::endl << "countRecvFromBelow:" << countRecvFromBelow_ << std::endl;
    }
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif





    TIMER_STOP(SendToRight_RecvFromLeft);



#ifndef _RELEASE_
    PushCallStack("SendToCrossDiagonal / RecvFromCrossDiagonal");
#endif

    TIMER_START(SendToCD);

    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end(); si++ ){
          Int isup = *si;
          Int isupProcRow = PROW( isup, grid_ );
          Int isupProcCol = PCOL( isup, grid_ );
          if( isup > ksup ){
            if( MYROW( grid_ ) == isupProcRow){ 
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<"["<<ksup<<"] should send "<<isup<<" to "<<PNUM(PROW(ksup,grid_),isupProcCol,grid_)<<std::endl;
#endif
              isSendToCrossDiagonal_(grid_->numProcCol, ksup ) = true;
              isSendToCrossDiagonal_(isupProcCol, ksup ) = true;
            }
          }
        } // for (si)
      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )

      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localRowBlockColIdx[ LBi(ksup, grid_) ].begin();
            si != localRowBlockColIdx[ LBi(ksup, grid_) ].end(); si++ ){
          Int jsup = *si;
          Int jsupProcCol = PCOL( jsup, grid_ );
          Int jsupProcRow = PROW( jsup, grid_ );
          if( jsup > ksup){
            if( MYCOL(grid_) == jsupProcCol ){
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<"["<<ksup<<"] should receive "<<jsup<<" from "<<PNUM(jsupProcRow,PCOL(ksup,grid_),grid_)<<std::endl;
#endif
              isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) = true;
              isRecvFromCrossDiagonal_(jsupProcRow, ksup ) = true;
            }
          }
        } // for (si)
      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )

    } // for (ksup)




#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToCrossDiagonal:" << std::endl;
    for(int j =0; j < isSendToCrossDiagonal_.n();j++){
      if(isSendToCrossDiagonal_(grid_->numProcCol,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isSendToCrossDiagonal_.m()-1;i++){
          if(isSendToCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(PROW(j,grid_),i,grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }

    statusOFS << std::endl << "isRecvFromCrossDiagonal:" << std::endl;
    for(int j =0; j < isRecvFromCrossDiagonal_.n();j++){
      if(isRecvFromCrossDiagonal_(grid_->numProcRow,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isRecvFromCrossDiagonal_.m()-1;i++){
          if(isRecvFromCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(i,PCOL(j,grid_),grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }
#endif


    TIMER_STOP(SendToCD);

#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PopCallStack();
#endif



    //Build the list of supernodes based on the elimination tree from SuperLU
    GetWorkSet(snodeEtree,this->WorkingSet());

    TIMER_STOP(ConstructCommunicationPattern_Collectives);
    return ;
  } 		// -----  end of method PMatrix::ConstructCommunicationPattern_Collectives  ----- 


#endif


  void PMatrix::getMaxCommunicatorSizes(){
    IntNumVec maxSizes(this->WorkingSet().size());
    SetValue(maxSizes,I_ONE);
    for (Int lidx=0; lidx<this->WorkingSet().size() ; lidx++){
      Int stepSuper = this->WorkingSet()[lidx].size(); 
      for(Int supidx = 0;supidx<stepSuper;supidx++){
        Int ksup = this->WorkingSet()[lidx][supidx];

        Int count= std::count(commSendToRightMask_[ksup].begin(), commSendToRightMask_[ksup].end(), true);
        maxSizes[lidx] = std::max(maxSizes[lidx],count);

        count= std::count(commSendToBelowMask_[ksup].begin(), commSendToBelowMask_[ksup].end(), true);
        maxSizes[lidx] = std::max(maxSizes[lidx],count);

        count= std::count(commRecvFromBelowMask_[ksup].begin(), commRecvFromBelowMask_[ksup].end(), true);
        maxSizes[lidx] = std::max(maxSizes[lidx],count);
      }
    }


    maxCommSizes_.Resize(maxSizes.m());
    SetValue(maxCommSizes_,I_ONE);
    MPI_Allreduce( maxSizes.Data(), maxCommSizes_.Data(), maxSizes.m(), MPI_INT, MPI_MAX, grid_->comm );
  }






inline  void PMatrix::SelInv_lookup_indexes(const Int ksup, std::vector<LBlock> & LcolRecv,  std::vector<UBlock> & UrowRecv, NumMat<Scalar> & AinvBuf,NumMat<Scalar> & UBuf,NumMat<Scalar> & LUpdateBuf){
    TIMER_START(Compute_Sinv_LT_Lookup_Indexes);

    // rowPtr[ib] gives the row index in LUpdateBuf for the first
    // nonzero row in LcolRecv[ib]. The total number of rows in
    // LUpdateBuf is given by rowPtr[end]-1
    std::vector<Int> rowPtr(LcolRecv.size() + 1);
    // colPtr[jb] gives the column index in UBuf for the first
    // nonzero column in UrowRecv[jb]. The total number of rows in
    // UBuf is given by colPtr[end]-1
    std::vector<Int> colPtr(UrowRecv.size() + 1);

    rowPtr[0] = 0;
    for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
      rowPtr[ib+1] = rowPtr[ib] + LcolRecv[ib].numRow;
    }
    colPtr[0] = 0;
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      colPtr[jb+1] = colPtr[jb] + UrowRecv[jb].numCol;
    }

    Int numRowAinvBuf = *rowPtr.rbegin();
    Int numColAinvBuf = *colPtr.rbegin();

    // Allocate for the computational storage
    AinvBuf.Resize( numRowAinvBuf, numColAinvBuf );

    LUpdateBuf.Resize( numRowAinvBuf, SuperSize( ksup, super_ ) );
    UBuf.Resize( SuperSize( ksup, super_ ), numColAinvBuf );
//    SetValue( AinvBuf, SCALAR_ZERO );
//    SetValue( LUpdateBuf, SCALAR_ZERO );
//    SetValue( UBuf, SCALAR_ZERO );

    // Fill UBuf first.  Make the transpose later in the Gemm phase.
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      UBlock& UB = UrowRecv[jb];
      if( UB.numRow != SuperSize(ksup, super_) ){
        throw std::logic_error( "The size of UB is not right.  Something is seriously wrong." );
      }
      lapack::Lacpy( 'A', UB.numRow, UB.numCol, UB.nzval.Data(),
          UB.numRow, UBuf.VecData( colPtr[jb] ), UBuf.m() );	
    }

    // Calculate the relative indices for (isup, jsup)
    // Fill AinvBuf with the information in L or U block.
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
        LBlock& LB = LcolRecv[ib];
        UBlock& UB = UrowRecv[jb];
        Int isup = LB.blockIdx;
        Int jsup = UB.blockIdx;
        Scalar* nzvalAinv = &AinvBuf( rowPtr[ib], colPtr[jb] );
        Int     ldAinv    = AinvBuf.m();

        // Pin down the corresponding block in the part of Sinv.
        if( isup >= jsup ){
          std::vector<LBlock>&  LcolSinv = this->L( LBj(jsup, grid_ ) );
          bool isBlockFound = false;
          for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( LcolSinv[ibSinv].blockIdx == isup ){
              LBlock& SinvB = LcolSinv[ibSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int* rowsLBPtr    = LB.rows.Data();
              Int* rowsSinvBPtr = SinvB.rows.Data();
              for( Int i = 0; i < LB.numRow; i++ ){
                bool isRowFound = false;
                for( Int i1 = 0; i1 < SinvB.numRow; i1++ ){
                  if( rowsLBPtr[i] == rowsSinvBPtr[i1] ){
                    isRowFound = true;
                    relRows[i] = i1;
                    break;
                  }
                }
                if( isRowFound == false ){
                  std::ostringstream msg;
                  msg << "Row " << rowsLBPtr[i] << 
                    " in LB cannot find the corresponding row in SinvB" << std::endl
                    << "LB.rows    = " << LB.rows << std::endl
                    << "SinvB.rows = " << SinvB.rows << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }

              // Column relative indicies
              std::vector<Int> relCols( UB.numCol );
              Int SinvColsSta = FirstBlockCol( jsup, super_ );
              for( Int j = 0; j < UB.numCol; j++ ){
                relCols[j] = UB.cols[j] - SinvColsSta;
              }

              // Transfer the values from Sinv to AinvBlock
              Scalar* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < UB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }

              isBlockFound = true;
              break;
            }	
          } // for (ibSinv )
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup 
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in L
        else{
          std::vector<UBlock>&   UrowSinv = this->U( LBi( isup, grid_ ) );
          bool isBlockFound = false;
          for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( UrowSinv[jbSinv].blockIdx == jsup ){
              UBlock& SinvB = UrowSinv[jbSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int SinvRowsSta = FirstBlockCol( isup, super_ );
              for( Int i = 0; i < LB.numRow; i++ ){
                relRows[i] = LB.rows[i] - SinvRowsSta;
              }

              // Column relative indices
              std::vector<Int> relCols( UB.numCol );
              Int* colsUBPtr    = UB.cols.Data();
              Int* colsSinvBPtr = SinvB.cols.Data();
              for( Int j = 0; j < UB.numCol; j++ ){
                bool isColFound = false;
                for( Int j1 = 0; j1 < SinvB.numCol; j1++ ){
                  if( colsUBPtr[j] == colsSinvBPtr[j1] ){
                    isColFound = true;
                    relCols[j] = j1;
                    break;
                  }
                }
                if( isColFound == false ){
                  std::ostringstream msg;
                  msg << "Col " << colsUBPtr[j] << 
                    " in UB cannot find the corresponding row in SinvB" << std::endl
                    << "UB.cols    = " << UB.cols << std::endl
                    << "UinvB.cols = " << SinvB.cols << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }


              // Transfer the values from Sinv to AinvBlock
              Scalar* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < UB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }

              isBlockFound = true;
              break;
            }
          } // for (jbSinv)
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup 
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in U

      } // for( ib )
    } // for ( jb )


    TIMER_STOP(Compute_Sinv_LT_Lookup_Indexes);
  }

inline  void PMatrix::SelInv_lookup_indexes(SuperNodeBufferType & snode, std::vector<LBlock> & LcolRecv,  std::vector<UBlock> & UrowRecv, NumMat<Scalar> & AinvBuf,NumMat<Scalar> & UBuf){
    TIMER_START(Compute_Sinv_LT_Lookup_Indexes);

    TIMER_START(Build_colptr_rowptr);
    // rowPtr[ib] gives the row index in snode.LUpdateBuf for the first
    // nonzero row in LcolRecv[ib]. The total number of rows in
    // snode.LUpdateBuf is given by rowPtr[end]-1
    std::vector<Int> rowPtr(LcolRecv.size() + 1);
    // colPtr[jb] gives the column index in UBuf for the first
    // nonzero column in UrowRecv[jb]. The total number of rows in
    // UBuf is given by colPtr[end]-1
    std::vector<Int> colPtr(UrowRecv.size() + 1);

    rowPtr[0] = 0;
    for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
      rowPtr[ib+1] = rowPtr[ib] + LcolRecv[ib].numRow;
    }
    colPtr[0] = 0;
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      colPtr[jb+1] = colPtr[jb] + UrowRecv[jb].numCol;
    }

    Int numRowAinvBuf = *rowPtr.rbegin();
    Int numColAinvBuf = *colPtr.rbegin();
    TIMER_STOP(Build_colptr_rowptr);

    TIMER_START(Allocate_lookup);
    // Allocate for the computational storage
    AinvBuf.Resize( numRowAinvBuf, numColAinvBuf );
    UBuf.Resize( SuperSize( snode.Index, super_ ), numColAinvBuf );
//    TIMER_START(SetValue_lookup);
//    SetValue( AinvBuf, SCALAR_ZERO );
    //SetValue( snode.LUpdateBuf, SCALAR_ZERO );
//    SetValue( UBuf, SCALAR_ZERO );
//    TIMER_STOP(SetValue_lookup);
    TIMER_STOP(Allocate_lookup);

    TIMER_START(Fill_UBuf);
    // Fill UBuf first.  Make the transpose later in the Gemm phase.
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      UBlock& UB = UrowRecv[jb];
      if( UB.numRow != SuperSize(snode.Index, super_) ){
        throw std::logic_error( "The size of UB is not right.  Something is seriously wrong." );
      }
      lapack::Lacpy( 'A', UB.numRow, UB.numCol, UB.nzval.Data(),
          UB.numRow, UBuf.VecData( colPtr[jb] ), SuperSize( snode.Index, super_ ) );	
    }
    TIMER_STOP(Fill_UBuf);

    // Calculate the relative indices for (isup, jsup)
    // Fill AinvBuf with the information in L or U block.
    TIMER_START(JB_Loop);

#ifdef STDFIND
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){

        UBlock& UB = UrowRecv[jb];
        Int jsup = UB.blockIdx;
        Int SinvColsSta = FirstBlockCol( jsup, super_ );

              // Column relative indicies
              std::vector<Int> relCols( UB.numCol );
              for( Int j = 0; j < UB.numCol; j++ ){
                relCols[j] = UB.cols[j] - SinvColsSta;
              }




      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
        LBlock& LB = LcolRecv[ib];
        Int isup = LB.blockIdx;
        Int SinvRowsSta = FirstBlockCol( isup, super_ );
        Scalar* nzvalAinv = &AinvBuf( rowPtr[ib], colPtr[jb] );
        Int     ldAinv    = numRowAinvBuf;

        // Pin down the corresponding block in the part of Sinv.
        if( isup >= jsup ){
          std::vector<LBlock>&  LcolSinv = this->L( LBj(jsup, grid_ ) );
          bool isBlockFound = false;
          for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( LcolSinv[ibSinv].blockIdx == isup ){
              LBlock& SinvB = LcolSinv[ibSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int* rowsLBPtr    = LB.rows.Data();
              Int* rowsSinvBPtr = SinvB.rows.Data();

              TIMER_START(STDFIND_ROW);
              Int * pos =&rowsSinvBPtr[0];
              Int * last =&rowsSinvBPtr[SinvB.numRow];
              for( Int i = 0; i < LB.numRow; i++ ){
//                pos = std::find(pos, &rowsSinvBPtr[SinvB.numRow-1], rowsLBPtr[i]);
                pos = std::find(rowsSinvBPtr, last, rowsLBPtr[i]);
                if(pos != last){
                  relRows[i] = (Int)(pos - rowsSinvBPtr);
                }
                else{
                  std::ostringstream msg;
                  msg << "Row " << rowsLBPtr[i] << 
                    " in LB cannot find the corresponding row in SinvB" << std::endl
                    << "LB.rows    = " << LB.rows << std::endl
                    << "SinvB.rows = " << SinvB.rows << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }
              TIMER_STOP(STDFIND_ROW);

              TIMER_START(Copy_Sinv_to_Ainv);
              // Transfer the values from Sinv to AinvBlock
              Scalar* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < UB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }
              TIMER_STOP(Copy_Sinv_to_Ainv);

              isBlockFound = true;
              break;
            }	
          } // for (ibSinv )
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup 
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in L
        else{
              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int SinvRowsSta = FirstBlockCol( isup, super_ );
              for( Int i = 0; i < LB.numRow; i++ ){
                relRows[i] = LB.rows[i] - SinvRowsSta;
              }
          std::vector<UBlock>&   UrowSinv = this->U( LBi( isup, grid_ ) );
          bool isBlockFound = false;
          for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( UrowSinv[jbSinv].blockIdx == jsup ){
              UBlock& SinvB = UrowSinv[jbSinv];

              

              // Column relative indices
              std::vector<Int> relCols( UB.numCol );
              Int* colsUBPtr    = UB.cols.Data();
              Int* colsSinvBPtr = SinvB.cols.Data();
              TIMER_START(STDFIND_COL);
              Int * pos =&colsSinvBPtr[0];
              Int * last =&colsSinvBPtr[SinvB.numCol];
              for( Int j = 0; j < UB.numCol; j++ ){
                //colsUB is sorted
                pos = std::find(colsSinvBPtr, last, colsUBPtr[j]);
                if(pos !=last){
                    relCols[j] = (Int)(pos - colsSinvBPtr);
                }
                else{
                  std::ostringstream msg;
                  msg << "Col " << colsUBPtr[j] << 
                    " in UB cannot find the corresponding row in SinvB" << std::endl
                    << "UB.cols    = " << UB.cols << std::endl
                    << "UinvB.cols = " << SinvB.cols << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }
              TIMER_STOP(STDFIND_COL);


              TIMER_START(Copy_Sinv_to_Ainv);
              // Transfer the values from Sinv to AinvBlock
              Scalar* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < UB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }
              TIMER_STOP(Copy_Sinv_to_Ainv);

              isBlockFound = true;
              break;
            }
          } // for (jbSinv)
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup 
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in U

      } // for( ib )
    } // for ( jb )
#else
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
        LBlock& LB = LcolRecv[ib];
        UBlock& UB = UrowRecv[jb];
        Int isup = LB.blockIdx;
        Int jsup = UB.blockIdx;
        Scalar* nzvalAinv = &AinvBuf( rowPtr[ib], colPtr[jb] );
        Int     ldAinv    = AinvBuf.m();

        // Pin down the corresponding block in the part of Sinv.
        if( isup >= jsup ){
          std::vector<LBlock>&  LcolSinv = this->L( LBj(jsup, grid_ ) );
          bool isBlockFound = false;
          for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( LcolSinv[ibSinv].blockIdx == isup ){
              LBlock& SinvB = LcolSinv[ibSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int* rowsLBPtr    = LB.rows.Data();
              Int* rowsSinvBPtr = SinvB.rows.Data();
              for( Int i = 0; i < LB.numRow; i++ ){
                bool isRowFound = false;
                for( Int i1 = 0; i1 < SinvB.numRow; i1++ ){
                  if( rowsLBPtr[i] == rowsSinvBPtr[i1] ){
                    isRowFound = true;
                    relRows[i] = i1;
                    break;
                  }
                }
                if( isRowFound == false ){
                  std::ostringstream msg;
                  msg << "Row " << rowsLBPtr[i] << 
                    " in LB cannot find the corresponding row in SinvB" << std::endl
                    << "LB.rows    = " << LB.rows << std::endl
                    << "SinvB.rows = " << SinvB.rows << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }

              // Column relative indicies
              std::vector<Int> relCols( UB.numCol );
              Int SinvColsSta = FirstBlockCol( jsup, super_ );
              for( Int j = 0; j < UB.numCol; j++ ){
                relCols[j] = UB.cols[j] - SinvColsSta;
              }

              // Transfer the values from Sinv to AinvBlock
              Scalar* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < UB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }

              isBlockFound = true;
              break;
            }	
          } // for (ibSinv )
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup 
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in L
        else{
          std::vector<UBlock>&   UrowSinv = this->U( LBi( isup, grid_ ) );
          bool isBlockFound = false;
          for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( UrowSinv[jbSinv].blockIdx == jsup ){
              UBlock& SinvB = UrowSinv[jbSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int SinvRowsSta = FirstBlockCol( isup, super_ );
              for( Int i = 0; i < LB.numRow; i++ ){
                relRows[i] = LB.rows[i] - SinvRowsSta;
              }

              // Column relative indices
              std::vector<Int> relCols( UB.numCol );
              Int* colsUBPtr    = UB.cols.Data();
              Int* colsSinvBPtr = SinvB.cols.Data();
              for( Int j = 0; j < UB.numCol; j++ ){
                bool isColFound = false;
                for( Int j1 = 0; j1 < SinvB.numCol; j1++ ){
                  if( colsUBPtr[j] == colsSinvBPtr[j1] ){
                    isColFound = true;
                    relCols[j] = j1;
                    break;
                  }
                }
                if( isColFound == false ){
                  std::ostringstream msg;
                  msg << "Col " << colsUBPtr[j] << 
                    " in UB cannot find the corresponding row in SinvB" << std::endl
                    << "UB.cols    = " << UB.cols << std::endl
                    << "UinvB.cols = " << SinvB.cols << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }


              // Transfer the values from Sinv to AinvBlock
              Scalar* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < UB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }

              isBlockFound = true;
              break;
            }
          } // for (jbSinv)
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup 
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in U

      } // for( ib )
    } // for ( jb )


#endif
    TIMER_STOP(JB_Loop);


    TIMER_STOP(Compute_Sinv_LT_Lookup_Indexes);
  }


      inline void PMatrix::SendRecvCD_UpdateU(std::vector<SuperNodeBufferType> & arrSuperNodes, Int stepSuper)
      {

      TIMER_START(Send_CD_Update_U);
      //compute the number of requests
      Int sendCount = 0;
      Int recvCount = 0;
      Int sendOffset[stepSuper];
      Int recvOffset[stepSuper];
      Int recvIdx=0;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        sendOffset[supidx]=sendCount;
        recvOffset[supidx]=recvCount;
        sendCount+= CountSendToCrossDiagonal(snode.Index);
        recvCount+= CountRecvFromCrossDiagonal(snode.Index);
      }


      std::vector<MPI_Request > arrMpiReqsSendCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeSendCD(sendCount, MPI_REQUEST_NULL );

      std::vector<MPI_Request > arrMpiReqsRecvCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeRecvCD(recvCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrLcolSendCD(sendCount);
      std::vector<int > arrSstrLcolSizeSendCD(sendCount);
      std::vector<std::vector<char> > arrSstrLcolRecvCD(recvCount);
      std::vector<int > arrSstrLcolSizeRecvCD(recvCount);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        // Send LUpdateBufReduced to the cross diagonal blocks. 
        // NOTE: This assumes square processor grid

        TIMER_START(Send_L_CrossDiag);

        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) && isSendToCrossDiagonal_(grid_->numProcCol, snode.Index ) ){

          Int sendIdx = 0;
          for(Int dstCol = 0; dstCol<grid_->numProcCol; dstCol++){
            if(isSendToCrossDiagonal_(dstCol,snode.Index) ){
              Int dest = PNUM(PROW(snode.Index,grid_),dstCol,grid_);

              if( MYPROC( grid_ ) != dest	){
                MPI_Request & mpiReqSizeSend = arrMpiReqsSizeSendCD[sendOffset[supidx]+sendIdx];
                MPI_Request & mpiReqSend = arrMpiReqsSendCD[sendOffset[supidx]+sendIdx];


                std::stringstream sstm;
                std::vector<char> & sstrLcolSend = arrSstrLcolSendCD[sendOffset[supidx]+sendIdx];
                Int & sstrSize = arrSstrLcolSizeSendCD[sendOffset[supidx]+sendIdx];

                serialize( snode.RowLocalPtr, sstm, NO_MASK );
                serialize( snode.BlockIdxLocal, sstm, NO_MASK );
                serialize( snode.LUpdateBuf, sstm, NO_MASK );

                sstrLcolSend.resize( Size(sstm) );
                sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
                sstrSize = sstrLcolSend.size();


                MPI_Isend( &sstrSize, 1, MPI_INT, dest, IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE), grid_->comm, &mpiReqSizeSend );
                MPI_Isend( (void*)&sstrLcolSend[0], sstrSize, MPI_BYTE, dest, IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT), grid_->comm, &mpiReqSend );
                sendIdx++;
              }
            }
          }


        } // sender
        TIMER_STOP(Send_L_CrossDiag);
      }


      //Do Irecv for sizes
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        //If I'm a receiver
        if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) && isRecvFromCrossDiagonal_(grid_->numProcRow, snode.Index ) ){
          Int recvIdx=0;
          for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
            if(isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,grid_),grid_);
              if( MYPROC( grid_ ) != src ){
                Int & sstrSize = arrSstrLcolSizeRecvCD[recvOffset[supidx]+recvIdx];
                MPI_Request & mpiReqSizeRecv = arrMpiReqsSizeRecvCD[recvOffset[supidx]+recvIdx];
                MPI_Irecv( &sstrSize, 1, MPI_INT, src, IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE), grid_->comm, &mpiReqSizeRecv );
                recvIdx++;
              }
            }
          }
        }//end if I'm a receiver
      }

      //waitall sizes
      mpi::Waitall(arrMpiReqsSizeRecvCD);
      //Allocate content and do Irecv
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        //If I'm a receiver
        if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) && isRecvFromCrossDiagonal_(grid_->numProcRow, snode.Index ) ){
          Int recvIdx=0;
          for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
            if(isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,grid_),grid_);
              if( MYPROC( grid_ ) != src ){
                Int & sstrSize = arrSstrLcolSizeRecvCD[recvOffset[supidx]+recvIdx];
                std::vector<char> & sstrLcolRecv = arrSstrLcolRecvCD[recvOffset[supidx]+recvIdx];
                MPI_Request & mpiReqRecv = arrMpiReqsRecvCD[recvOffset[supidx]+recvIdx];
                sstrLcolRecv.resize( sstrSize);
                MPI_Irecv( (void*)&sstrLcolRecv[0], sstrSize, MPI_BYTE, src, IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT), grid_->comm, &mpiReqRecv );
//statusOFS<<"P"<<MYPROC(grid_)<<" received "<<sstrSize<<" bytes of L/U from CD P"<<src<<std::endl;               
                recvIdx++;
              }
            }
          }
        }//end if I'm a receiver
      }

      //waitall content
      mpi::Waitall(arrMpiReqsRecvCD);
      //Do the work
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        // Send LUpdateBufReduced to the cross diagonal blocks. 
        // NOTE: This assumes square processor grid
        if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) && isRecvFromCrossDiagonal_(grid_->numProcRow, snode.Index ) ){

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl <<  " ["<<snode.Index<<"] "<<  "Update the upper triangular block" << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode.Index<<"] "<<   "blockIdxLocal:" << snode.BlockIdxLocal << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode.Index<<"] "<<   "rowLocalPtr:" << snode.RowLocalPtr << std::endl << std::endl; 
#endif

          std::vector<UBlock>&  Urow = this->U( LBi( snode.Index, grid_ ) );
          std::vector<bool> isBlockFound(Urow.size(),false);

          recvIdx=0;
          for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
            if(isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,grid_),grid_);
              TIMER_START(Recv_L_CrossDiag);

              std::vector<Int> rowLocalPtrRecv;
              std::vector<Int> blockIdxLocalRecv;
              NumMat<Scalar> UUpdateBuf;

              if( MYPROC( grid_ ) != src ){




                std::stringstream sstm;
                Int & sstrSize = arrSstrLcolSizeRecvCD[recvOffset[supidx]+recvIdx];
                std::vector<char> & sstrLcolRecv = arrSstrLcolRecvCD[recvOffset[supidx]+recvIdx];
                sstm.write( &sstrLcolRecv[0], sstrSize );

                deserialize( rowLocalPtrRecv, sstm, NO_MASK );
                deserialize( blockIdxLocalRecv, sstm, NO_MASK );
                deserialize( UUpdateBuf, sstm, NO_MASK );	

                recvIdx++;

              } // sender is not the same as receiver
              else{

                rowLocalPtrRecv   = snode.RowLocalPtr;
                blockIdxLocalRecv = snode.BlockIdxLocal;
                UUpdateBuf = snode.LUpdateBuf;
              } // sender is the same as receiver



              TIMER_STOP(Recv_L_CrossDiag);

#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<" ["<<snode.Index<<"] P"<<MYPROC(grid_)<<" ("<<MYROW(grid_)<<","<<MYCOL(grid_)<<") <--- LBj("<<snode.Index<<") <--- P"<<src<<std::endl;
              statusOFS << std::endl << " ["<<snode.Index<<"] "<<   "rowLocalPtrRecv:" << rowLocalPtrRecv << std::endl << std::endl; 
              statusOFS << std::endl << " ["<<snode.Index<<"] "<<   "blockIdxLocalRecv:" << blockIdxLocalRecv << std::endl << std::endl; 
#endif


              // Update U
              for( Int ib = 0; ib < blockIdxLocalRecv.size(); ib++ ){
                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  UBlock& UB = Urow[jb];
                  if( UB.blockIdx == blockIdxLocalRecv[ib] ){
                    NumMat<Scalar> Ltmp ( UB.numCol, UB.numRow );
                    lapack::Lacpy( 'A', Ltmp.m(), Ltmp.n(), 
                        &UUpdateBuf( rowLocalPtrRecv[ib], 0 ),
                        UUpdateBuf.m(), Ltmp.Data(), Ltmp.m() );
                    isBlockFound[jb] = true;
                    Transpose( Ltmp, UB.nzval );
                    break;
                  }
                }
              }
            }
          }


          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock& UB = Urow[jb];
            if( !isBlockFound[jb] ){
              throw std::logic_error( "UBlock cannot find its update. Something is seriously wrong." );
            }
          }
        } // receiver
      }

      TIMER_STOP(Send_CD_Update_U);

      mpi::Waitall(arrMpiReqsSizeSendCD);
      mpi::Waitall(arrMpiReqsSendCD);
     }; 








inline void PMatrix::UnpackData(SuperNodeBufferType & snode, std::vector<LBlock> & LcolRecv, std::vector<UBlock> & UrowRecv){

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Unpack the received data for processors participate in Gemm. " << std::endl << std::endl; 
#endif
  // U part
//statusOFS<<"Deserializing U"<<std::endl;
  if( MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
    std::stringstream sstm;
    sstm.write( &snode.SstrUrowRecv[0], snode.SstrUrowRecv.size() );
//statusOFS<<"sstm written"<<std::endl;
    std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
    Int numUBlock;
    deserialize( numUBlock, sstm, NO_MASK );
    UrowRecv.resize( numUBlock );
    for( Int jb = 0; jb < numUBlock; jb++ ){
      deserialize( UrowRecv[jb], sstm, mask );
    } 
  } // sender is not the same as receiver
  else{
    // U is obtained locally, just make a copy. Include everything
    // (there is no diagonal block)
    //Is it a copy ?
//    UrowRecv = this->U( LBi( snode.Index, grid_ ) );
    UrowRecv.resize(this->U( LBi( snode.Index, grid_ ) ).size());
    std::copy(this->U( LBi( snode.Index, grid_ ) ).begin(),this->U( LBi( snode.Index, grid_ )).end(),UrowRecv.begin());
  } // sender is the same as receiver


//statusOFS<<"Deserializing L"<<std::endl;
  //L part
  if( MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) ){
    std::stringstream     sstm;
    sstm.write( &snode.SstrLcolRecv[0], snode.SstrLcolRecv.size() );
    std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
    mask[LBlockMask::NZVAL] = 0; // nzval is excluded
    Int numLBlock;
    deserialize( numLBlock, sstm, NO_MASK );
    LcolRecv.resize( numLBlock );
    for( Int ib = 0; ib < numLBlock; ib++ ){
      deserialize( LcolRecv[ib], sstm, mask );
    }
  } // sender is not the same as receiver
  else{
    // L is obtained locally, just make a copy. 
    // Do not include the diagonal block
    std::vector<LBlock>& Lcol =  this->L( LBj( snode.Index, grid_ ) );
    Int startIdx = ( MYROW( grid_ ) == PROW( snode.Index, grid_ ) )?1:0;
      LcolRecv.resize( Lcol.size() - startIdx );
//      for( Int ib = 0; ib < Lcol.size() - startIdx; ib++ ){
//        LcolRecv[ib] = Lcol[ib+startIdx];
//      }
    std::copy(Lcol.begin()+startIdx,Lcol.end(),LcolRecv.begin());

//    if( MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
//      LcolRecv.resize( Lcol.size() );
//      for( Int ib = 0; ib < Lcol.size(); ib++ ){
//        LcolRecv[ib] = Lcol[ib];
//      }
//    }
//    else{
//      LcolRecv.resize( Lcol.size() - 1 );
//      for( Int ib = 0; ib < Lcol.size() - 1; ib++ ){
//        LcolRecv[ib] = Lcol[ib+1];
//      }
//    }
  } // sender is the same as receiver
//statusOFS<<"Done Deserializing"<<std::endl;


////statusOFS<<"LcolRecv"<<std::endl;
////for(Int i =0;i<LcolRecv.size();i++){statusOFS<<LcolRecv[i].blockIdx<<" ";}
////statusOFS<<std::endl;
////
////statusOFS<<"UrowRecv"<<std::endl;
////for(Int i =0;i<UrowRecv.size();i++){statusOFS<<UrowRecv[i].blockIdx<<" ";}
////statusOFS<<std::endl;

};

inline void PMatrix::ComputeDiagUpdate(SuperNodeBufferType & snode){

        //---------Computing  Diagonal block, all processors in the column are participating to all pipelined supernodes
        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){
#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "<<   "Updating the diagonal block" << std::endl << std::endl; 
#endif
          std::vector<LBlock>&  Lcol = this->L( LBj( snode.Index, grid_ ) );

          //Allocate DiagBuf even if Lcol.size() == 0
          snode.DiagBuf.Resize(SuperSize( snode.Index, super_ ), SuperSize( snode.Index, super_ ));
          SetValue(snode.DiagBuf, SCALAR_ZERO);

          // Do I own the diagonal block ?
          Int startIb = (MYROW( grid_ ) == PROW( snode.Index, grid_ ))?1:0;
          for( Int ib = startIb; ib < Lcol.size(); ib++ ){
            blas::Gemm( 'T', 'N', snode.DiagBuf.m(), snode.DiagBuf.n(), Lcol[ib].numRow, 
                SCALAR_MINUS_ONE, &snode.LUpdateBuf( snode.RowLocalPtr[ib-startIb], 0 ), snode.LUpdateBuf.m(),
                Lcol[ib].nzval.Data(), Lcol[ib].nzval.m(), SCALAR_ONE, snode.DiagBuf.Data(), snode.DiagBuf.m() );
          } 

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "<<   "Updated the diagonal block" << std::endl << std::endl; 
#endif
        }
};




  inline void PMatrix::SelInvIntra_P2p(Int lidx)
  {

#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
    Real begin_SendULWaitContentFirst, end_SendULWaitContentFirst, time_SendULWaitContentFirst = 0;
#endif
    Int numSuper = this->NumSuper(); 
    std::vector<std::vector<Int> > & superList = this->WorkingSet();
    Int numSteps = superList.size();
      Int stepSuper = superList[lidx].size(); 




      TIMER_START(AllocateBuffer);

      //This is required to send the size and content of U/L
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendToBelow;
      arrMpireqsSendToBelow.resize( stepSuper, std::vector<MPI_Request>( 2 * grid_->numProcRow, MPI_REQUEST_NULL ));
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendToRight;
      arrMpireqsSendToRight.resize(stepSuper, std::vector<MPI_Request>( 2 * grid_->numProcCol, MPI_REQUEST_NULL ));

      //This is required to reduce L
      std::vector<MPI_Request>  arrMpireqsSendToLeft;
      arrMpireqsSendToLeft.resize(stepSuper, MPI_REQUEST_NULL );
      //This is required to reduce D
      std::vector<MPI_Request>  arrMpireqsSendToAbove;
      arrMpireqsSendToAbove.resize(stepSuper, MPI_REQUEST_NULL );
    
      //This is required to receive the size and content of U/L
      std::vector<MPI_Request>   arrMpireqsRecvSizeFromAny;
      arrMpireqsRecvSizeFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);
      std::vector<MPI_Request>   arrMpireqsRecvContentFromAny;
      arrMpireqsRecvContentFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);


      //allocate the buffers for this supernode
      std::vector<SuperNodeBufferType> arrSuperNodes(stepSuper);
      for (Int supidx=0; supidx<stepSuper; supidx++){ arrSuperNodes[supidx].Index = superList[lidx][supidx];  }
      
      NumMat<Scalar> AinvBuf, UBuf;






      TIMER_STOP(AllocateBuffer);


#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_P2p::UpdateL");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Communication to the Schur complement." << std::endl << std::endl; 
#endif


{
      // Senders
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        std::vector<MPI_Request> & mpireqsSendToBelow = arrMpireqsSendToBelow[supidx];
        std::vector<MPI_Request> & mpireqsSendToRight = arrMpireqsSendToRight[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<  "["<<snode.Index<<"] "<< "Communication for the U part." << std::endl << std::endl; 
#endif


        // Communication for the U part.
        if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) ){
          // Pack the data in U
TIMER_START(Serialize_UL);
          std::stringstream sstm;
          std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
          std::vector<UBlock>&  Urow = this->U( LBi(snode.Index, grid_) );
          // All blocks are to be sent down.
          serialize( (Int)Urow.size(), sstm, NO_MASK );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            serialize( Urow[jb], sstm, mask );
          }
          snode.SstrUrowSend.resize( Size( sstm ) );
          sstm.read( &snode.SstrUrowSend[0], snode.SstrUrowSend.size() );
          snode.SizeSstrUrowSend = snode.SstrUrowSend.size();
TIMER_STOP(Serialize_UL);

          for( Int iProcRow = 0; iProcRow < grid_->numProcRow; iProcRow++ ){
            if( MYROW( grid_ ) != iProcRow &&
                isSendToBelow_( iProcRow,snode.Index ) == true ){
              // Use Isend to send to multiple targets
              MPI_Isend( &snode.SizeSstrUrowSend, 1, MPI_INT,  
                  iProcRow, IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE), grid_->colComm, &mpireqsSendToBelow[2*iProcRow] );
              MPI_Isend( (void*)&snode.SstrUrowSend[0], snode.SizeSstrUrowSend, MPI_BYTE, 
                  iProcRow, IDX_TO_TAG(supidx,SELINV_TAG_U_CONTENT), 
                  grid_->colComm, &mpireqsSendToBelow[2*iProcRow+1] );
#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending U " << snode.SizeSstrUrowSend << " BYTES"<< std::endl <<  std::endl; 
#endif
            } // Send 
          } // for (iProcRow)
        } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "<< "Communication for the L part." << std::endl << std::endl; 
#endif



        // Communication for the L part.
        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){

          TIMER_START(Serialize_UL);
          // Pack the data in L 
          std::stringstream sstm;
          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

          std::vector<LBlock>&  Lcol = this->L( LBj(snode.Index, grid_) );
          // All blocks except for the diagonal block are to be sent right

          if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) )
            serialize( (Int)Lcol.size() - 1, sstm, NO_MASK );
          else
            serialize( (Int)Lcol.size(), sstm, NO_MASK );

          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            if( Lcol[ib].blockIdx > snode.Index ){
#if ( _DEBUGlevel_ >= 2 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Serializing Block index " << Lcol[ib].blockIdx << std::endl;
#endif
              serialize( Lcol[ib], sstm, mask );
            }
          }
          snode.SstrLcolSend.resize( Size( sstm ) );
          sstm.read( &snode.SstrLcolSend[0], snode.SstrLcolSend.size() );
          snode.SizeSstrLcolSend = snode.SstrLcolSend.size();
          TIMER_STOP(Serialize_UL);

          for( Int iProcCol = 0; iProcCol < grid_->numProcCol ; iProcCol++ ){
            if( MYCOL( grid_ ) != iProcCol &&
                isSendToRight_( iProcCol, snode.Index ) == true ){
              // Use Isend to send to multiple targets
              MPI_Isend( &snode.SizeSstrLcolSend, 1, MPI_INT,  
                  iProcCol, IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE), 
                  grid_->rowComm, &mpireqsSendToRight[2*iProcCol] );
              MPI_Isend( (void*)&snode.SstrLcolSend[0], snode.SizeSstrLcolSend, MPI_BYTE, 
                  iProcCol, IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT), 
                  grid_->rowComm, &mpireqsSendToRight[2*iProcCol+1] );
#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending L " << snode.SizeSstrLcolSend<< " BYTES"  << std::endl <<  std::endl; 
#endif
            } // Send 
          } // for (iProcCol)
        } // if I am the sender
      } //Senders

      //TODO Ideally, we should not receive data in sequence but in any order with ksup packed with the data
      // Receivers (Size)
      for (Int supidx=0; supidx<stepSuper ; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        MPI_Request * mpireqsRecvFromAbove = &arrMpireqsRecvSizeFromAny[supidx*2];
        MPI_Request * mpireqsRecvFromLeft = &arrMpireqsRecvSizeFromAny[supidx*2+1];

        // Receive the size first
        if( isRecvFromAbove_( snode.Index ) && 
            MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
          MPI_Irecv( &snode.SizeSstrUrowRecv, 1, MPI_INT, PROW( snode.Index, grid_ ), 
              IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE),
              grid_->colComm, mpireqsRecvFromAbove );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving U size on tag " << IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE)<< std::endl <<  std::endl; 
#endif
        } // if I need to receive from up


        if( isRecvFromLeft_( snode.Index ) &&
            MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) ){
          MPI_Irecv( &snode.SizeSstrLcolRecv, 1, MPI_INT, PCOL( snode.Index, grid_ ), 
              IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE),
              grid_->rowComm, mpireqsRecvFromLeft );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving L size on tag " << IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE)<< std::endl <<  std::endl; 
#endif
        } // if I need to receive from left
      }

      //Wait to receive all the sizes
        TIMER_START(WaitSize_UL);
        mpi::Waitall(arrMpireqsRecvSizeFromAny);
        TIMER_STOP(WaitSize_UL);

      // Receivers (Content)
      for (Int supidx=0; supidx<stepSuper ; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        MPI_Request * mpireqsRecvFromAbove = &arrMpireqsRecvContentFromAny[supidx*2];
        MPI_Request * mpireqsRecvFromLeft = &arrMpireqsRecvContentFromAny[supidx*2+1];

        if( isRecvFromAbove_( snode.Index ) && 
            MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
          snode.SstrUrowRecv.resize( snode.SizeSstrUrowRecv );
          MPI_Irecv( &snode.SstrUrowRecv[0], snode.SizeSstrUrowRecv, MPI_BYTE, 
              PROW( snode.Index, grid_ ), IDX_TO_TAG(supidx,SELINV_TAG_U_CONTENT), 
              grid_->colComm, mpireqsRecvFromAbove );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving U " << snode.SizeSstrUrowRecv << " BYTES"<< std::endl <<  std::endl; 
#endif
        } // if I need to receive from up

        if( isRecvFromLeft_( snode.Index ) &&
            MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) ){
          snode.SstrLcolRecv.resize( snode.SizeSstrLcolRecv );
          MPI_Irecv( &snode.SstrLcolRecv[0], snode.SizeSstrLcolRecv, MPI_BYTE, 
              PCOL( snode.Index, grid_ ), IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT), 
              grid_->rowComm,
              mpireqsRecvFromLeft );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving L " << snode.SizeSstrLcolRecv << " BYTES"<< std::endl <<  std::endl; 
#endif
        } // if I need to receive from left
      }

}







      TIMER_START(Compute_Sinv_LT);
{
      Int gemmProcessed = 0;
      Int gemmToDo = 0;
//      Int toRecvGemm = 0;
      //copy the list of supernodes we need to process
      std::vector<Int> readySupidx;
      //find local things to do
      for(Int supidx = 0;supidx<stepSuper;supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];



//        if( isRecvFromAbove_( snode.Index ) && 
//            MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
//            toRecvGemm++;
//        }
//
//        if( isRecvFromLeft_( snode.Index ) &&
//            MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) ){
//            toRecvGemm++;
//        }

        if( isRecvFromAbove_( snode.Index ) && isRecvFromLeft_( snode.Index )){
          gemmToDo++;
          if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){
            snode.isReady++;
          }

          if(  MYROW( grid_ ) == PROW( snode.Index, grid_ ) ){
            snode.isReady++;
          }

          if(snode.isReady==2){
            readySupidx.push_back(supidx);
#if ( _DEBUGlevel_ >= 1 )
            statusOFS<<std::endl<<"Locally processing ["<<snode.Index<<"]"<<std::endl;
#endif
          }
        }
        else if( (isRecvFromLeft_( snode.Index )  ) && MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) )
        {
          MPI_Request & mpireqsSendToLeft = arrMpireqsSendToLeft[supidx];
          //Dummy 0-b send If I was a receiver, I need to send my data to proc in column of snode.Index
          MPI_Isend( NULL, 0, MPI_BYTE, PCOL(snode.Index,grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_L_REDUCE), grid_->rowComm, &mpireqsSendToLeft );

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYPROC(grid_)<<" has sent "<< 0 << " bytes to " << PNUM(MYROW(grid_),PCOL(snode.Index,grid_),grid_) << std::endl;
#endif
        }// if( isRecvFromAbove_( snode.Index ) && isRecvFromLeft_( snode.Index ))
      }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<std::endl<<"gemmToDo ="<<gemmToDo<<std::endl;
//      statusOFS<<std::endl<<"toRecvGemm ="<<toRecvGemm<<std::endl;
#endif


#if defined (PROFILE) 
      end_SendULWaitContentFirst=0;
      begin_SendULWaitContentFirst=0;
#endif

      while(gemmProcessed<gemmToDo)
      {
        Int reqidx = MPI_UNDEFINED;
        Int supidx = -1;


        //while I don't have anything to do, wait for data to arrive 
        do        {


          int reqIndices[arrMpireqsRecvContentFromAny.size()];
          int numRecv = 0; 

          //then process with the remote ones

          TIMER_START(WaitContent_UL);
#if defined(PROFILE)
          if(begin_SendULWaitContentFirst==0){
            TIMER_START(WaitContent_UL_First);
          }
#endif


//          MPI_Waitany(2*stepSuper, &arrMpireqsRecvContentFromAny[0], &reqidx, MPI_STATUS_IGNORE);
//          TIMER_STOP(WaitContent_UL);
//          {
          numRecv = 0;
          MPI_Waitsome(2*stepSuper, &arrMpireqsRecvContentFromAny[0], &numRecv, reqIndices, MPI_STATUSES_IGNORE);
          TIMER_STOP(WaitContent_UL);

        for(int i =0;i<numRecv;i++){
          reqidx = reqIndices[i];
          //I've received something
          if(reqidx!=MPI_UNDEFINED)
          { 

            supidx = reqidx/2;
            SuperNodeBufferType & snode = arrSuperNodes[supidx];
            snode.isReady++;

#if ( _DEBUGlevel_ >= 1 )
            statusOFS<<std::endl<<"Received data for ["<<snode.Index<<"] reqidx%2="<<reqidx%2<<" is ready ?"<<snode.isReady<<std::endl;
#endif
            //if we received both L and U, the supernode is ready
            if(snode.isReady==2){
              readySupidx.push_back(supidx);

#if defined(PROFILE)
              if(end_SendULWaitContentFirst==0){
                TIMER_STOP(WaitContent_UL_First);
              }
#endif
            }
          }

      }//end for waitsome

        } while(readySupidx.size()==0);






        //If I have some work to do 
        //while(readySupidx.size()>0)
        if(readySupidx.size()>0)
        {
          supidx = readySupidx.back();
          readySupidx.pop_back();
          SuperNodeBufferType & snode = arrSuperNodes[supidx];


          // Only the processors received information participate in the Gemm 
          if( isRecvFromAbove_( snode.Index ) && isRecvFromLeft_( snode.Index ) ){

          std::vector<LBlock> LcolRecv;
          std::vector<UBlock> UrowRecv;
          // Save all the data to be updated for { L( isup, snode.Index ) | isup > snode.Index }.
          // The size will be updated in the Gemm phase and the reduce phase

            UnpackData(snode, LcolRecv, UrowRecv);

            //NumMat<Scalar> AinvBuf, UBuf;
            SelInv_lookup_indexes(snode,LcolRecv, UrowRecv,AinvBuf,UBuf);

            snode.LUpdateBuf.Resize( AinvBuf.m(), SuperSize( snode.Index, super_ ) );
            TIMER_START(Compute_Sinv_LT_GEMM);
            blas::Gemm( 'N', 'T', AinvBuf.m(), UBuf.m(), AinvBuf.n(), SCALAR_MINUS_ONE, 
                AinvBuf.Data(), AinvBuf.m(), 
                UBuf.Data(), UBuf.m(), SCALAR_ZERO,
                snode.LUpdateBuf.Data(), snode.LUpdateBuf.m() ); 
            TIMER_STOP(Compute_Sinv_LT_GEMM);


#if ( _DEBUGlevel_ >= 2 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "snode.LUpdateBuf: " << snode.LUpdateBuf << std::endl;
#endif
          } // if Gemm is to be done locally


          //If I was a receiver, I need to send my data to proc in column of snode.Index
        if( isRecvFromAbove_( snode.Index )  ){
          if( isRecvFromLeft_( snode.Index ) && MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) )
          {
            MPI_Request & mpireqsSendToLeft = arrMpireqsSendToLeft[supidx];

          MPI_Isend( snode.LUpdateBuf.Data(), snode.LUpdateBuf.ByteSize(), MPI_BYTE, PCOL(snode.Index,grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_L_REDUCE), grid_->rowComm, &mpireqsSendToLeft );

#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYCOL(grid_)<<" has sent "<< snode.LUpdateBuf.ByteSize() << " bytes to " << PCOL(snode.Index,grid_) << std::endl;
#endif

          }//Sender
        }
          gemmProcessed++;


#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<std::endl<<"gemmProcessed ="<<gemmProcessed<<"/"<<gemmToDo<<std::endl;
#endif


        }
      }

}
      TIMER_STOP(Compute_Sinv_LT);
      //Reduce Sinv L^T to the processors in PCOL(ksup,grid_)
      TIMER_START(Reduce_Sinv_LT);


      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];


        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){

          //determine the number of rows in LUpdateBufReduced
          Int numRowLUpdateBuf;
          std::vector<LBlock>&  Lcol = this->L( LBj( snode.Index, grid_ ) );
          if( MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
            snode.RowLocalPtr.resize( Lcol.size() + 1 );
            snode.BlockIdxLocal.resize( Lcol.size() );
            snode.RowLocalPtr[0] = 0;
            for( Int ib = 0; ib < Lcol.size(); ib++ ){
              snode.RowLocalPtr[ib+1] = snode.RowLocalPtr[ib] + Lcol[ib].numRow;
              snode.BlockIdxLocal[ib] = Lcol[ib].blockIdx;
            }
          } // I do not own the diagonal block
          else{
            snode.RowLocalPtr.resize( Lcol.size() );
            snode.BlockIdxLocal.resize( Lcol.size() - 1 );
            snode.RowLocalPtr[0] = 0;
            for( Int ib = 1; ib < Lcol.size(); ib++ ){
              snode.RowLocalPtr[ib] = snode.RowLocalPtr[ib-1] + Lcol[ib].numRow;
              snode.BlockIdxLocal[ib-1] = Lcol[ib].blockIdx;
            }
          } // I own the diagonal block, skip the diagonal block
          numRowLUpdateBuf = *snode.RowLocalPtr.rbegin();


          if( numRowLUpdateBuf > 0 ){
            if( snode.LUpdateBuf.m() == 0 && snode.LUpdateBuf.n() == 0 ){
              snode.LUpdateBuf.Resize( numRowLUpdateBuf,SuperSize( snode.Index, super_ ) );
              // Fill zero is important
              SetValue( snode.LUpdateBuf, SCALAR_ZERO );
            }
          }

#if ( _DEBUGlevel_ >= 2 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<<   "LUpdateBuf Before Reduction: " <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif

          Int totCountRecv = 0;
          Int numRecv = CountSendToRight(snode.Index);
          NumMat<Scalar>  LUpdateBufRecv(numRowLUpdateBuf,SuperSize( snode.Index, super_ ) );
          for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
            //Do the blocking recv
            MPI_Status stat;
            Int size = 0;
TIMER_START(L_RECV);
            MPI_Recv(LUpdateBufRecv.Data(), LUpdateBufRecv.ByteSize(), MPI_BYTE, MPI_ANY_SOURCE,IDX_TO_TAG(supidx,SELINV_TAG_L_REDUCE), grid_->rowComm,&stat);
TIMER_STOP(L_RECV);
            MPI_Get_count(&stat, MPI_BYTE, &size);
            //if the processor contributes
            if(size>0){

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYCOL(grid_)<<" has received "<< size << " bytes from " << stat.MPI_SOURCE << std::endl;
#endif
#if ( _DEBUGlevel_ >= 2 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<<   "LUpdateBufRecv: " <<  LUpdateBufRecv << std::endl << std::endl; 
#endif
              //do the sum
              blas::Axpy(snode.LUpdateBuf.Size(), SCALAR_ONE, LUpdateBufRecv.Data(), 1, snode.LUpdateBuf.Data(), 1 );
            }
          } // for (iProcCol)

#if ( _DEBUGlevel_ >= 2 ) 
          statusOFS << std::endl << "["<<snode.Index<<"] "<<   "LUpdateBuf After Reduction: " <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif
        } // Receiver
      }

      TIMER_STOP(Reduce_Sinv_LT);

        mpi::Waitall( arrMpireqsSendToLeft );


      //--------------------- End of reduce of LUpdateBuf-------------------------
#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_P2p::UpdateD");
#endif

      TIMER_START(Update_Diagonal);
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        ComputeDiagUpdate(snode);

        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){
          if( MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
            if(isSendToDiagonal_(snode.Index)){
              //send to above
              MPI_Request & mpireqsSendToAbove = arrMpireqsSendToAbove[supidx];
              MPI_Isend( snode.DiagBuf.Data(),  snode.DiagBuf.ByteSize(), MPI_BYTE,
                  PROW(snode.Index,grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_D_REDUCE), grid_->colComm, &mpireqsSendToAbove );

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYROW(grid_)<<" has sent "<< snode.DiagBuf.ByteSize() << " bytes of DiagBuf to " << PROW(snode.Index,grid_) << " isSendToDiagonal = "<< isSendToDiagonal_(snode.Index) <<  std::endl;
#endif
            }
          }
        }
      }

      TIMER_STOP(Update_Diagonal);



      TIMER_START(Reduce_Diagonal);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){
          if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) ){
            if(snode.DiagBuf.Size()==0){
              snode.DiagBuf.Resize( SuperSize( snode.Index, super_ ), SuperSize( snode.Index, super_ ));
              SetValue(snode.DiagBuf, SCALAR_ZERO);
            }
            //receive from below
            Int totCountRecv = 0;
            Int numRecv = CountRecvFromBelow(snode.Index);
            NumMat<Scalar>  DiagBufRecv(snode.DiagBuf.m(),snode.DiagBuf.n());

            for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
              //Do the blocking recv
              MPI_Status stat;
              Int size = 0;
TIMER_START(D_RECV);
              MPI_Recv(DiagBufRecv.Data(), DiagBufRecv.ByteSize(), MPI_BYTE, MPI_ANY_SOURCE,IDX_TO_TAG(supidx,SELINV_TAG_D_REDUCE), grid_->colComm,&stat);
TIMER_STOP(D_RECV);
              MPI_Get_count(&stat, MPI_BYTE, &size);
              //if the processor contributes
              if(size>0){
                // Add DiagBufRecv to diagonal block.
                blas::Axpy(snode.DiagBuf.Size(), SCALAR_ONE, DiagBufRecv.Data(),
                    1, snode.DiagBuf.Data(), 1 );
              }
            }
            LBlock&  LB = this->L( LBj( snode.Index, grid_ ) )[0];
            // Symmetrize LB
            blas::Axpy( LB.numRow * LB.numCol, SCALAR_ONE, snode.DiagBuf.Data(), 1, LB.nzval.Data(), 1 );
            Symmetrize( LB.nzval );
          }

        } 
      }


      TIMER_STOP(Reduce_Diagonal);

#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_P2p::UpdateU");
#endif

      SendRecvCD_UpdateU(arrSuperNodes, stepSuper);

#ifndef _RELEASE_
      PopCallStack();
#endif

#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_P2p::UpdateLFinal");
#endif

      TIMER_START(Update_L);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Finish updating the L part by filling LUpdateBufReduced back to L" << std::endl << std::endl; 
#endif



        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) && snode.LUpdateBuf.m() > 0 ){
          std::vector<LBlock>&  Lcol = this->L( LBj( snode.Index, grid_ ) );
          //Need to skip the diagonal block if present
          Int startBlock = (MYROW( grid_ ) == PROW( snode.Index, grid_ ))?1:0;
          for( Int ib = startBlock; ib < Lcol.size(); ib++ ){
            LBlock& LB = Lcol[ib];
            lapack::Lacpy( 'A', LB.numRow, LB.numCol, &snode.LUpdateBuf( snode.RowLocalPtr[ib-startBlock], 0 ),
                snode.LUpdateBuf.m(), LB.nzval.Data(), LB.numRow );
          }
        } // Finish updating L	
      } // for (snode.Index) : Main loop


      TIMER_STOP(Update_L);




#ifndef _RELEASE_
      PopCallStack();
#endif


      TIMER_START(Barrier);
        mpi::Waitall(arrMpireqsRecvContentFromAny);
      //Sync for reduce L
//      mpi::Waitall( arrMpireqsSendToLeft );
      //Sync for reduce D
      mpi::Waitall(arrMpireqsSendToAbove);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int ksup = superList[lidx][supidx];
        std::vector<MPI_Request> & mpireqsSendToRight = arrMpireqsSendToRight[supidx];
        std::vector<MPI_Request> & mpireqsSendToBelow = arrMpireqsSendToBelow[supidx];

        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          MPI_Barrier(grid_->colComm);
        }

        mpi::Waitall( mpireqsSendToRight );
        mpi::Waitall( mpireqsSendToBelow );

      }

      TIMER_STOP(Barrier);


      if (options_->maxPipelineDepth!=-1){
        MPI_Barrier(grid_->comm);
      }
  }


  inline void PMatrix::SelInvIntra_Collectives(Int lidx)
  {

#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
    Real begin_SendULWaitContentFirst, end_SendULWaitContentFirst, time_SendULWaitContentFirst = 0;
#endif
    Int numSuper = this->NumSuper(); 
    std::vector<std::vector<Int> > & superList = this->WorkingSet();
    Int numSteps = superList.size();
      Int stepSuper = superList[lidx].size(); 
      TIMER_START(AllocateBuffer);

      //allocate the buffers for this supernode
      std::vector<SuperNodeBufferType> arrSuperNodes(stepSuper);
      for (Int supidx=0; supidx<stepSuper; supidx++){
        arrSuperNodes[supidx].Index = superList[lidx][supidx];
      }

      NumMat<Scalar> AinvBuf, UBuf;
      TIMER_STOP(AllocateBuffer);

#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_Collectives::UpdateL");
#endif


      TIMER_START(WaitContent_UL);
      TIMER_START(WaitContent_UL_Bcast);

        for (Int supidx=0; supidx<stepSuper ; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];


          // Communication for the U part.
          if(countSendToBelow_(snode.Index)>1){
            MPI_Comm * colComm = &commSendToBelow_[supidx];
            Int root = commSendToBelowRoot_[snode.Index];
            bool isRecvFromAbove =  isRecvFromAbove_( snode.Index );
            if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) ){

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl <<"BCAST ["<<snode.Index<<"] Communication to the Schur complement. (Sending)" << std::endl << std::endl; 
#endif
              // Pack the data in U
              std::stringstream sstm;
              std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
              std::vector<UBlock>&  Urow = this->U( LBi(snode.Index, grid_) );
              // All blocks are to be sent down.
              serialize( (Int)Urow.size(), sstm, NO_MASK );
              for( Int jb = 0; jb < Urow.size(); jb++ ){
                serialize( Urow[jb], sstm, mask );
              }
              snode.SstrUrowSend.resize( Size( sstm ) );
              sstm.read( &snode.SstrUrowSend[0], snode.SstrUrowSend.size() );
              snode.SizeSstrUrowSend = snode.SstrUrowSend.size();

              //send size
              MPI_Bcast(&snode.SizeSstrUrowSend, 1 , MPI_INT, root, *colComm );
              //send content
              MPI_Bcast(&snode.SstrUrowSend[0], snode.SizeSstrUrowSend , MPI_BYTE, root, *colComm );
            }
            else if( isRecvFromAbove && 
                MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl <<"BCAST ["<<snode.Index<<"] Communication to the Schur complement. (Receiving)" << std::endl << std::endl; 
#endif
              //size
              MPI_Bcast(&snode.SizeSstrUrowRecv, 1 , MPI_INT, root, *colComm );
              //content
              snode.SstrUrowRecv.resize( snode.SizeSstrUrowRecv );
              MPI_Bcast(&snode.SstrUrowRecv[0], snode.SizeSstrUrowRecv , MPI_BYTE, root, *colComm );
            }
          }
          // Communication for the L part.
          if(countSendToRight_(snode.Index)>1){
            MPI_Comm * rowComm = &commSendToRight_[supidx];
            Int root = commSendToRightRoot_[snode.Index];
            bool isRecvFromLeft =  isRecvFromLeft_( snode.Index );
            if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl <<"BCAST ["<<snode.Index<<"] Communication to the Schur complement. (Sending)" << std::endl << std::endl; 
#endif
              // Pack the data in L 
              std::stringstream sstm;
              std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
              mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

              std::vector<LBlock>&  Lcol = this->L( LBj(snode.Index, grid_) );
              // All blocks except for the diagonal block are to be sent right
              if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) )
                serialize( (Int)Lcol.size() - 1, sstm, NO_MASK );
              else
                serialize( (Int)Lcol.size(), sstm, NO_MASK );

              for( Int ib = 0; ib < Lcol.size(); ib++ ){
                if( Lcol[ib].blockIdx > snode.Index ){
                  serialize( Lcol[ib], sstm, mask );
                }
              }
              snode.SstrLcolSend.resize( Size( sstm ) );
              sstm.read( &snode.SstrLcolSend[0], snode.SstrLcolSend.size() );
              snode.SizeSstrLcolSend = snode.SstrLcolSend.size();
              //send size
              MPI_Bcast(&snode.SizeSstrLcolSend, 1 , MPI_INT, root, *rowComm );
              //send content
              MPI_Bcast(&snode.SstrLcolSend[0], snode.SizeSstrLcolSend , MPI_BYTE, root, *rowComm );
            }
            else if( isRecvFromLeft && MYCOL( grid_ ) != PCOL( snode.Index, grid_ ) ){
#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl <<"BCAST ["<<snode.Index<<"] Communication to the Schur complement. (Receiving)" << std::endl << std::endl; 
#endif
              //size
              MPI_Bcast(&snode.SizeSstrLcolRecv, 1 , MPI_INT, root, *rowComm );
              //content
              snode.SstrLcolRecv.resize( snode.SizeSstrLcolRecv );
              MPI_Bcast(&snode.SstrLcolRecv[0], snode.SizeSstrLcolRecv , MPI_BYTE, root, *rowComm );
            }
          }

        }

      TIMER_STOP(WaitContent_UL_Bcast);
      TIMER_STOP(WaitContent_UL);

      TIMER_START(Compute_Sinv_LT);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];
        // Overlap the communication with computation.  All processors move
        // to Gemm phase when ready 
        // Only the processors received information participate in the Gemm 
        if( isRecvFromAbove_( snode.Index ) && isRecvFromLeft_( snode.Index ) ){

        std::vector<LBlock> LcolRecv;
        std::vector<UBlock> UrowRecv;
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "BCAST ["<<snode.Index<<"] "<<  "Main work: Gemm" << std::endl << std::endl; 
#endif

        // Save all the data to be updated for { L( isup, snode.Index ) | isup > snode.Index }.
        // The size will be updated in the Gemm phase and the reduce phase
          UnpackData(snode, LcolRecv, UrowRecv);

          SelInv_lookup_indexes(snode,LcolRecv, UrowRecv,AinvBuf,UBuf);

          snode.LUpdateBuf.Resize( AinvBuf.m(), SuperSize( snode.Index, super_ ) );

          TIMER_START(Compute_Sinv_LT_GEMM);


#if ( _DEBUGlevel_ >= 2 )
          statusOFS << std::endl << "BCAST ["<<snode.Index<<"] "<<  "AinvBuf: " << AinvBuf << std::endl;
          statusOFS << std::endl << "BCAST ["<<snode.Index<<"] "<<  "UBuf: " << UBuf << std::endl;
#endif
          // Gemm for snode.LUpdateBuf = -AinvBuf * UBuf^T
          blas::Gemm( 'N', 'T', AinvBuf.m(), UBuf.m(), AinvBuf.n(), SCALAR_MINUS_ONE, 
              AinvBuf.Data(), AinvBuf.m(), 
              UBuf.Data(), UBuf.m(), SCALAR_ZERO,
              snode.LUpdateBuf.Data(), snode.LUpdateBuf.m() ); 

          TIMER_STOP(Compute_Sinv_LT_GEMM);

#if ( _DEBUGlevel_ >= 2 )
          statusOFS << std::endl << "BCAST ["<<snode.Index<<"] "<<  "snode.LUpdateBuf: " << snode.LUpdateBuf << std::endl;
#endif
        } // if Gemm is to be done locally

      }

      TIMER_STOP(Compute_Sinv_LT);


      //Reduce Sinv L^T to the processors in PCOL(snode.Index,grid_)
      TIMER_START(Reduce_Sinv_LT);
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        // Processor column of snode.Index collects the symbolic data for snode.LUpdateBuf.
        Int numRowLUpdateBuf;


        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){
          std::vector<LBlock>&  Lcol = this->L( LBj( snode.Index, grid_ ) );
          if( MYROW( grid_ ) != PROW( snode.Index, grid_ ) ){
            snode.RowLocalPtr.resize( Lcol.size() + 1 );
            snode.BlockIdxLocal.resize( Lcol.size() );
            snode.RowLocalPtr[0] = 0;
            for( Int ib = 0; ib < Lcol.size(); ib++ ){
              snode.RowLocalPtr[ib+1] = snode.RowLocalPtr[ib] + Lcol[ib].numRow;
              snode.BlockIdxLocal[ib] = Lcol[ib].blockIdx;
            }
          } // I do not own the diagonal block
          else{
            snode.RowLocalPtr.resize( Lcol.size() );
            snode.BlockIdxLocal.resize( Lcol.size() - 1 );
            snode.RowLocalPtr[0] = 0;
            for( Int ib = 1; ib < Lcol.size(); ib++ ){
              snode.RowLocalPtr[ib] = snode.RowLocalPtr[ib-1] + Lcol[ib].numRow;
              snode.BlockIdxLocal[ib-1] = Lcol[ib].blockIdx;
            }
          } // I owns the diagonal block, skip the diagonal block
          numRowLUpdateBuf = *snode.RowLocalPtr.rbegin();
          if( numRowLUpdateBuf > 0 ){
            if( snode.LUpdateBuf.m() == 0 && snode.LUpdateBuf.n() == 0 ){
              snode.LUpdateBuf.Resize( numRowLUpdateBuf, SuperSize( snode.Index, super_ ) );
              // Fill zero is important
              SetValue( snode.LUpdateBuf, SCALAR_ZERO );
            }
          }
        } 


        if(countSendToRight_(snode.Index)>1){

          MPI_Comm * rowComm = &commSendToRight_[supidx];
          Int root = commSendToRightRoot_[snode.Index];

          // Processor column sends the total row dimension to all processors
          // in the same row to prepare for reduce
          MPI_Bcast( &numRowLUpdateBuf, 1, MPI_INT, root, *rowComm );

          // If LUpdatebuf has not been constructed, resize and fill with zero
          if( numRowLUpdateBuf > 0 ){
            if( snode.LUpdateBuf.m() == 0 && snode.LUpdateBuf.n() == 0 ){
              snode.LUpdateBuf.Resize( numRowLUpdateBuf, SuperSize( snode.Index, super_ ) );
              // Fill zero is important
              SetValue( snode.LUpdateBuf, SCALAR_ZERO );
            }

        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){
            mpi::Reduce( (Scalar*)MPI_IN_PLACE, snode.LUpdateBuf.Data(),
                numRowLUpdateBuf * SuperSize( snode.Index, super_ ), MPI_SUM, 
                root, *rowComm );
        }
        else{
            mpi::Reduce( snode.LUpdateBuf.Data(), NULL,
                numRowLUpdateBuf * SuperSize( snode.Index, super_ ), MPI_SUM, 
                root, *rowComm );
        }

            if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){

#if ( _DEBUGlevel_ >= 2 ) 
              statusOFS << std::endl << "BCAST ["<<snode.Index<<"] "<<   "LUpdateBuf: " <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif

            }

          } // Perform reduce for nonzero block rows in the column of snode.Index


        }


      }

      TIMER_STOP(Reduce_Sinv_LT);
      //--------------------- End of reduce of snode.LUpdateBuf-------------------------

#ifndef _RELEASE_
      PopCallStack();
#endif

#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_Collectives::UpdateD");
#endif

      TIMER_START(Update_Diagonal);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        ComputeDiagUpdate(snode);

      } //end for
      TIMER_STOP(Update_Diagonal);


      TIMER_START(Reduce_Diagonal);


      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) ){


#if ( _DEBUGlevel_ >= 2 )
          statusOFS << std::endl << "BCAST ["<<snode.Index<<"] "<<   "snode.DiagBuf: " << snode.DiagBuf << std::endl << std::endl; 
#endif



          if(countRecvFromBelow_(snode.Index)>1)
          {
            MPI_Comm * colComm = &commRecvFromBelow_[supidx];
            Int root = commRecvFromBelowRoot_[snode.Index];

            if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) ){
              if(snode.DiagBuf.m()==0 && snode.DiagBuf.n()==0){
                snode.DiagBuf.Resize( SuperSize( snode.Index, super_ ), SuperSize( snode.Index, super_ ));
                SetValue(snode.DiagBuf, SCALAR_ZERO);
              }

              mpi::Reduce((Scalar*) MPI_IN_PLACE, snode.DiagBuf.Data(), 
                  SuperSize( snode.Index, super_ ) * SuperSize( snode.Index, super_ ),
                  MPI_SUM, root, *colComm );
            }
            else{

              mpi::Reduce( snode.DiagBuf.Data(), NULL, 
                  SuperSize( snode.Index, super_ ) * SuperSize( snode.Index, super_ ),
                  MPI_SUM, root, *colComm );
            }

          }

          // Add DiagBufReduced to diagonal block.
          if( MYROW( grid_ ) == PROW( snode.Index, grid_ ) ){

            if(snode.DiagBuf.m()>0 && snode.DiagBuf.n()>0){
              LBlock&  LB = this->L( LBj( snode.Index, grid_ ) )[0];
              // Symmetrize LB
              blas::Axpy( LB.numRow * LB.numCol, SCALAR_ONE, snode.DiagBuf.Data(),
                  1, LB.nzval.Data(), 1 );
              Symmetrize( LB.nzval );
            }
          }
        } 
      }

      TIMER_STOP(Reduce_Diagonal);

#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_Collectives::UpdateU");
#endif

      SendRecvCD_UpdateU(arrSuperNodes, stepSuper);

#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_Collectives::UpdateLFinal");
#endif

      TIMER_START(Update_L);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferType & snode = arrSuperNodes[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Finish updating the L part by filling LUpdateBufReduced back to L" << std::endl << std::endl; 
#endif
        if( MYCOL( grid_ ) == PCOL( snode.Index, grid_ ) && snode.LUpdateBuf.m() > 0 ){
          std::vector<LBlock>&  Lcol = this->L( LBj( snode.Index, grid_ ) );
          //Need to skip the diagonal block if present
          Int startBlock = (MYROW( grid_ ) == PROW( snode.Index, grid_ ))?1:0;
          for( Int ib = startBlock; ib < Lcol.size(); ib++ ){
            LBlock& LB = Lcol[ib];
            lapack::Lacpy( 'A', LB.numRow, LB.numCol, &snode.LUpdateBuf( snode.RowLocalPtr[ib-startBlock], 0 ),
                snode.LUpdateBuf.m(), LB.nzval.Data(), LB.numRow );
          }
        } // Finish updating L	
      } // for (snode.Index) : Main loop


      TIMER_STOP(Update_L);

      if (options_->maxPipelineDepth!=-1){
        MPI_Barrier(grid_->comm);
      }

#ifndef _RELEASE_
      PopCallStack();
#endif
  }

        void PMatrix::ConstructCommunicators_Collectives(Int lidx){

          Int numSuper = this->NumSuper(); 
          std::vector<std::vector<Int> > & superList = this->WorkingSet();
          Int numSteps = superList.size();
          Int stepSuper = superList[lidx].size();

      MPI_Group rowCommGroup;  
      MPI_Comm_group(grid_->rowComm,&rowCommGroup);
      MPI_Group colCommGroup; 
      MPI_Comm_group(grid_->colComm,&colCommGroup);


        //create the communicators
        commSendToRight_.resize(stepSuper,MPI_COMM_NULL);
        commSendToBelow_.resize(stepSuper,MPI_COMM_NULL);
        commRecvFromBelow_.resize(stepSuper,MPI_COMM_NULL);
        for(Int supidx = 0;supidx<stepSuper;supidx++){
          Int ksup = superList[lidx][supidx];

          Int count= std::count(commSendToRightMask_[ksup].begin(), commSendToRightMask_[ksup].end(), true);


          if(count>1){
            MPI_Group commGroup;  

            Int color = commSendToRightMask_[ksup][MYCOL(grid_)];

            MPI_Comm_split(grid_->rowComm, color  ,MYCOL(grid_) , &commSendToRight_[supidx]);
            MPI_Comm_group(commSendToRight_[supidx],&commGroup);
            //now for each supernode, we need to store the pointer to the communnicator and the rank of the root
            Int curRoot = PCOL(ksup,grid_);
            Int newRank = -1;
            if(color>0){
              MPI_Group_translate_ranks(rowCommGroup, 1,&curRoot,commGroup, &newRank);

              if(newRank==MPI_UNDEFINED){
                statusOFS<<"["<<ksup<<"] Root COL "<<curRoot<<" has no matching rank"<<std::endl;
              }
            }

            commSendToRightRoot_[ksup] = newRank;
          }

          count= std::count(commSendToBelowMask_[ksup].begin(), commSendToBelowMask_[ksup].end(), true);

          if(count>1){
            MPI_Group commGroup;  

            Int color = commSendToBelowMask_[ksup][MYROW(grid_)];

            MPI_Comm_split(grid_->colComm, color  ,MYROW(grid_) , &commSendToBelow_[supidx]);
            MPI_Comm_group(commSendToBelow_[supidx],&commGroup);
            //now for each supernode, we need to store the pointer to the communnicator and the rank of the root
            Int curRoot = PROW(ksup,grid_);
            Int newRank = -1;
            if(color>0){
              MPI_Group_translate_ranks(colCommGroup, 1,&curRoot,commGroup, &newRank);

              if(newRank==MPI_UNDEFINED){
                statusOFS<<"["<<ksup<<"] Root ROW "<<curRoot<<" has no matching rank"<<std::endl;
              }
            }

            commSendToBelowRoot_[ksup] = newRank;
          }

          count= std::count(commRecvFromBelowMask_[ksup].begin(), commRecvFromBelowMask_[ksup].end(), true);
          if(count>1){
            MPI_Group commGroup;  

            Int color = commRecvFromBelowMask_[ksup][MYROW(grid_)];

            MPI_Comm_split(grid_->colComm, color  ,MYROW(grid_) , &commRecvFromBelow_[supidx]);
            MPI_Comm_group(commRecvFromBelow_[supidx],&commGroup);
            //now for each supernode, we need to store the pointer to the communnicator and the rank of the root
            Int curRoot = PROW(ksup,grid_);
            Int newRank = -1;
            if(color>0){
              MPI_Group_translate_ranks(colCommGroup, 1,&curRoot,commGroup, &newRank);

              if(newRank==MPI_UNDEFINED){
                statusOFS<<"["<<ksup<<"] Root ROW "<<curRoot<<" has no matching rank"<<std::endl;
              }
            }

            commRecvFromBelowRoot_[ksup] = newRank;
          }

        }
        };

  void PMatrix::SelInv_Hybrid	(Int threshold = BCAST_THRESHOLD)
  {

#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
//    TAU_PROFILE_SET_CONTEXT(grid_->comm);
#endif


    TIMER_START(SelInv_Hybrid);

#ifndef _RELEASE_
    PushCallStack("PMatrix::SelInv_Hybrid");
#endif



    Int numSuper = this->NumSuper(); 

      commSendToRightRoot_.resize(numSuper);
      commSendToBelowRoot_.resize(numSuper);
      commRecvFromBelowRoot_.resize(numSuper);

    getMaxCommunicatorSizes();

    // Main loop
    std::vector<std::vector<Int> > & superList = this->WorkingSet();
    Int numSteps = superList.size();

    for (Int lidx=0; lidx<numSteps ; lidx++){
      Int stepSuper = superList[lidx].size(); 

      Int logMaxSize = std::ceil(log2(maxCommSizes_[lidx]));
//      statusOFS<<maxCommSizes_[lidx]<<" vs "<<logMaxSize*BCAST_THRESHOLD<<std::endl;

      if(maxCommSizes_[lidx]>1 &&  maxCommSizes_[lidx]  > logMaxSize*threshold){ 
//        statusOFS<<"BCAST VARIANT"<<std::endl;
        ConstructCommunicators_Collectives(lidx);
        SelInvIntra_Collectives(lidx);
        DestructCommunicators_Collectives	(  );
      }
      else
      {
//        statusOFS<<"PIPELINE VARIANT"<<std::endl;
        SelInvIntra_P2p(lidx);
      }
    }

#ifndef _RELEASE_
    PopCallStack();
#endif

    TIMER_STOP(SelInv_Hybrid);


    return ;
  } 		// -----  end of method PMatrix::SelInv_Hybrid  ----- 



  void PMatrix::SelInv_Collectives	(  )
  {

#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
//    TAU_PROFILE_SET_CONTEXT(grid_->comm);
#endif


    TIMER_START(SelInv_Collectives);

#ifndef _RELEASE_
    PushCallStack("PMatrix::SelInv_Collectives");
#endif



    Int numSuper = this->NumSuper(); 

      commSendToRightRoot_.resize(numSuper);
      commSendToBelowRoot_.resize(numSuper);
      commRecvFromBelowRoot_.resize(numSuper);


    // Main loop
    std::vector<std::vector<Int> > & superList = this->WorkingSet();
    Int numSteps = superList.size();

    for (Int lidx=0; lidx<numSteps ; lidx++){
      Int stepSuper = superList[lidx].size(); 
        ConstructCommunicators_Collectives(lidx);
        SelInvIntra_Collectives(lidx);
        DestructCommunicators_Collectives	(  );
    }

#ifndef _RELEASE_
    PopCallStack();
#endif

    TIMER_STOP(SelInv_Collectives);


    return ;
  } 		// -----  end of method PMatrix::SelInv_Collectives  ----- 

  void PMatrix::ConstructCommunicationPattern	(  )
  {
#ifndef _RELEASE_
    PushCallStack("PMatrix::ConstructCommunicationPattern");
#endif
    Int numSuper = this->NumSuper();
#ifndef _RELEASE_
    PushCallStack( "Initialize the communication pattern" );
#endif
    isSendToBelow_.Resize(grid_->numProcRow, numSuper);
    isSendToRight_.Resize(grid_->numProcCol, numSuper);
    SetValue( isSendToBelow_, false );
    SetValue( isSendToRight_, false );

    isSendToCrossDiagonal_.Resize(grid_->numProcCol+1, numSuper );
    SetValue( isSendToCrossDiagonal_, false );
    isRecvFromCrossDiagonal_.Resize(grid_->numProcRow+1, numSuper );
    SetValue( isRecvFromCrossDiagonal_, false );

    isRecvFromAbove_.Resize( numSuper );
    isRecvFromLeft_.Resize( numSuper );
    SetValue( isRecvFromAbove_, false );
    SetValue( isRecvFromLeft_, false );
#ifndef _RELEASE_
    PopCallStack();
#endif

      std::vector<Int> snodeEtree(this->NumSuper());
      GetEtree(snodeEtree);

#ifndef _RELEASE_
    PushCallStack( "Local column communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local column communication" << std::endl;
#endif
    // localColBlockRowIdx stores the nonzero block indices for each local block column.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockCol x numNonzeroBlock
    std::vector<std::set<Int> >   localColBlockRowIdx;

    localColBlockRowIdx.resize( this->NumLocalBlockCol() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        std::vector<Int>  tBlockRowIdx;
        tBlockRowIdx.clear();

        // L part
        std::vector<LBlock>& Lcol = this->L( LBj(ksup, grid_) );
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          tBlockRowIdx.push_back( Lcol[ib].blockIdx );
        }

        // U part
        for( Int ib = 0; ib < this->NumLocalBlockRow(); ib++ ){
          std::vector<UBlock>& Urow = this->U(ib);
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            if( Urow[jb].blockIdx == ksup ){
              tBlockRowIdx.push_back( GBi( ib, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockRowIdx;
        mpi::Allgatherv( tBlockRowIdx, tAllBlockRowIdx, grid_->colComm );

        localColBlockRowIdx[LBj( ksup, grid_ )].insert(
            tAllBlockRowIdx.begin(), tAllBlockRowIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Column block " << ksup 
          << " has the following nonzero block rows" << std::endl;
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for(ksup)


#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack( "Local row communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "Local row communication" << std::endl;
#endif
    // localRowBlockColIdx stores the nonzero block indices for each local block row.
    // The nonzero block indices including contribution from both L and U.
    // Dimension: numLocalBlockRow x numNonzeroBlock
    std::vector<std::set<Int> >   localRowBlockColIdx;

    localRowBlockColIdx.resize( this->NumLocalBlockRow() );

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      // All block columns perform independently
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        std::vector<Int>  tBlockColIdx;
        tBlockColIdx.clear();

        // U part
        std::vector<UBlock>& Urow = this->U( LBi(ksup, grid_) );
        for( Int jb = 0; jb < Urow.size(); jb++ ){
          tBlockColIdx.push_back( Urow[jb].blockIdx );
        }

        // L part
        for( Int jb = 0; jb < this->NumLocalBlockCol(); jb++ ){
          std::vector<LBlock>& Lcol = this->L(jb);
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            if( Lcol[ib].blockIdx == ksup ){
              tBlockColIdx.push_back( GBj( jb, grid_ ) );
            }
          }
        }

        // Communication
        std::vector<Int> tAllBlockColIdx;
        mpi::Allgatherv( tBlockColIdx, tAllBlockColIdx, grid_->rowComm );

        localRowBlockColIdx[LBi( ksup, grid_ )].insert(
            tAllBlockColIdx.begin(), tAllBlockColIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS 
          << " Row block " << ksup 
          << " has the following nonzero block columns" << std::endl;
        for( std::set<Int>::iterator si = localRowBlockColIdx[LBi( ksup, grid_ )].begin();
            si != localRowBlockColIdx[LBi( ksup, grid_ )].end();
            si++ ){
          statusOFS << *si << "  ";
        }
        statusOFS << std::endl; 
#endif

      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
    } // for(ksup)

#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack("SendToBelow / RecvFromAbove");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes to the right of ksup


      Int jsup = snodeEtree[ksup];
      while(jsup<numSuper){
        Int jsupLocalBlockCol = LBj( jsup, grid_ );
        Int jsupProcCol = PCOL( jsup, grid_ );
        if( MYCOL( grid_ ) == jsupProcCol ){

          // SendToBelow / RecvFromAbove only if (ksup, jsup) is nonzero.
          if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ) {
            for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
                si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
              Int isup = *si;
              Int isupProcRow = PROW( isup, grid_ );
              if( isup > ksup ){
                if( MYROW( grid_ ) == isupProcRow ){
                  isRecvFromAbove_(ksup) = true;
                }
                if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
                  isSendToBelow_( isupProcRow, ksup ) = true;
                }
              } // if( isup > ksup )
            } // for (si)
          } // if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 )

        } // if( MYCOL( grid_ ) == PCOL( jsup, grid_ ) )

        jsup = snodeEtree[jsup];

      } // for(jsup)
    } // for(ksup)

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToBelow:" << std::endl;
    for(int j = 0;j< isSendToBelow_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToBelow_.m();i++){
        statusOFS<< isSendToBelow_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromAbove:" << std::endl;
    for(int j = 0;j< isRecvFromAbove_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromAbove_(j)<<std::endl;
    }
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif

#ifndef _RELEASE_
    PushCallStack("SendToRight / RecvFromLeft");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      // Loop over all the supernodes below ksup

      Int isup = snodeEtree[ksup];
      while(isup<numSuper){
        Int isupLocalBlockRow = LBi( isup, grid_ );
        Int isupProcRow       = PROW( isup, grid_ );
        if( MYROW( grid_ ) == isupProcRow ){
          // SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.
          if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
            for( std::set<Int>::iterator si = localRowBlockColIdx[isupLocalBlockRow].begin();
                si != localRowBlockColIdx[isupLocalBlockRow].end(); si++ ){
              Int jsup = *si;
              Int jsupProcCol = PCOL( jsup, grid_ );
              if( jsup > ksup ){

                if( MYCOL( grid_ ) == jsupProcCol ){
                  isRecvFromLeft_(ksup) = true;
                }
                if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
                  isSendToRight_( jsupProcCol, ksup ) = true;
                }
              }
            } // for (si)
          } // if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 )
        } // if( MYROW( grid_ ) == isupProcRow )
        isup = snodeEtree[isup];

      } // for (isup)
    }	 // for (ksup)

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToRight:" << std::endl;
    for(int j = 0;j< isSendToRight_.n();j++){
      statusOFS << "["<<j<<"] ";
      for(int i =0; i < isSendToRight_.m();i++){
        statusOFS<< isSendToRight_(i,j) << " ";
      }
      statusOFS<<std::endl;
    }

    statusOFS << std::endl << "isRecvFromLeft:" << std::endl;
    for(int j = 0;j< isRecvFromLeft_.m();j++){
      statusOFS << "["<<j<<"] "<< isRecvFromLeft_(j)<<std::endl;
    }
#endif

#ifndef _RELEASE_
    PopCallStack();
#endif

#ifndef _RELEASE_
    PushCallStack("SendToCrossDiagonal / RecvFromCrossDiagonal");
#endif
    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
            si != localColBlockRowIdx[LBj( ksup, grid_ )].end(); si++ ){
          Int isup = *si;
          Int isupProcRow = PROW( isup, grid_ );
          Int isupProcCol = PCOL( isup, grid_ );
          if( isup > ksup && MYROW( grid_ ) == isupProcRow ){
            isSendToCrossDiagonal_(grid_->numProcCol, ksup ) = true;
            isSendToCrossDiagonal_(isupProcCol, ksup ) = true;
          }
        } // for (si)
      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for (ksup)

    for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
        for( std::set<Int>::iterator si = localRowBlockColIdx[ LBi(ksup, grid_) ].begin();
            si != localRowBlockColIdx[ LBi(ksup, grid_) ].end(); si++ ){
          Int jsup = *si;
          Int jsupProcCol = PCOL( jsup, grid_ );
          Int jsupProcRow = PROW( jsup, grid_ );
          if( jsup > ksup && MYCOL(grid_) == jsupProcCol ){
            isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) = true;
            isRecvFromCrossDiagonal_(jsupProcRow, ksup ) = true;
          }
        } // for (si)
      } // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
    } // for (ksup)
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "isSendToCrossDiagonal:" << std::endl;
    for(int j =0; j < isSendToCrossDiagonal_.n();j++){
      if(isSendToCrossDiagonal_(grid_->numProcCol,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isSendToCrossDiagonal_.m()-1;i++){
          if(isSendToCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(PROW(j,grid_),i,grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }

    statusOFS << std::endl << "isRecvFromCrossDiagonal:" << std::endl;
    for(int j =0; j < isRecvFromCrossDiagonal_.n();j++){
      if(isRecvFromCrossDiagonal_(grid_->numProcRow,j)){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < isRecvFromCrossDiagonal_.m()-1;i++){
          if(isRecvFromCrossDiagonal_(i,j))
          {
            statusOFS<< PNUM(i,PCOL(j,grid_),grid_)<<" ";
          }
        }
        statusOFS<<std::endl;
      }
    }


#endif

#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PopCallStack();
#endif


    return ;
  } 		// -----  end of method PMatrix::ConstructCommunicationPattern  ----- 



  void PMatrix::SelInv	(  )
  {
#ifndef _RELEASE_
    PushCallStack("PMatrix::SelInv");
#endif
    Int numSuper = this->NumSuper(); 


#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
//    TAU_PROFILE_SET_CONTEXT(grid_->comm);
#endif

    TIMER_START(SelInv);


    for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv::UpdateL");
#endif


      // Communication for the U part.
      std::vector<MPI_Request> mpireqsSendToBelow( 2 * grid_->numProcRow, MPI_REQUEST_NULL );
      std::vector<char> sstrUrowSend;

      // Senders
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<"ORIGINAL ["<<ksup<<"] Communication to the Schur complement. (Sending)" << std::endl << std::endl; 
#endif
        // Pack the data in U
        std::stringstream sstm;
        std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
        std::vector<UBlock>&  Urow = this->U( LBi(ksup, grid_) );
        // All blocks are to be sent down.
        serialize( (Int)Urow.size(), sstm, NO_MASK );
        for( Int jb = 0; jb < Urow.size(); jb++ ){
          serialize( Urow[jb], sstm, mask );
        }
        sstrUrowSend.resize( Size( sstm ) );
        sstm.read( &sstrUrowSend[0], sstrUrowSend.size() );
        for( Int iProcRow = 0; iProcRow < grid_->numProcRow; iProcRow++ ){
          if( MYROW( grid_ ) != iProcRow &&
              isSendToBelow_( iProcRow,ksup ) == true ){
            // Use Isend to send to multiple targets
            Int sizeStm = sstrUrowSend.size();
            MPI_Isend( &sizeStm, 1, MPI_INT,  
                iProcRow, SELINV_TAG_U_SIZE, 
                grid_->colComm, &mpireqsSendToBelow[grid_->numProcRow + iProcRow] );
            MPI_Isend( (void*)&sstrUrowSend[0], sizeStm, MPI_BYTE, 
                iProcRow, SELINV_TAG_U_CONTENT, 
                grid_->colComm, &mpireqsSendToBelow[iProcRow] );
          } // Send 
        } // for (iProcRow)
      } // if I am the sender

      // Communication for the L part.

      std::vector<MPI_Request> mpireqsSendToRight( 2 * grid_->numProcCol, MPI_REQUEST_NULL ); 
      std::vector<char> sstrLcolSend;

      // Senders
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        // Pack the data in L 
        std::stringstream sstm;
        std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
        mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

        std::vector<LBlock>&  Lcol = this->L( LBj(ksup, grid_) );
        // All blocks except for the diagonal block are to be sent right
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
          serialize( (Int)Lcol.size() - 1, sstm, NO_MASK );
        else
          serialize( (Int)Lcol.size(), sstm, NO_MASK );

        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          if( Lcol[ib].blockIdx > ksup ){
            serialize( Lcol[ib], sstm, mask );
          }
        }
        sstrLcolSend.resize( Size( sstm ) );
        sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
        for( Int iProcCol = 0; iProcCol < grid_->numProcCol ; iProcCol++ ){
          if( MYCOL( grid_ ) != iProcCol &&
              isSendToRight_( iProcCol, ksup ) == true ){
            // Use Isend to send to multiple targets
            Int sizeStm = sstrLcolSend.size();
            MPI_Isend( &sizeStm, 1, MPI_INT,  
                iProcCol, SELINV_TAG_L_SIZE, 
                grid_->rowComm, &mpireqsSendToRight[grid_->numProcCol + iProcCol] );
            MPI_Isend( (void*)&sstrLcolSend[0], sizeStm, MPI_BYTE, 
                iProcCol, SELINV_TAG_L_CONTENT, 
                grid_->rowComm, &mpireqsSendToRight[iProcCol] );
          } // Send 
        } // for (iProcCol)
      } // if I am the sender

      // Receive
      std::vector<MPI_Request> mpireqsRecvFromAbove( 2, MPI_REQUEST_NULL ); 
      std::vector<MPI_Request> mpireqsRecvFromLeft( 2, MPI_REQUEST_NULL ); 
      Int sizeStmFromAbove, sizeStmFromLeft;

      std::vector<char>     sstrLcolRecv;
      std::vector<char>     sstrUrowRecv;

      // Receive the size first

      if( isRecvFromAbove_( ksup ) && 
          MYROW( grid_ ) != PROW( ksup, grid_ ) ){
        MPI_Irecv( &sizeStmFromAbove, 1, MPI_INT, PROW( ksup, grid_ ), 
            SELINV_TAG_U_SIZE,
            grid_->colComm, &mpireqsRecvFromAbove[0] );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<"ORIGINAL ["<<ksup<<"] Communication to the Schur complement. (Receiving)" << std::endl << std::endl; 
#endif
      } // if I need to receive from up


      if( isRecvFromLeft_( ksup ) &&
          MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){
        MPI_Irecv( &sizeStmFromLeft, 1, MPI_INT, PCOL( ksup, grid_ ), 
            SELINV_TAG_L_SIZE,
            grid_->rowComm, &mpireqsRecvFromLeft[0] );

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<"ORIGINAL ["<<ksup<<"] Communication to the Schur complement. (Receiving)" << std::endl << std::endl; 
#endif
      } // if I need to receive from left





      TIMER_START(WaitSize_UL);
      // Wait to obtain size information
      mpi::Wait( mpireqsRecvFromAbove[0] );
      mpi::Wait( mpireqsRecvFromLeft[0] );

      TIMER_STOP(WaitSize_UL);

      if( isRecvFromAbove_( ksup ) && 
          MYROW( grid_ ) != PROW( ksup, grid_ ) ){

        sstrUrowRecv.resize( sizeStmFromAbove );
        MPI_Irecv( &sstrUrowRecv[0], sizeStmFromAbove, MPI_BYTE, 
            PROW( ksup, grid_ ), SELINV_TAG_U_CONTENT, 
            grid_->colComm, &mpireqsRecvFromAbove[1] );
      } // if I need to receive from up

      if( isRecvFromLeft_( ksup ) &&
          MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){

        sstrLcolRecv.resize( sizeStmFromLeft );
        MPI_Irecv( &sstrLcolRecv[0], sizeStmFromLeft, MPI_BYTE, 
            PCOL( ksup, grid_ ), SELINV_TAG_L_CONTENT, 
            grid_->rowComm,
            &mpireqsRecvFromLeft[1] );
      } // if I need to receive from left

      // Wait for all communication to finish
      // Wait to obtain packed information in a string and then write into stringstream
      TIMER_START(WaitContent_UL);

      mpi::Wait( mpireqsRecvFromAbove[1] );
      mpi::Wait( mpireqsRecvFromLeft[1] );
      TIMER_STOP(WaitContent_UL);




      // Overlap the communication with computation.  All processors move
      // to Gemm phase when ready 
      TIMER_START(Compute_Sinv_LT);



      std::vector<LBlock>   LcolRecv;
      std::vector<UBlock>   UrowRecv;


      if( isRecvFromAbove_( ksup ) && isRecvFromLeft_( ksup ) ){

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "Unpack the received data for processors participate in Gemm. " << std::endl << std::endl; 
#endif

        // U part
        if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
          std::stringstream sstm;
          sstm.write( &sstrUrowRecv[0], sizeStmFromAbove );
          std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
          Int numUBlock;
          deserialize( numUBlock, sstm, NO_MASK );
          UrowRecv.resize( numUBlock );
          for( Int jb = 0; jb < numUBlock; jb++ ){
            deserialize( UrowRecv[jb], sstm, mask );
          } 
        } // sender is not the same as receiver
        else{
          // U is obtained locally, just make a copy. Include everything
          // (there is no diagonal block)
          UrowRecv = this->U( LBi( ksup, grid_ ) );
        } // sender is the same as receiver

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "UrowRecv: " <<UrowRecv.size()<< std::endl << std::endl; 
        for(Int jb=0;jb<UrowRecv.size();jb++){ statusOFS<<"jb="<<jb<<" " << UrowRecv[jb] << std::endl; }
#endif

        if( MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){
          std::stringstream     sstm;
          sstm.write( &sstrLcolRecv[0], sizeStmFromLeft );
          std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
          mask[LBlockMask::NZVAL] = 0; // nzval is excluded
          Int numLBlock;
          deserialize( numLBlock, sstm, NO_MASK );
          LcolRecv.resize( numLBlock );
          for( Int ib = 0; ib < numLBlock; ib++ ){
            deserialize( LcolRecv[ib], sstm, mask );
          }
        } // sender is not the same as receiver
        else{
          // L is obtained locally, just make a copy. 
          // Do not include the diagonal block
          std::vector<LBlock>& Lcol =  this->L( LBj( ksup, grid_ ) );
          if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
            LcolRecv.resize( Lcol.size() );
            for( Int ib = 0; ib < Lcol.size(); ib++ ){
              LcolRecv[ib] = Lcol[ib];
            }
          }
          else{
            LcolRecv.resize( Lcol.size() - 1 );
            for( Int ib = 0; ib < Lcol.size() - 1; ib++ ){
              LcolRecv[ib] = Lcol[ib+1];
            }
          }
        } // sender is the same as receiver

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "LcolRecv: " <<LcolRecv.size()<< std::endl << std::endl; 
        for(Int ib=0;ib<LcolRecv.size();ib++){ statusOFS<<"ib="<<ib<<" " << LcolRecv[ib] << std::endl; }
#endif
      } // if I am a receiver

      // Save all the data to be updated for { L( isup, ksup ) | isup > ksup }.
      // The size will be updated in the Gemm phase and the reduce phase
      NumMat<Scalar>  LUpdateBuf;

      // Only the processors received information participate in the Gemm 
      if( isRecvFromAbove_( ksup ) && isRecvFromLeft_( ksup ) ){

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "Main work: Gemm" << std::endl << std::endl; 
#endif


        NumMat<Scalar> AinvBuf, UBuf;
        SelInv_lookup_indexes(ksup,LcolRecv, UrowRecv,AinvBuf,UBuf,LUpdateBuf);

        TIMER_START(Compute_Sinv_LT_GEMM);

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "AinvBuf: " << AinvBuf << std::endl;
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "UBuf: " << UBuf << std::endl;
#endif
        // Gemm for LUpdateBuf = AinvBuf * UBuf^T
        blas::Gemm( 'N', 'T', AinvBuf.m(), UBuf.m(), AinvBuf.n(), SCALAR_MINUS_ONE, 
            AinvBuf.Data(), AinvBuf.m(), 
            UBuf.Data(), UBuf.m(), SCALAR_ZERO,
            LUpdateBuf.Data(), LUpdateBuf.m() ); 


        TIMER_STOP(Compute_Sinv_LT_GEMM);

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<  "LUpdateBuf: " << LUpdateBuf << std::endl;
#endif


      } // if Gemm is to be done locally
      TIMER_STOP(Compute_Sinv_LT);




      // Now all the Isend / Irecv should have finished.
      mpi::Waitall( mpireqsSendToRight );
      mpi::Waitall( mpireqsSendToBelow );


      // Reduce LUpdateBuf across all the processors in the same processor row.

      TIMER_START(Reduce_Sinv_LT);
      NumMat<Scalar> LUpdateBufReduced;

      // Processor column of ksup collects the symbolic data for LUpdateBuf.
      std::vector<Int>  rowLocalPtr;
      std::vector<Int>  blockIdxLocal;
      Int numRowLUpdateBuf;
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
        if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
          rowLocalPtr.resize( Lcol.size() + 1 );
          blockIdxLocal.resize( Lcol.size() );
          rowLocalPtr[0] = 0;
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            rowLocalPtr[ib+1] = rowLocalPtr[ib] + Lcol[ib].numRow;
            blockIdxLocal[ib] = Lcol[ib].blockIdx;
          }
        } // I do not own the diaogonal block
        else{
          rowLocalPtr.resize( Lcol.size() );
          blockIdxLocal.resize( Lcol.size() - 1 );
          rowLocalPtr[0] = 0;
          for( Int ib = 1; ib < Lcol.size(); ib++ ){
            rowLocalPtr[ib] = rowLocalPtr[ib-1] + Lcol[ib].numRow;
            blockIdxLocal[ib-1] = Lcol[ib].blockIdx;
          }
        } // I owns the diagonal block, skip the diagonal block
        numRowLUpdateBuf = *rowLocalPtr.rbegin();
        if( numRowLUpdateBuf > 0 ){
          LUpdateBufReduced.Resize( numRowLUpdateBuf, SuperSize( ksup, super_ ) );
          SetValue( LUpdateBufReduced, SCALAR_ZERO );
        }
      } 

      // Processor column sends the total row dimension to all processors
      // in the same row to prepare for reduce
      MPI_Bcast( &numRowLUpdateBuf, 1, MPI_INT, PCOL( ksup, grid_ ), grid_->rowComm );

      // If LUpdatebuf has not been constructed, resize and fill with zero
      if( numRowLUpdateBuf > 0 ){
        if( LUpdateBuf.m() == 0 && LUpdateBuf.n() == 0 ){
          LUpdateBuf.Resize( numRowLUpdateBuf, SuperSize( ksup, super_ ) );
          // Fill zero is important
          SetValue( LUpdateBuf, SCALAR_ZERO );
        }


        mpi::Reduce( LUpdateBuf.Data(), LUpdateBufReduced.Data(),
            numRowLUpdateBuf * SuperSize( ksup, super_ ), MPI_SUM, 
            PCOL( ksup, grid_ ), grid_->rowComm );


      } // Perform reduce for nonzero block rows in the column of ksup


      TIMER_STOP(Reduce_Sinv_LT);


#if ( _DEBUGlevel_ >= 1 ) 
      if( numRowLUpdateBuf > 0 ){
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<   "LUpdateBufReduced: " <<  LUpdateBufReduced << std::endl << std::endl; 
        }
      }
#endif

#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv::UpdateD");
#endif


      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<   "Update the diagonal block" << std::endl << std::endl; 
#endif
        TIMER_START(Update_Diagonal);
        NumMat<Scalar> DiagBuf( SuperSize( ksup, super_ ), SuperSize( ksup, super_ ) );
        SetValue( DiagBuf, SCALAR_ZERO );
        std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
        if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            blas::Gemm( 'T', 'N', SuperSize( ksup, super_ ), SuperSize( ksup, super_ ), Lcol[ib].numRow, 
                SCALAR_MINUS_ONE, &LUpdateBufReduced( rowLocalPtr[ib], 0 ), LUpdateBufReduced.m(),
                Lcol[ib].nzval.Data(), Lcol[ib].nzval.m(), 
                SCALAR_ONE, DiagBuf.Data(), DiagBuf.m() );
          }
        } // I do not own the diagonal block
        else{
          for( Int ib = 1; ib < Lcol.size(); ib++ ){
            blas::Gemm( 'T', 'N', SuperSize( ksup, super_ ), SuperSize( ksup, super_ ), Lcol[ib].numRow, 
                SCALAR_MINUS_ONE, &LUpdateBufReduced( rowLocalPtr[ib-1], 0 ), LUpdateBufReduced.m(),	
                Lcol[ib].nzval.Data(), Lcol[ib].nzval.m(), 
                SCALAR_ONE, DiagBuf.Data(), DiagBuf.m() );
          }
        } // I own the diagonal block, skip the diagonal block

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<   "DiagBuf: " <<  DiagBuf << std::endl << std::endl; 
#endif

        TIMER_STOP(Update_Diagonal);


        TIMER_START(Reduce_Diagonal);
        NumMat<Scalar> DiagBufReduced( SuperSize( ksup, super_ ), SuperSize( ksup, super_ ) );

        if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
          SetValue( DiagBufReduced, SCALAR_ZERO );

        mpi::Reduce( DiagBuf.Data(), DiagBufReduced.Data(), 
            SuperSize( ksup, super_ ) * SuperSize( ksup, super_ ),
            MPI_SUM, PROW( ksup, grid_ ), grid_->colComm );


        // Add DiagBufReduced to diagonal block.
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          LBlock&  LB = this->L( LBj( ksup, grid_ ) )[0];

#if ( _DEBUGlevel_ >= 1 )
          statusOFS<< std::endl << "ORIGINAL ["<<ksup<<"] "<<   "DiagBufReduced: " << DiagBufReduced << std::endl << std::endl; 
#endif
          // Symmetrize LB
          blas::Axpy( LB.numRow * LB.numCol, SCALAR_ONE, DiagBufReduced.Data(),
              1, LB.nzval.Data(), 1 );

          Symmetrize( LB.nzval );

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "ORIGINAL ["<<ksup<<"] "<<   "Diag of Ainv: " << LB.nzval << std::endl << std::endl; 
#endif
        }

        TIMER_STOP(Reduce_Diagonal);
      } // Update the diagonal in the processor column of ksup. All processors participate


#ifndef _RELEASE_
      PopCallStack();
#endif



#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv::UpdateU");
#endif


      //compute the number of requests
      Int sendCount= CountSendToCrossDiagonal(ksup);
      Int recvCount= CountRecvFromCrossDiagonal(ksup);

      std::vector<MPI_Request > mpiReqsSendCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > mpiReqsSizeSendCD(sendCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > sstrLcolSendCD(sendCount);
      std::vector<int > sstrLcolSizeSendCD(sendCount);

      std::vector<MPI_Request > mpiReqsRecvCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > mpiReqsSizeRecvCD(recvCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > sstrLcolRecvCD(recvCount);
      std::vector<int > sstrLcolSizeRecvCD(recvCount);


      TIMER_START(Update_U);

      // Send LUpdateBufReduced to the cross diagonal blocks. 
      // NOTE: This assumes square processor grid
      TIMER_START(Send_L_CrossDiag);

      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) && isSendToCrossDiagonal_(grid_->numProcCol, ksup ) ){

        Int sendIdx = 0;
        for(Int dstCol = 0; dstCol<grid_->numProcCol; dstCol++){
          if(isSendToCrossDiagonal_(dstCol,ksup) ){
            Int dest = PNUM(PROW(ksup,grid_),dstCol,grid_);

            if( MYPROC( grid_ ) != dest	){
              std::stringstream sstm;
              std::vector<char> & sstrLcolSend = sstrLcolSendCD[sendIdx];
              Int & sstrSize = sstrLcolSizeSendCD[sendIdx];
              MPI_Request & mpiReqSizeSend = mpiReqsSizeSendCD[sendIdx];
              MPI_Request & mpiReqSend = mpiReqsSendCD[sendIdx];


              serialize( rowLocalPtr, sstm, NO_MASK );
              serialize( blockIdxLocal, sstm, NO_MASK );
              serialize( LUpdateBufReduced, sstm, NO_MASK );

              sstrLcolSend.resize( Size(sstm) );
              sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
              sstrSize = sstrLcolSend.size();

#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<"BCAST ["<<ksup<<"] P"<<MYPROC(grid_)<<" ("<<MYROW(grid_)<<","<<MYCOL(grid_)<<") ---> LBj("<<ksup<<") ---> P"<<dest<<std::endl;
              statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "rowLocalPtr:" << rowLocalPtr << std::endl << std::endl; 
              statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "blockIdxLocal:" << blockIdxLocal << std::endl << std::endl; 
#endif
              MPI_Isend( &sstrSize, 1, MPI_INT, dest, SELINV_TAG_L_SIZE, grid_->comm, &mpiReqSizeSend );
              MPI_Isend( (void*)&sstrLcolSend[0], sstrSize, MPI_BYTE, dest, SELINV_TAG_L_CONTENT, grid_->comm, &mpiReqSend );


              sendIdx++;
            }
          }
        }

      } // sender
      TIMER_STOP(Send_L_CrossDiag);






      //Do Irecv for sizes
      //If I'm a receiver
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) && isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) ){
        Int recvIdx=0;
        for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
          if(isRecvFromCrossDiagonal_(srcRow,ksup) ){
            Int src = PNUM(srcRow,PCOL(ksup,grid_),grid_);
            if( MYPROC( grid_ ) != src ){
              Int & sstrSize = sstrLcolSizeRecvCD[recvIdx];
              MPI_Request & mpiReqSizeRecv = mpiReqsSizeRecvCD[recvIdx];
              MPI_Irecv( &sstrSize, 1, MPI_INT, src, SELINV_TAG_L_SIZE, grid_->comm, &mpiReqSizeRecv );
              recvIdx++;
            }
          }
        }
      }//end if I'm a receiver

      //waitall sizes
      mpi::Waitall(mpiReqsSizeRecvCD);

      //Allocate content and do Irecv
      //If I'm a receiver
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) && isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) ){
        Int recvIdx=0;
        for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
          if(isRecvFromCrossDiagonal_(srcRow,ksup) ){
            Int src = PNUM(srcRow,PCOL(ksup,grid_),grid_);
            if( MYPROC( grid_ ) != src ){
              Int & sstrSize = sstrLcolSizeRecvCD[recvIdx];
              std::vector<char> & sstrLcolRecv = sstrLcolRecvCD[recvIdx];
              MPI_Request & mpiReqRecv = mpiReqsRecvCD[recvIdx];
              sstrLcolRecv.resize( sstrSize);
              MPI_Irecv( (void*)&sstrLcolRecv[0], sstrSize, MPI_BYTE, src, SELINV_TAG_L_CONTENT, grid_->comm, &mpiReqRecv );
              recvIdx++;
            }
          }
        }
      }//end if I'm a receiver


      //waitall content
      mpi::Waitall(mpiReqsRecvCD);


      //Do the work
      if( MYROW( grid_ ) == PROW( ksup, grid_ ) && isRecvFromCrossDiagonal_(grid_->numProcRow, ksup ) ){

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl <<  "BCAST ["<<ksup<<"] "<<  "Update the upper triangular block" << std::endl << std::endl; 
        statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "blockIdxLocal:" << blockIdxLocal << std::endl << std::endl; 
        statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "rowLocalPtr:" << rowLocalPtr << std::endl << std::endl; 
#endif


        std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
        std::vector<bool> isBlockFound(Urow.size(),false);

        Int recvIdx=0;
        for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
          if(isRecvFromCrossDiagonal_(srcRow,ksup) ){
            Int src = PNUM(srcRow,PCOL(ksup,grid_),grid_);
            TIMER_START(Recv_L_CrossDiag);
            std::vector<Int> rowLocalPtrRecv;
            std::vector<Int> blockIdxLocalRecv;
            NumMat<Scalar> UUpdateBuf;


            if( MYPROC( grid_ ) != src ){
              std::stringstream sstm;
              Int & sstrSize = sstrLcolSizeRecvCD[recvIdx];
              std::vector<char> & sstrLcolRecv = sstrLcolRecvCD[recvIdx];
              sstm.write( &sstrLcolRecv[0], sstrSize );

              deserialize( rowLocalPtrRecv, sstm, NO_MASK );
              deserialize( blockIdxLocalRecv, sstm, NO_MASK );
              deserialize( UUpdateBuf, sstm, NO_MASK );	
              recvIdx++;
            } // sender is not the same as receiver
            else{
              rowLocalPtrRecv   = rowLocalPtr;
              blockIdxLocalRecv = blockIdxLocal;
              UUpdateBuf = LUpdateBufReduced;
            } // sender is the same as receiver



            TIMER_STOP(Recv_L_CrossDiag);

#if ( _DEBUGlevel_ >= 1 )
            statusOFS<<"BCAST ["<<ksup<<"] P"<<MYPROC(grid_)<<" ("<<MYROW(grid_)<<","<<MYCOL(grid_)<<") <--- LBj("<<ksup<<") <--- P"<<src<<std::endl;
            statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "rowLocalPtrRecv:" << rowLocalPtrRecv << std::endl << std::endl; 
            //              statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "UUpdateBuf:" << UUpdateBuf << std::endl << std::endl; 
            statusOFS << std::endl << "BCAST ["<<ksup<<"] "<<   "blockIdxLocalRecv:" << blockIdxLocalRecv << std::endl << std::endl; 
#endif




            // Update U
            for( Int ib = 0; ib < blockIdxLocalRecv.size(); ib++ ){
              for( Int jb = 0; jb < Urow.size(); jb++ ){
                UBlock& UB = Urow[jb];
                if( UB.blockIdx == blockIdxLocalRecv[ib] ){
                  NumMat<Scalar> Ltmp ( UB.numCol, UB.numRow );
                  lapack::Lacpy( 'A', Ltmp.m(), Ltmp.n(), 
                      &UUpdateBuf( rowLocalPtrRecv[ib], 0 ),
                      UUpdateBuf.m(), Ltmp.Data(), Ltmp.m() );
                  isBlockFound[jb] = true;
                  Transpose( Ltmp, UB.nzval );
                  break;
                }
              }
            }
          }
        }


        for( Int jb = 0; jb < Urow.size(); jb++ ){
          UBlock& UB = Urow[jb];
          if( !isBlockFound[jb] ){
            throw std::logic_error( "UBlock cannot find its update. Something is seriously wrong." );
          }
        }
      } // receiver
      TIMER_STOP(Update_U);

      mpi::Waitall(mpiReqsSizeSendCD);
      mpi::Waitall(mpiReqsSendCD);


#ifndef _RELEASE_
      PopCallStack();
#endif

#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv::UpdateLFinal");
#endif

      TIMER_START(Update_L);



      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) && numRowLUpdateBuf > 0 ){
        std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
        if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            LBlock& LB = Lcol[ib];
            lapack::Lacpy( 'A', LB.numRow, LB.numCol, &LUpdateBufReduced( rowLocalPtr[ib], 0 ),
                LUpdateBufReduced.m(), LB.nzval.Data(), LB.numRow );
          }
        } // I do not own the diagonal block
        else{
          for( Int ib = 1; ib < Lcol.size(); ib++ ){
            LBlock& LB = Lcol[ib];
            lapack::Lacpy( 'A', LB.numRow, LB.numCol, &LUpdateBufReduced( rowLocalPtr[ib-1], 0 ),
                LUpdateBufReduced.m(), LB.nzval.Data(), LB.numRow );
          }
        } // I owns the diagonal block
      } // Finish updating L	


      TIMER_STOP(Update_L);


#ifndef _RELEASE_
      PopCallStack();
#endif

      MPI_Barrier( grid_-> comm );

    } // for (ksup) : Main loop

    TIMER_STOP(SelInv);


#ifndef _RELEASE_
    PopCallStack();
#endif

    return ;
  } 		// -----  end of method PMatrix::SelInv  ----- 

  void PMatrix::PreSelInv	(  )
  {
#ifndef _RELEASE_
    PushCallStack("PMatrix::PreSelInv");
#endif

    Int numSuper = this->NumSuper(); 

#ifndef _RELEASE_
    PushCallStack("L(i,k) <- L(i,k) * L(k,k)^{-1}");
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "L(i,k) <- L(i,k) * L(k,k)^{-1}" << std::endl << std::endl; 
#endif
    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        // Broadcast the diagonal L block
        NumMat<Scalar> nzvalLDiag;
        std::vector<LBlock>& Lcol = this->L( LBj( ksup, grid_ ) );
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          nzvalLDiag = Lcol[0].nzval;
          if( nzvalLDiag.m() != SuperSize(ksup, super_) ||
              nzvalLDiag.n() != SuperSize(ksup, super_) ){
            throw std::runtime_error( "The size of the diagonal block of L is wrong." );
          }
        } // Owns the diagonal block
        else
        {
          nzvalLDiag.Resize( SuperSize(ksup, super_), SuperSize(ksup, super_) );
        }
        MPI_Bcast( (void*)nzvalLDiag.Data(), 
            sizeof(Scalar) * SuperSize(ksup, super_) * SuperSize(ksup, super_),
            MPI_BYTE, PROW( ksup, grid_ ), grid_->colComm );

        // Triangular solve
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          LBlock& LB = Lcol[ib];
          if( LB.blockIdx > ksup ){
#if ( _DEBUGlevel_ >= 2 )
            // Check the correctness of the triangular solve for the first local column
            if( LBj( ksup, grid_ ) == 0 ){
              statusOFS << "Diag   L(" << ksup << ", " << ksup << "): " << nzvalLDiag << std::endl;
              statusOFS << "Before solve L(" << LB.blockIdx << ", " << ksup << "): " << LB.nzval << std::endl;
            }
#endif
            blas::Trsm( 'R', 'L', 'N', 'U', LB.numRow, LB.numCol, SCALAR_ONE,
                nzvalLDiag.Data(), LB.numCol, LB.nzval.Data(), LB.numRow );
#if ( _DEBUGlevel_ >= 2 )
            // Check the correctness of the triangular solve for the first local column
            if( LBj( ksup, grid_ ) == 0 ){
              statusOFS << "After solve  L(" << LB.blockIdx << ", " << ksup << "): " << LB.nzval << std::endl;
            }
#endif
          }
        }
      } // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
    } // for (ksup)


#ifndef _RELEASE_
    PopCallStack();
#endif


#ifndef _RELEASE_
    PushCallStack("U(k,i) <- L(i,k)");
#endif
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "U(k,i) <- L(i,k)" << std::endl << std::endl; 
#endif

    for( Int ksup = 0; ksup < numSuper; ksup++ ){
      Int ksupProcRow = PROW( ksup, grid_ );
      Int ksupProcCol = PCOL( ksup, grid_ );

      Int sendCount = CountSendToCrossDiagonal(ksup);
      Int recvCount = CountRecvFromCrossDiagonal(ksup);

      std::vector<MPI_Request > arrMpiReqsSend(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeSend(sendCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrLcolSend(sendCount);
      std::vector<Int > arrSstrLcolSizeSend(sendCount);

      std::vector<MPI_Request > arrMpiReqsRecv(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeRecv(recvCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrLcolRecv(recvCount);
      std::vector<Int > arrSstrLcolSizeRecv(recvCount);



      // Sender
      if( isSendToCrossDiagonal_(grid_->numProcCol,ksup) ){
#if ( _DEBUGlevel_ >= 1 )
        statusOFS<<"["<<ksup<<"] P"<<MYPROC(grid_)<<" should send to "<<CountSendToCrossDiagonal(ksup)<<" processors"<<std::endl;
#endif

        Int sendIdx = 0;
        for(Int dstCol = 0; dstCol<grid_->numProcCol; dstCol++){
          if(isSendToCrossDiagonal_(dstCol,ksup) ){
            Int dst = PNUM(PROW(ksup,grid_),dstCol,grid_);
            if(MYPROC(grid_)!= dst){
              // Pack L data
              std::stringstream sstm;
              std::vector<char> & sstrLcolSend = arrSstrLcolSend[sendIdx];
              Int & sstrSize = arrSstrLcolSizeSend[sendIdx];
              MPI_Request & mpiReqSend = arrMpiReqsSend[sendIdx];
              MPI_Request & mpiReqSizeSend = arrMpiReqsSizeSend[sendIdx];

              std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
              std::vector<LBlock>&  Lcol = this->L( LBj(ksup, grid_) );
              // All blocks except for the diagonal block are to be sent right
              //TODO not true > this is a scatter operation ! Can we know the destination ?

              Int count = 0;
              if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
                for( Int ib = 1; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > ksup &&  (Lcol[ib].blockIdx % grid_->numProcCol) == dstCol  ){
                    count++;
                  }
                }
              }
              else{

                for( Int ib = 0; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > ksup &&  (Lcol[ib].blockIdx % grid_->numProcCol) == dstCol  ){
                    count++;
                  }
                }
              }

              serialize( (Int)count, sstm, NO_MASK );

              for( Int ib = 0; ib < Lcol.size(); ib++ ){
                if( Lcol[ib].blockIdx > ksup &&  (Lcol[ib].blockIdx % grid_->numProcCol) == dstCol  ){ 
                  //                if( Lcol[ib].blockIdx > ksup ){ 
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS<<"["<<ksup<<"] SEND contains "<<Lcol[ib].blockIdx<< " which corresponds to "<<GBj(ib,grid_)<<std::endl;
#endif
                  serialize( Lcol[ib], sstm, mask );
                }
                }

                sstrLcolSend.resize( Size(sstm) );
                sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
                sstrSize = sstrLcolSend.size();



                // Send/Recv is possible here due to the one to one correspondence
                // in the case of square processor grid

#if ( _DEBUGlevel_ >= 1 )
                statusOFS<<"["<<ksup<<"] P"<<MYPROC(grid_)<<" ("<<MYROW(grid_)<<","<<MYCOL(grid_)<<") ---> LBj("<<ksup<<")="<<LBj(ksup,grid_)<<" ---> P"<<dst<<" ("<<ksupProcRow<<","<<dstCol<<")"<<std::endl;
#endif
                MPI_Isend( &sstrSize, 1, MPI_INT, dst, SELINV_TAG_D_SIZE, grid_->comm, &mpiReqSizeSend );
                MPI_Isend( (void*)&sstrLcolSend[0], sstrSize, MPI_BYTE, dst, SELINV_TAG_D_CONTENT, grid_->comm, &mpiReqSend );
                //mpi::Send( sstm, dst,SELINV_TAG_D_SIZE, SELINV_TAG_D_CONTENT, grid_->comm );

                sendIdx++;
              } // if I am a sender
            }
          }
        }





        // Receiver
        if( isRecvFromCrossDiagonal_(grid_->numProcRow,ksup) ){


#if ( _DEBUGlevel_ >= 1 )
          statusOFS<<"["<<ksup<<"] P"<<MYPROC(grid_)<<" should receive from "<<CountRecvFromCrossDiagonal(ksup)<<" processors"<<std::endl;
#endif


          std::vector<UBlock>& Urow = this->U( LBi( ksup, grid_ ) );
          std::vector<bool> isBlockFound(Urow.size(),false);


          Int recvIdx = 0;
          //receive size first
          for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
            if(isRecvFromCrossDiagonal_(srcRow,ksup) ){
              std::vector<LBlock> LcolRecv;
              Int src = PNUM(srcRow,PCOL(ksup,grid_),grid_);
              if(MYPROC(grid_)!= src){
                MPI_Request & mpiReqSizeRecv = arrMpiReqsSizeRecv[recvIdx];
                Int & sstrSize = arrSstrLcolSizeRecv[recvIdx];

                MPI_Irecv( &sstrSize, 1, MPI_INT, src, SELINV_TAG_D_SIZE, grid_->comm, &mpiReqSizeRecv );

                recvIdx++;
              }
            }
          }

          mpi::Waitall(arrMpiReqsSizeRecv);



          //receive content
          recvIdx = 0;
          for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
            if(isRecvFromCrossDiagonal_(srcRow,ksup) ){
              std::vector<LBlock> LcolRecv;
              Int src = PNUM(srcRow,PCOL(ksup,grid_),grid_);
              if(MYPROC(grid_)!= src){
                MPI_Request & mpiReqRecv = arrMpiReqsRecv[recvIdx];
                Int & sstrSize = arrSstrLcolSizeRecv[recvIdx];
                std::vector<char> & sstrLcolRecv = arrSstrLcolRecv[recvIdx];
                sstrLcolRecv.resize(sstrSize);

                MPI_Irecv( &sstrLcolRecv[0], sstrSize, MPI_BYTE, src, SELINV_TAG_D_CONTENT, grid_->comm, &mpiReqRecv );

                recvIdx++;
              }
            }
          }

          mpi::Waitall(arrMpiReqsRecv);



          //Process the content
          recvIdx = 0;
          for(Int srcRow = 0; srcRow<grid_->numProcRow; srcRow++){
            if(isRecvFromCrossDiagonal_(srcRow,ksup) ){
              std::vector<LBlock> LcolRecv;
              Int src = PNUM(srcRow,PCOL(ksup,grid_),grid_);
              if(MYPROC(grid_)!= src){

                Int & sstrSize = arrSstrLcolSizeRecv[recvIdx];
                std::vector<char> & sstrLcolRecv = arrSstrLcolRecv[recvIdx];
                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS<<"["<<ksup<<"] P"<<MYPROC(grid_)<<" ("<<MYROW(grid_)<<","<<MYCOL(grid_)<<") <--- LBj("<<ksup<<") <--- P"<<src<<" ("<<srcRow<<","<<ksupProcCol<<")"<<std::endl;
#endif


                sstm.write( &sstrLcolRecv[0], sstrSize );

                // Unpack L data.  
                Int numLBlock;
                std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
                deserialize( numLBlock, sstm, NO_MASK );
                LcolRecv.resize(numLBlock);
                for( Int ib = 0; ib < numLBlock; ib++ ){
                  deserialize( LcolRecv[ib], sstm, mask );
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS<<"["<<ksup<<"] RECV contains "<<LcolRecv[ib].blockIdx<< " which corresponds to "<< ib * grid_->numProcRow + srcRow; // <<std::endl;
                  //                  statusOFS<<" L is on row "<< srcRow <<" whereas U is on col "<<((ib * grid_->numProcRow + srcRow)/grid_->numProcCol)%grid_->numProcCol <<std::endl;
                  statusOFS<<" L is on row "<< srcRow <<" whereas U is on col "<< LcolRecv[ib].blockIdx % grid_->numProcCol <<std::endl;
#endif
                }


                recvIdx++;

              } // sender is not the same as receiver
              else{
                // L is obtained locally, just make a copy. Do not include the diagonal block
                std::vector<LBlock>& Lcol = this->L( LBj( ksup, grid_ ) );
                if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
                  LcolRecv.resize( Lcol.size() );
                  for( Int ib = 0; ib < Lcol.size(); ib++ ){
                    LcolRecv[ib] = Lcol[ib];
                  }
                }
                else{
                  LcolRecv.resize( Lcol.size() - 1 );
                  for( Int ib = 0; ib < Lcol.size() - 1; ib++ ){
                    LcolRecv[ib] = Lcol[ib+1];
                  }
                }
              } // sender is the same as receiver

              // Update U
              // Make sure that the size of L and the corresponding U blocks match.
              for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
                LBlock& LB = LcolRecv[ib];
                if( LB.blockIdx <= ksup ){
                  throw std::logic_error( "LcolRecv contains the wrong blocks." );
                }
                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  UBlock&  UB = Urow[jb];
                  if( LB.blockIdx == UB.blockIdx ){
                    // Compare size
                    if( LB.numRow != UB.numCol || LB.numCol != UB.numRow ){
                      std::ostringstream msg;
                      msg << "LB(" << LB.blockIdx << ", " << ksup << ") and UB(" 
                        << ksup << ", " << UB.blockIdx << ")	do not share the same size." << std::endl
                        << "LB: " << LB.numRow << " x " << LB.numCol << std::endl
                        << "UB: " << UB.numRow << " x " << UB.numCol << std::endl;
                      throw std::runtime_error( msg.str().c_str() );
                    }

                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    UB.cols = LB.rows;
                    Transpose( LB.nzval, UB.nzval );

#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<ksup<<"] USING "<<LB.blockIdx<< std::endl;
#endif
                    isBlockFound[jb] = true;
                    break;
                  } // if( LB.blockIdx == UB.blockIdx )
                } // for (jb)
              } // for (ib)
            }
          }

          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock& UB = Urow[jb];
            if( !isBlockFound[jb] ){
              throw std::logic_error( "UBlock cannot find its update. Something is seriously wrong." );
            }
          }



        } // if I am a receiver


        //Wait until every receiver is done
        mpi::Waitall(arrMpiReqsSizeSend);
        mpi::Waitall(arrMpiReqsSend);


      } // for (ksup)

#ifndef _RELEASE_
      PopCallStack();
#endif

#ifndef _RELEASE_
      PushCallStack("L(i,i) <- [L(k,k) * U(k,k)]^{-1} ");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "L(i,i) <- [L(k,k) * U(k,k)]^{-1}" << std::endl << std::endl; 
#endif

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) &&
            MYCOL( grid_ ) == PCOL( ksup, grid_ )	){
          IntNumVec ipiv( SuperSize( ksup, super_ ) );
          // Note that the pivoting vector ipiv should follow the FORTRAN
          // notation by adding the +1
          for(Int i = 0; i < SuperSize( ksup, super_ ); i++){
            ipiv[i] = i + 1;
          }
          LBlock& LB = (this->L( LBj( ksup, grid_ ) ))[0];
#if ( _DEBUGlevel_ >= 2 )
          // Check the correctness of the matrix inversion for the first local column
          statusOFS << "Factorized A (" << ksup << ", " << ksup << "): " << LB.nzval << std::endl;
#endif
          lapack::Getri( SuperSize( ksup, super_ ), LB.nzval.Data(), 
              SuperSize( ksup, super_ ), ipiv.Data() );

          // Symmetrize the diagonal block
          Symmetrize( LB.nzval );
#if ( _DEBUGlevel_ >= 2 )
          // Check the correctness of the matrix inversion for the first local column
          statusOFS << "Inversed   A (" << ksup << ", " << ksup << "): " << LB.nzval << std::endl;
#endif
        } // if I need to inverse the diagonal block
      } // for (ksup)


#ifndef _RELEASE_
      PopCallStack();
#endif



#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::PreSelInv  ----- 


#ifdef SANITY_CHECK
    void PMatrix::GetColumn	( Int colIdx,  NumVec<Scalar>& col )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::GetColumn");
#endif
      Int numSuper = this->NumSuper(); 

      Int numRow = this->NumCol();
      const IntNumVec& permInv = super_->permInv;
      NumVec<Scalar> colLocal( numRow );
      SetValue( colLocal, SCALAR_ZERO );

      col.Resize( numRow );
      SetValue( col, SCALAR_ZERO );

      //find the supernode containing colIdx
      Int ksup = BlockIdx(colIdx,super_);

      Int localColIdx = colIdx - FirstBlockCol( ksup, super_ );

      if ( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
        std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
          statusOFS<<Lcol[ib]<<std::endl;
          Int j = localColIdx;
          for( Int i = 0; i < Lcol[ib].numRow; i++ ){
            colLocal( permInv( Lcol[ib].rows(i) ) ) = Lcol[ib].nzval( i, j );
          }         
        }
      }

      mpi::Allreduce( colLocal.Data(), col.Data(), numRow, MPI_SUM, grid_->comm );

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::GetColumn  ----- 
#endif





    void PMatrix::GetDiagonal	( NumVec<Scalar>& diag )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::GetDiagonal");
#endif
      Int numSuper = this->NumSuper(); 

      Int numCol = this->NumCol();
      const IntNumVec& permInv = super_->permInv;

      NumVec<Scalar> diagLocal( numCol );
      SetValue( diagLocal, SCALAR_ZERO );

      diag.Resize( numCol );
      SetValue( diag, SCALAR_ZERO );


      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // I own the diagonal block	
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) &&
            MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          LBlock& LB = this->L( LBj( ksup, grid_ ) )[0];
          for( Int i = 0; i < LB.numRow; i++ ){
            diagLocal( permInv( LB.rows(i) ) ) = LB.nzval( i, i );
          }
        }
      }

      // All processors own diag
      mpi::Allreduce( diagLocal.Data(), diag.Data(), numCol, MPI_SUM, grid_->comm );

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::GetDiagonal  ----- 


#ifdef SANITY_CHECK
    void PMatrix::CompareDiagonal	( PMatrix & Ref, SelInvErrors & errors)
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::CompareDiagonal");
#endif
      Int numSuper = this->NumSuper(); 

      Int numCol = this->NumCol();

      SelInvErrors localErrors;

      Int maxIb, maxI, maxJ, maxK =0;


      for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
        // I own the diagonal block	
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) &&
            MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){



          std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
          std::vector<LBlock>&  LcolRef = Ref.L( LBj( ksup, grid_ ) );

          if( MYROW( grid_ ) == PROW( ksup, grid_ )){
            for( Int ib = 0; ib < std::min(1,(int)Lcol.size()); ib++ ){
              Real fNormRef = 0.0;
              Real fNorm = 0.0;


              for( Int j = 0; j < Lcol[ib].numCol; j++ ){
                Real colNorm = 0.0;
                Real colNormRef = 0.0;


                for( Int i = 0; i < Lcol[ib].numRow; i++ ){
                  //Compute the 2-norm of current columns 
                  colNorm += abs( pow(Lcol[ib].nzval(i,j),2) );
                  colNormRef += abs( pow(LcolRef[ib].nzval(i,j),2) );
                }

                Int gj = FirstBlockCol(ksup,super_) + j;


                //Compute max row-wise errors
                Real cwiseAbsError = abs(colNorm - colNormRef);
                Real cwiseRelError = cwiseAbsError / abs(colNormRef);
                if( cwiseRelError > localErrors.MaxCwiseRelError.Value){localErrors.MaxCwiseRelError.Set(cwiseRelError,ksup,ib,-1,gj); localErrors.CorrCwiseAbsError.Set(cwiseAbsError,ksup,ib,-1,gj); }
                if( cwiseAbsError > localErrors.MaxCwiseAbsError.Value){localErrors.MaxCwiseAbsError.Set(cwiseAbsError,ksup,ib,-1,gj); }



              }
              for( Int i = 0; i < Lcol[ib].numRow; i++ ){
                Real rowNorm = 0.0;
                Real rowNormRef = 0.0;
                Int gi = FirstBlockRow(ksup,super_) + i;

                for( Int j = 0; j < Lcol[ib].numCol; j++ ){

                  Int gj = FirstBlockCol(ksup,super_) + j;
                  //Compute the 2-norm of current row      
                  rowNorm += abs( pow(Lcol[ib].nzval(i,j),2) );
                  rowNormRef += abs( pow(LcolRef[ib].nzval(i,j),2) );

                  std::stringstream msg;
                  Real absError = abs(Lcol[ib].nzval(i,j)-LcolRef[ib].nzval(i,j));
                  Real relError = absError/abs(LcolRef[ib].nzval(i,j));
                  msg<< "["<<ksup<< "] D Block "<<ib<< " Row "<<i<<" Col "<<j<<" is wrong : "<< Lcol[ib].nzval(i,j) << " vs "<<LcolRef[ib].nzval(i,j)<< " relative error is "<< relError <<std::endl; 

                  if( relError > localErrors.MaxRelError.Value){localErrors.MaxRelError.Set(relError,ksup,ib,gi,gj); localErrors.CorrAbsError.Set(absError,ksup,ib,gi,gj); }
                  if( absError > localErrors.MaxAbsError.Value){localErrors.MaxAbsError.Set(absError,ksup,ib,gi,gj); }

                  if(relError > SANITY_PRECISION){
                    statusOFS<<msg;
                  }
                }

                //Compute max row-wise errors
                Real rwiseAbsError = abs(rowNorm - rowNormRef);
                Real rwiseRelError = rwiseAbsError / abs(rowNormRef);

                if( rwiseRelError > localErrors.MaxRwiseRelError.Value){localErrors.MaxRwiseRelError.Set(rwiseRelError,ksup,ib,gi,-1); localErrors.CorrRwiseAbsError.Set(rwiseAbsError,ksup,ib,gi,-1); }
                if( rwiseAbsError > localErrors.MaxRwiseAbsError.Value){localErrors.MaxRwiseAbsError.Set(rwiseAbsError,ksup,ib,gi,-1); }
                //Compute the Frobenius norm
                fNorm += rowNorm;
                fNormRef += rowNormRef;
              }


              Real nwiseAbsError = abs(fNorm - fNormRef);
              Real nwiseRelError = nwiseAbsError / abs(fNormRef);

              if( nwiseRelError > localErrors.MaxNwiseRelError.Value){localErrors.MaxNwiseRelError.Set(nwiseRelError,ksup,ib,-1,-1); localErrors.CorrNwiseAbsError.Set(nwiseAbsError,ksup,ib,-1,-1); }
              if( nwiseAbsError > localErrors.MaxNwiseAbsError.Value){localErrors.MaxNwiseAbsError.Set(nwiseAbsError,ksup,ib,-1,-1); }
            }
          }
        }
      }

      struct{ double val; int rank; } lmaxRelError,gmaxRelError;

      //Get max element-wise errors
      lmaxRelError.val = localErrors.MaxRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrAbsError = localErrors.CorrAbsError; errors.MaxRelError = localErrors.MaxRelError; }
      MPI_Bcast(&errors.CorrAbsError, sizeof(errors.CorrAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxRelError, sizeof(errors.MaxRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxAbsError = localErrors.MaxAbsError; }
      MPI_Bcast(&errors.MaxAbsError, sizeof(errors.MaxAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );


      //Get max norm-wise errors
      lmaxRelError.val = localErrors.MaxNwiseRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrNwiseAbsError = localErrors.CorrNwiseAbsError; errors.MaxNwiseRelError = localErrors.MaxNwiseRelError; }
      MPI_Bcast(&errors.CorrNwiseAbsError, sizeof(errors.CorrNwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxNwiseRelError, sizeof(errors.MaxNwiseRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxNwiseAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxNwiseAbsError = localErrors.MaxNwiseAbsError; }
      MPI_Bcast(&errors.MaxNwiseAbsError, sizeof(errors.MaxNwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      //Get max row-wise errors
      lmaxRelError.val = localErrors.MaxRwiseRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrRwiseAbsError = localErrors.CorrRwiseAbsError; errors.MaxRwiseRelError = localErrors.MaxRwiseRelError; }
      MPI_Bcast(&errors.CorrRwiseAbsError, sizeof(errors.CorrRwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxRwiseRelError, sizeof(errors.MaxRwiseRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxRwiseAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxRwiseAbsError = localErrors.MaxRwiseAbsError; }
      MPI_Bcast(&errors.MaxRwiseAbsError, sizeof(errors.MaxRwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      //Get max column-wise errors
      lmaxRelError.val = localErrors.MaxCwiseRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrCwiseAbsError = localErrors.CorrCwiseAbsError; errors.MaxCwiseRelError = localErrors.MaxCwiseRelError; }
      MPI_Bcast(&errors.CorrCwiseAbsError, sizeof(errors.CorrCwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxCwiseRelError, sizeof(errors.MaxCwiseRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxCwiseAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxCwiseAbsError = localErrors.MaxCwiseAbsError; }
      MPI_Bcast(&errors.MaxCwiseAbsError, sizeof(errors.MaxCwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );


#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::CompareDiagonal  ----- 








    void PMatrix::CompareOffDiagonal	( PMatrix & Ref,SelInvErrors & errors)
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::CompareOffDiagonal");
#endif
      Int numSuper = this->NumSuper(); 

      Int numCol = this->NumCol();
      Real maxRelError = 0;
      Real corrAbsError = 0;
      Real maxAbsError = 0;

      SelInvErrors localErrors;

      Int maxIb, maxI, maxJ, maxK =0;


      for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
        if( MYROW( grid_ ) == PROW( ksup, grid_ )){
          std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
          std::vector<UBlock>&  UrowRef = Ref.U( LBi( ksup, grid_ ) );
          for( Int ib = 0; ib < Urow.size(); ib++ ){
            Real fNormRef = 0.0;
            Real fNorm = 0.0;
            for( Int j = 0; j < Urow[ib].numCol; j++ ){
              Int gj = FirstBlockCol(ksup,super_) + j;
              Real colNorm = 0.0;
              Real colNormRef = 0.0;

              for( Int i = 0; i < Urow[ib].numRow; i++ ){
                //Compute the 2-norm of current columns 
                colNorm += abs( pow(Urow[ib].nzval(i,j),2) );
                colNormRef += abs( pow(UrowRef[ib].nzval(i,j),2) );
              }

              //Compute max row-wise errors
              Real cwiseAbsError = abs(colNorm - colNormRef);
              Real cwiseRelError = cwiseAbsError / abs(colNormRef);
              if( cwiseRelError > localErrors.MaxCwiseRelError.Value){
                localErrors.MaxCwiseRelError.Set(cwiseRelError,ksup,ib,-1,gj);
                localErrors.CorrCwiseAbsError.Set(cwiseAbsError,ksup,ib,-1,gj);
              }
              if( cwiseAbsError > localErrors.MaxCwiseAbsError.Value){
                localErrors.MaxCwiseAbsError.Set(cwiseAbsError,ksup,ib,-1,gj);
              }
            }

            for( Int i = 0; i < Urow[ib].numRow; i++ ){
              Int gi = FirstBlockRow(ksup,super_) + i;
              Real rowNorm = 0.0;
              Real rowNormRef = 0.0;
              for( Int j = 0; j < Urow[ib].numCol; j++ ){
                Int gj = FirstBlockCol(ksup,super_) + j;
                //Compute the 2-norm of current row      
                rowNorm += abs( pow(Urow[ib].nzval(i,j),2) );
                rowNormRef += abs( pow(UrowRef[ib].nzval(i,j),2) );
                std::stringstream msg;
                Real absError = abs(Urow[ib].nzval(i,j)-UrowRef[ib].nzval(i,j));
                Real relError = absError/abs(UrowRef[ib].nzval(i,j));

                msg<< "["<<ksup<< "] U Block "<<ib<< " Row "<<i<<" Col "<<j<<" is wrong : "<< Urow[ib].nzval(i,j) << " vs "<<UrowRef[ib].nzval(i,j)<< " relative error is "<< relError <<std::endl; 
                if( relError > localErrors.MaxRelError.Value){
                  localErrors.MaxRelError.Set(relError,ksup,ib,gi,gj); 
                  localErrors.CorrAbsError.Set(absError,ksup,ib,gi,gj); 
                }
                if( absError > localErrors.MaxAbsError.Value){
                  localErrors.MaxAbsError.Set(absError,ksup,ib,gi,gj); 
                }
                if(relError > SANITY_PRECISION){
                  statusOFS<<msg;
                }
              }
              //Compute max row-wise errors
              Real rwiseAbsError = abs(rowNorm - rowNormRef);
              Real rwiseRelError = rwiseAbsError / abs(rowNormRef);
              if( rwiseRelError > localErrors.MaxRwiseRelError.Value){
                localErrors.MaxRwiseRelError.Set(rwiseRelError,ksup,ib,gi,-1);
                localErrors.CorrRwiseAbsError.Set(rwiseAbsError,ksup,ib,gi,-1); 
              }
              if( rwiseAbsError > localErrors.MaxRwiseAbsError.Value){
                localErrors.MaxRwiseAbsError.Set(rwiseAbsError,ksup,ib,gi,-1); 
              }
              //Compute the Frobenius norm
              fNorm += rowNorm;
              fNormRef += rowNormRef;
            }         
            //Compute max norm-wise errors
            Real nwiseAbsError = abs(fNorm - fNormRef);
            Real nwiseRelError = nwiseAbsError / abs(fNormRef);
            if( nwiseRelError > localErrors.MaxNwiseRelError.Value){
              localErrors.MaxNwiseRelError.Set(nwiseRelError,ksup,ib,-1,-1); 
              localErrors.CorrNwiseAbsError.Set(nwiseAbsError,ksup,ib,-1,-1);
            }
            if( nwiseAbsError > localErrors.MaxNwiseAbsError.Value){
              localErrors.MaxNwiseAbsError.Set(nwiseAbsError,ksup,ib,-1,-1); 
            }
          }
        }

        if ( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
          std::vector<LBlock>&  LcolRef = Ref.L( LBj( ksup, grid_ ) );
          Int firstBlk = ( MYROW( grid_ ) == PROW( ksup, grid_ )) ? 1:0;
          for( Int ib = firstBlk; ib < Lcol.size(); ib++ ){
            Real fNormRef = 0.0;
            Real fNorm = 0.0;

            for( Int j = 0; j < Lcol[ib].numCol; j++ ){
              Int gj = FirstBlockCol(ksup,super_) + j;
              Real colNorm = 0.0;
              Real colNormRef = 0.0;

              for( Int i = 0; i < Lcol[ib].numRow; i++ ){
                //Compute the 2-norm of current columns 
                colNorm += abs( pow(Lcol[ib].nzval(i,j),2) );
                colNormRef += abs( pow(LcolRef[ib].nzval(i,j),2) );
              }

              //Compute max row-wise errors
              Real cwiseAbsError = abs(colNorm - colNormRef);
              Real cwiseRelError = cwiseAbsError / abs(colNormRef);
              if( cwiseRelError > localErrors.MaxCwiseRelError.Value){
                localErrors.MaxCwiseRelError.Set(cwiseRelError,ksup,ib,-1,gj);
                localErrors.CorrCwiseAbsError.Set(cwiseAbsError,ksup,ib,-1,gj);
              }
              if( cwiseAbsError > localErrors.MaxCwiseAbsError.Value){
                localErrors.MaxCwiseAbsError.Set(cwiseAbsError,ksup,ib,-1,gj);
              }
            }

            for( Int i = 0; i < Lcol[ib].numRow; i++ ){
              Int gi = FirstBlockRow(ksup,super_) + i;
              Real rowNorm = 0.0;
              Real rowNormRef = 0.0;
              for( Int j = 0; j < Lcol[ib].numCol; j++ ){
                Int gj = FirstBlockCol(ksup,super_) + j;
                //Compute the 2-norm of current row      
                rowNorm += abs( pow(Lcol[ib].nzval(i,j),2) );
                rowNormRef += abs( pow(LcolRef[ib].nzval(i,j),2) );
                std::stringstream msg;
                Real absError = abs(Lcol[ib].nzval(i,j)-LcolRef[ib].nzval(i,j));
                Real relError = absError/abs(LcolRef[ib].nzval(i,j));

                msg<< "["<<ksup<< "] L Block "<<ib<< " Row "<<i<<" Col "<<j<<" is wrong : "<< Lcol[ib].nzval(i,j) << " vs "<<LcolRef[ib].nzval(i,j)<< " relative error is "<< relError <<std::endl; 
                if( relError > localErrors.MaxRelError.Value){
                  localErrors.MaxRelError.Set(relError,ksup,ib,gi,gj); 
                  localErrors.CorrAbsError.Set(absError,ksup,ib,gi,gj); 
                }
                if( absError > localErrors.MaxAbsError.Value){
                  localErrors.MaxAbsError.Set(absError,ksup,ib,gi,gj); 
                }
                if(relError > SANITY_PRECISION){
                  statusOFS<<msg;
                }
              }
              //Compute max row-wise errors
              Real rwiseAbsError = abs(rowNorm - rowNormRef);
              Real rwiseRelError = rwiseAbsError / abs(rowNormRef);
              if( rwiseRelError > localErrors.MaxRwiseRelError.Value){
                localErrors.MaxRwiseRelError.Set(rwiseRelError,ksup,ib,gi,-1);
                localErrors.CorrRwiseAbsError.Set(rwiseAbsError,ksup,ib,gi,-1); 
              }
              if( rwiseAbsError > localErrors.MaxRwiseAbsError.Value){
                localErrors.MaxRwiseAbsError.Set(rwiseAbsError,ksup,ib,gi,-1); 
              }
              //Compute the Frobenius norm
              fNorm += rowNorm;
              fNormRef += rowNormRef;
            }         
            //Compute max norm-wise errors
            Real nwiseAbsError = abs(fNorm - fNormRef);
            Real nwiseRelError = nwiseAbsError / abs(fNormRef);
            if( nwiseRelError > localErrors.MaxNwiseRelError.Value){
              localErrors.MaxNwiseRelError.Set(nwiseRelError,ksup,ib,-1,-1); 
              localErrors.CorrNwiseAbsError.Set(nwiseAbsError,ksup,ib,-1,-1);
            }
            if( nwiseAbsError > localErrors.MaxNwiseAbsError.Value){
              localErrors.MaxNwiseAbsError.Set(nwiseAbsError,ksup,ib,-1,-1); 
            }

          }
        }
      }

      struct{ double val; int rank; } lmaxRelError,gmaxRelError;

      //Get max element-wise errors
      lmaxRelError.val = localErrors.MaxRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrAbsError = localErrors.CorrAbsError; errors.MaxRelError = localErrors.MaxRelError; }
      MPI_Bcast(&errors.CorrAbsError, sizeof(errors.CorrAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxRelError, sizeof(errors.MaxRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxAbsError = localErrors.MaxAbsError; }
      MPI_Bcast(&errors.MaxAbsError, sizeof(errors.MaxAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );


      //Get max norm-wise errors
      lmaxRelError.val = localErrors.MaxNwiseRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrNwiseAbsError = localErrors.CorrNwiseAbsError; errors.MaxNwiseRelError = localErrors.MaxNwiseRelError; }
      MPI_Bcast(&errors.CorrNwiseAbsError, sizeof(errors.CorrNwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxNwiseRelError, sizeof(errors.MaxNwiseRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxNwiseAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxNwiseAbsError = localErrors.MaxNwiseAbsError; }
      MPI_Bcast(&errors.MaxNwiseAbsError, sizeof(errors.MaxNwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      //Get max row-wise errors
      lmaxRelError.val = localErrors.MaxRwiseRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrRwiseAbsError = localErrors.CorrRwiseAbsError; errors.MaxRwiseRelError = localErrors.MaxRwiseRelError; }
      MPI_Bcast(&errors.CorrRwiseAbsError, sizeof(errors.CorrRwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxRwiseRelError, sizeof(errors.MaxRwiseRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxRwiseAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxRwiseAbsError = localErrors.MaxRwiseAbsError; }
      MPI_Bcast(&errors.MaxRwiseAbsError, sizeof(errors.MaxRwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      //Get max column-wise errors
      lmaxRelError.val = localErrors.MaxCwiseRelError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.CorrCwiseAbsError = localErrors.CorrCwiseAbsError; errors.MaxCwiseRelError = localErrors.MaxCwiseRelError; }
      MPI_Bcast(&errors.CorrCwiseAbsError, sizeof(errors.CorrCwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );
      MPI_Bcast(&errors.MaxCwiseRelError, sizeof(errors.MaxCwiseRelError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );

      lmaxRelError.val = localErrors.MaxCwiseAbsError.Value;
      lmaxRelError.rank = MYPROC(grid_);
      MPI_Allreduce( &lmaxRelError, &gmaxRelError, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid_->comm );
      if(MYPROC(grid_) == gmaxRelError.rank){ errors.MaxCwiseAbsError = localErrors.MaxCwiseAbsError; }
      MPI_Bcast(&errors.MaxCwiseAbsError, sizeof(errors.MaxCwiseAbsError) , MPI_BYTE, gmaxRelError.rank, grid_->comm );


#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::CompareDiagonal  ----- 

#endif


    void PMatrix::PMatrixToDistSparseMatrix	( DistSparseMatrix<Scalar>& A )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::PMatrixToDistSparseMatrix");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Converting PMatrix to DistSparseMatrix." << std::endl;
#endif
      Int mpirank = grid_->mpirank;
      Int mpisize = grid_->mpisize;

      std::vector<Int>     rowSend( mpisize );
      std::vector<Int>     colSend( mpisize );
      std::vector<Scalar>  valSend( mpisize );
      std::vector<Int>     sizeSend( mpisize, 0 );
      std::vector<Int>     displsSend( mpisize, 0 );

      std::vector<Int>     rowRecv( mpisize );
      std::vector<Int>     colRecv( mpisize );
      std::vector<Scalar>  valRecv( mpisize );
      std::vector<Int>     sizeRecv( mpisize, 0 );
      std::vector<Int>     displsRecv( mpisize, 0 );

      Int numSuper = this->NumSuper();
      const IntNumVec& permInv = super_->permInv;

      // The number of local columns in DistSparseMatrix format for the
      // processor with rank 0.  This number is the same for processors
      // with rank ranging from 0 to mpisize - 2, and may or may not differ
      // from the number of local columns for processor with rank mpisize -
      // 1.
      Int numColFirst = this->NumCol() / mpisize;

      // Count the size first.
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // L blocks
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            for( Int j = 0; j < Lcol[ib].numCol; j++ ){
              Int jcol = permInv( j + FirstBlockCol( ksup, super_ ) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              sizeSend[dest] += Lcol[ib].numRow;
            }
          }
        } // I own the column of ksup 

        // U blocks
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            IntNumVec& cols = Urow[jb].cols;
            for( Int j = 0; j < cols.m(); j++ ){
              Int jcol = permInv( cols(j) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              sizeSend[dest] += Urow[jb].numRow;
            }
          }
        } // I own the row of ksup
      } // for (ksup)

      // All-to-all exchange of size information
      MPI_Alltoall( 
          &sizeSend[0], 1, MPI_INT,
          &sizeRecv[0], 1, MPI_INT, grid_->comm );



      // Reserve the space
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( ip == 0 ){
          displsSend[ip] = 0;
        }
        else{
          displsSend[ip] = displsSend[ip-1] + sizeSend[ip-1];
        }

        if( ip == 0 ){
          displsRecv[ip] = 0;
        }
        else{
          displsRecv[ip] = displsRecv[ip-1] + sizeRecv[ip-1];
        }
      }
      Int sizeSendTotal = displsSend[mpisize-1] + sizeSend[mpisize-1];
      Int sizeRecvTotal = displsRecv[mpisize-1] + sizeRecv[mpisize-1];

      rowSend.resize( sizeSendTotal );
      colSend.resize( sizeSendTotal );
      valSend.resize( sizeSendTotal );

      rowRecv.resize( sizeRecvTotal );
      colRecv.resize( sizeRecvTotal );
      valRecv.resize( sizeRecvTotal );

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "displsSend = " << displsSend << std::endl;
      statusOFS << "displsRecv = " << displsRecv << std::endl;
#endif

      // Put (row, col, val) to the sending buffer
      std::vector<Int>   cntSize( mpisize, 0 );

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // L blocks
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            IntNumVec&  rows = Lcol[ib].rows;
            NumMat<Scalar>& nzval = Lcol[ib].nzval;
            for( Int j = 0; j < Lcol[ib].numCol; j++ ){
              Int jcol = permInv( j + FirstBlockCol( ksup, super_ ) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              for( Int i = 0; i < rows.m(); i++ ){
                rowSend[displsSend[dest] + cntSize[dest]] = permInv( rows(i) );
                colSend[displsSend[dest] + cntSize[dest]] = jcol;
                valSend[displsSend[dest] + cntSize[dest]] = nzval( i, j );
                cntSize[dest]++;
              }
            }
          }
        } // I own the column of ksup 

        // U blocks
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            IntNumVec& cols = Urow[jb].cols;
            NumMat<Scalar>& nzval = Urow[jb].nzval;
            for( Int j = 0; j < cols.m(); j++ ){
              Int jcol = permInv( cols(j) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              for( Int i = 0; i < Urow[jb].numRow; i++ ){
                rowSend[displsSend[dest] + cntSize[dest]] = 
                  permInv( i + FirstBlockCol( ksup, super_ ) );
                colSend[displsSend[dest] + cntSize[dest]] = jcol;
                valSend[displsSend[dest] + cntSize[dest]] = nzval( i, j );
                cntSize[dest]++;
              }
            }
          }
        } // I own the row of ksup
      }

      // Check sizes match
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( cntSize[ip] != sizeSend[ip] )
          throw std::runtime_error( "Sizes of the sending information do not match." );
      }


      // Alltoallv to exchange information
      mpi::Alltoallv( 
          &rowSend[0], &sizeSend[0], &displsSend[0],
          &rowRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );
      mpi::Alltoallv( 
          &colSend[0], &sizeSend[0], &displsSend[0],
          &colRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );
      mpi::Alltoallv( 
          &valSend[0], &sizeSend[0], &displsSend[0],
          &valRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Alltoallv communication finished." << std::endl;
#endif

      //#if ( _DEBUGlevel_ >= 1 )
      //	for( Int ip = 0; ip < mpisize; ip++ ){
      //		statusOFS << "rowSend[" << ip << "] = " << rowSend[ip] << std::endl;
      //		statusOFS << "rowRecv[" << ip << "] = " << rowRecv[ip] << std::endl;
      //		statusOFS << "colSend[" << ip << "] = " << colSend[ip] << std::endl;
      //		statusOFS << "colRecv[" << ip << "] = " << colRecv[ip] << std::endl;
      //		statusOFS << "valSend[" << ip << "] = " << valSend[ip] << std::endl;
      //		statusOFS << "valRecv[" << ip << "] = " << valRecv[ip] << std::endl;
      //	}
      //#endif

      // Organize the received message.
      Int firstCol = mpirank * numColFirst;
      Int numColLocal;
      if( mpirank == mpisize-1 )
        numColLocal = this->NumCol() - numColFirst * (mpisize-1);
      else
        numColLocal = numColFirst;

      std::vector<std::vector<Int> > rows( numColLocal );
      std::vector<std::vector<Scalar> > vals( numColLocal );

      for( Int ip = 0; ip < mpisize; ip++ ){
        Int*     rowRecvCur = &rowRecv[displsRecv[ip]];
        Int*     colRecvCur = &colRecv[displsRecv[ip]];
        Scalar*  valRecvCur = &valRecv[displsRecv[ip]];
        for( Int i = 0; i < sizeRecv[ip]; i++ ){
          rows[colRecvCur[i]-firstCol].push_back( rowRecvCur[i] );
          vals[colRecvCur[i]-firstCol].push_back( valRecvCur[i] );
        } // for (i)
      } // for (ip)

      // Sort the rows
      std::vector<std::vector<Int> > sortIndex( numColLocal );
      for( Int j = 0; j < numColLocal; j++ ){
        sortIndex[j].resize( rows[j].size() );
        for( Int i = 0; i < sortIndex[j].size(); i++ )
          sortIndex[j][i] = i;
        std::sort( sortIndex[j].begin(), sortIndex[j].end(),
            IndexComp<std::vector<Int>& > ( rows[j] ) );
      } // for (j)

      // Form DistSparseMatrix according to the received message	
      // NOTE: for indicies,  DistSparseMatrix follows the FORTRAN
      // convention (1 based) while PMatrix follows the C convention (0
      // based)
      A.size = this->NumCol();
      A.nnzLocal  = 0;
      A.colptrLocal.Resize( numColLocal + 1 );
      // Note that 1 is important since the index follows the FORTRAN convention
      A.colptrLocal(0) = 1;
      for( Int j = 0; j < numColLocal; j++ ){
        A.nnzLocal += rows[j].size();
        A.colptrLocal(j+1) = A.colptrLocal(j) + rows[j].size();
      }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "nnzLocal = " << A.nnzLocal << std::endl;
      statusOFS << "nnz      = " << A.Nnz()      << std::endl;
#endif


      A.rowindLocal.Resize( A.nnzLocal );
      A.nzvalLocal.Resize(  A.nnzLocal );
      A.comm = grid_->comm;

      Int*     rowPtr = A.rowindLocal.Data();
      Scalar*  nzvalPtr = A.nzvalLocal.Data();
      for( Int j = 0; j < numColLocal; j++ ){
        std::vector<Int>& rowsCur = rows[j];
        std::vector<Int>& sortIndexCur = sortIndex[j];
        std::vector<Scalar>& valsCur = vals[j];
        for( Int i = 0; i < rows[j].size(); i++ ){
          // Note that 1 is important since the index follows the FORTRAN convention
          *(rowPtr++)   = rowsCur[sortIndexCur[i]] + 1;
          *(nzvalPtr++) = valsCur[sortIndexCur[i]]; 
        }
      }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "A.colptrLocal[end]   = " << A.colptrLocal(numColLocal) << std::endl;
      statusOFS << "A.rowindLocal.size() = " << A.rowindLocal.m() << std::endl;
      statusOFS << "A.rowindLocal[end]   = " << A.rowindLocal(A.nnzLocal-1) << std::endl;
      statusOFS << "A.nzvalLocal[end]    = " << A.nzvalLocal(A.nnzLocal-1) << std::endl;
#endif


#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::PMatrixToDistSparseMatrix  ----- 



    void PMatrix::PMatrixToDistSparseMatrix	( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B	)
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::PMatrixToDistSparseMatrix");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Converting PMatrix to DistSparseMatrix (2nd format)." << std::endl;
#endif
      Int mpirank = grid_->mpirank;
      Int mpisize = grid_->mpisize;

      std::vector<Int>     rowSend( mpisize );
      std::vector<Int>     colSend( mpisize );
      std::vector<Scalar>  valSend( mpisize );
      std::vector<Int>     sizeSend( mpisize, 0 );
      std::vector<Int>     displsSend( mpisize, 0 );

      std::vector<Int>     rowRecv( mpisize );
      std::vector<Int>     colRecv( mpisize );
      std::vector<Scalar>  valRecv( mpisize );
      std::vector<Int>     sizeRecv( mpisize, 0 );
      std::vector<Int>     displsRecv( mpisize, 0 );

      Int numSuper = this->NumSuper();
      const IntNumVec& permInv = super_->permInv;

      // The number of local columns in DistSparseMatrix format for the
      // processor with rank 0.  This number is the same for processors
      // with rank ranging from 0 to mpisize - 2, and may or may not differ
      // from the number of local columns for processor with rank mpisize -
      // 1.
      Int numColFirst = this->NumCol() / mpisize;

      // Count the size first.
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // L blocks
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            for( Int j = 0; j < Lcol[ib].numCol; j++ ){
              Int jcol = permInv( j + FirstBlockCol( ksup, super_ ) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              sizeSend[dest] += Lcol[ib].numRow;
            }
          }
        } // I own the column of ksup 

        // U blocks
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            IntNumVec& cols = Urow[jb].cols;
            for( Int j = 0; j < cols.m(); j++ ){
              Int jcol = permInv( cols(j) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              sizeSend[dest] += Urow[jb].numRow;
            }
          }
        } // I own the row of ksup
      } // for (ksup)

      // All-to-all exchange of size information
      MPI_Alltoall( 
          &sizeSend[0], 1, MPI_INT,
          &sizeRecv[0], 1, MPI_INT, grid_->comm );



      // Reserve the space
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( ip == 0 ){
          displsSend[ip] = 0;
        }
        else{
          displsSend[ip] = displsSend[ip-1] + sizeSend[ip-1];
        }

        if( ip == 0 ){
          displsRecv[ip] = 0;
        }
        else{
          displsRecv[ip] = displsRecv[ip-1] + sizeRecv[ip-1];
        }
      }
      Int sizeSendTotal = displsSend[mpisize-1] + sizeSend[mpisize-1];
      Int sizeRecvTotal = displsRecv[mpisize-1] + sizeRecv[mpisize-1];

      rowSend.resize( sizeSendTotal );
      colSend.resize( sizeSendTotal );
      valSend.resize( sizeSendTotal );

      rowRecv.resize( sizeRecvTotal );
      colRecv.resize( sizeRecvTotal );
      valRecv.resize( sizeRecvTotal );

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "displsSend = " << displsSend << std::endl;
      statusOFS << "displsRecv = " << displsRecv << std::endl;
#endif

      // Put (row, col, val) to the sending buffer
      std::vector<Int>   cntSize( mpisize, 0 );


      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // L blocks
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            IntNumVec&  rows = Lcol[ib].rows;
            NumMat<Scalar>& nzval = Lcol[ib].nzval;
            for( Int j = 0; j < Lcol[ib].numCol; j++ ){
              Int jcol = permInv( j + FirstBlockCol( ksup, super_ ) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              for( Int i = 0; i < rows.m(); i++ ){
                rowSend[displsSend[dest] + cntSize[dest]] = permInv( rows(i) );
                colSend[displsSend[dest] + cntSize[dest]] = jcol;
                valSend[displsSend[dest] + cntSize[dest]] = nzval( i, j );
                cntSize[dest]++;
              }
            }
          }
        } // I own the column of ksup 


        // U blocks
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            IntNumVec& cols = Urow[jb].cols;
            NumMat<Scalar>& nzval = Urow[jb].nzval;
            for( Int j = 0; j < cols.m(); j++ ){
              Int jcol = permInv( cols(j) );
              Int dest = std::min( jcol / numColFirst, mpisize - 1 );
              for( Int i = 0; i < Urow[jb].numRow; i++ ){
                rowSend[displsSend[dest] + cntSize[dest]] = 
                  permInv( i + FirstBlockCol( ksup, super_ ) );
                colSend[displsSend[dest] + cntSize[dest]] = jcol;
                valSend[displsSend[dest] + cntSize[dest]] = nzval( i, j );
                cntSize[dest]++;
              }
            }
          }
        } // I own the row of ksup
      }



      // Check sizes match
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( cntSize[ip] != sizeSend[ip] )
          throw std::runtime_error( "Sizes of the sending information do not match." );
      }

      // Alltoallv to exchange information
      mpi::Alltoallv( 
          &rowSend[0], &sizeSend[0], &displsSend[0],
          &rowRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );
      mpi::Alltoallv( 
          &colSend[0], &sizeSend[0], &displsSend[0],
          &colRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );
      mpi::Alltoallv( 
          &valSend[0], &sizeSend[0], &displsSend[0],
          &valRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Alltoallv communication finished." << std::endl;
#endif

      //#if ( _DEBUGlevel_ >= 1 )
      //	for( Int ip = 0; ip < mpisize; ip++ ){
      //		statusOFS << "rowSend[" << ip << "] = " << rowSend[ip] << std::endl;
      //		statusOFS << "rowRecv[" << ip << "] = " << rowRecv[ip] << std::endl;
      //		statusOFS << "colSend[" << ip << "] = " << colSend[ip] << std::endl;
      //		statusOFS << "colRecv[" << ip << "] = " << colRecv[ip] << std::endl;
      //		statusOFS << "valSend[" << ip << "] = " << valSend[ip] << std::endl;
      //		statusOFS << "valRecv[" << ip << "] = " << valRecv[ip] << std::endl;
      //	}
      //#endif

      // Organize the received message.
      Int firstCol = mpirank * numColFirst;
      Int numColLocal;
      if( mpirank == mpisize-1 )
        numColLocal = this->NumCol() - numColFirst * (mpisize-1);
      else
        numColLocal = numColFirst;

      std::vector<std::vector<Int> > rows( numColLocal );
      std::vector<std::vector<Scalar> > vals( numColLocal );

      for( Int ip = 0; ip < mpisize; ip++ ){
        Int*     rowRecvCur = &rowRecv[displsRecv[ip]];
        Int*     colRecvCur = &colRecv[displsRecv[ip]];
        Scalar*  valRecvCur = &valRecv[displsRecv[ip]];
        for( Int i = 0; i < sizeRecv[ip]; i++ ){
          rows[colRecvCur[i]-firstCol].push_back( rowRecvCur[i] );
          vals[colRecvCur[i]-firstCol].push_back( valRecvCur[i] );
        } // for (i)
      } // for (ip)

      // Sort the rows
      std::vector<std::vector<Int> > sortIndex( numColLocal );
      for( Int j = 0; j < numColLocal; j++ ){
        sortIndex[j].resize( rows[j].size() );
        for( Int i = 0; i < sortIndex[j].size(); i++ )
          sortIndex[j][i] = i;
        std::sort( sortIndex[j].begin(), sortIndex[j].end(),
            IndexComp<std::vector<Int>& > ( rows[j] ) );
      } // for (j)

      // Form DistSparseMatrix according to the received message	
      // NOTE: for indicies,  DistSparseMatrix follows the FORTRAN
      // convention (1 based) while PMatrix follows the C convention (0
      // based)
      if( A.size != this->NumCol() ){
        throw std::runtime_error( "The DistSparseMatrix providing the pattern has a different size from PMatrix." );
      }
      if( A.colptrLocal.m() != numColLocal + 1 ){
        throw std::runtime_error( "The DistSparseMatrix providing the pattern has a different number of local columns from PMatrix." );
      }

      B.size = A.size;
      B.nnz  = A.nnz;
      B.nnzLocal = A.nnzLocal;
      B.colptrLocal = A.colptrLocal;
      B.rowindLocal = A.rowindLocal;
      B.nzvalLocal.Resize( B.nnzLocal );
      SetValue( B.nzvalLocal, SCALAR_ZERO );
      // Make sure that the communicator of A and B are the same.
      // FIXME Find a better way to compare the communicators
      //			if( grid_->comm != A.comm ){
      //				throw std::runtime_error( "The DistSparseMatrix providing the pattern has a different communicator from PMatrix." );
      //			}
      B.comm = grid_->comm;

      Int*     rowPtr = B.rowindLocal.Data();
      Scalar*  nzvalPtr = B.nzvalLocal.Data();
      for( Int j = 0; j < numColLocal; j++ ){
        std::vector<Int>& rowsCur = rows[j];
        std::vector<Int>& sortIndexCur = sortIndex[j];
        std::vector<Scalar>& valsCur = vals[j];
        std::vector<Int>  rowsCurSorted( rowsCur.size() );
        // Note that 1 is important since the index follows the FORTRAN convention
        for( Int i = 0; i < rowsCurSorted.size(); i++ ){
          rowsCurSorted[i] = rowsCur[sortIndexCur[i]] + 1;
        }

        // Search and match the indices
        std::vector<Int>::iterator it;
        for( Int i = B.colptrLocal(j) - 1; 
            i < B.colptrLocal(j+1) - 1; i++ ){
          it = std::lower_bound( rowsCurSorted.begin(), rowsCurSorted.end(),
              *(rowPtr++) );
          if( it == rowsCurSorted.end() ){
            // Did not find the row, set it to zero
            *(nzvalPtr++) = SCALAR_ZERO;
          }
          else{
            // Found the row, set it according to the received value
            *(nzvalPtr++) = valsCur[ sortIndexCur[it-rowsCurSorted.begin()] ];
          }
        } // for (i)	
      } // for (j)

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "B.colptrLocal[end]   = " << B.colptrLocal(numColLocal) << std::endl;
      statusOFS << "B.rowindLocal.size() = " << B.rowindLocal.m() << std::endl;
      statusOFS << "B.rowindLocal[end]   = " << B.rowindLocal(B.nnzLocal-1) << std::endl;
      statusOFS << "B.nzvalLocal[end]    = " << B.nzvalLocal(B.nnzLocal-1) << std::endl;
#endif


#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::PMatrixToDistSparseMatrix  ----- 


    // A (maybe) more memory efficient way for converting the PMatrix to a
    // DistSparseMatrix structure.
    //
    // FIXME NOTE: This routine assumes the matrix to be symmetric!
    void PMatrix::PMatrixToDistSparseMatrix2 ( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::PMatrixToDistSparseMatrix2");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Converting PMatrix to DistSparseMatrix (2nd format)." << std::endl;
#endif
      Int mpirank = grid_->mpirank;
      Int mpisize = grid_->mpisize;

      std::vector<Int>     rowSend( mpisize );
      std::vector<Int>     colSend( mpisize );
      std::vector<Scalar>  valSend( mpisize );
      std::vector<Int>     sizeSend( mpisize, 0 );
      std::vector<Int>     displsSend( mpisize, 0 );

      std::vector<Int>     rowRecv( mpisize );
      std::vector<Int>     colRecv( mpisize );
      std::vector<Scalar>  valRecv( mpisize );
      std::vector<Int>     sizeRecv( mpisize, 0 );
      std::vector<Int>     displsRecv( mpisize, 0 );

      Int numSuper = this->NumSuper();
      const IntNumVec& perm    = super_->perm;
      const IntNumVec& permInv = super_->permInv;


      // Count the sizes from the A matrix first
      Int numColFirst = this->NumCol() / mpisize;
      Int firstCol = mpirank * numColFirst;
      Int numColLocal;
      if( mpirank == mpisize-1 )
        numColLocal = this->NumCol() - numColFirst * (mpisize-1);
      else
        numColLocal = numColFirst;

      Int*     rowPtr = A.rowindLocal.Data();
      Int*     colPtr = A.colptrLocal.Data();

      for( Int j = 0; j < numColLocal; j++ ){
        Int col         = perm( firstCol + j );
        Int blockColIdx = BlockIdx( col, super_ );
        Int procCol     = PCOL( blockColIdx, grid_ );
        for( Int i = colPtr[j] - 1; i < colPtr[j+1] - 1; i++ ){
          Int row         = perm( *(rowPtr++) - 1 );
          Int blockRowIdx = BlockIdx( row, super_ );
          Int procRow     = PROW( blockRowIdx, grid_ );
          Int dest = PNUM( procRow, procCol, grid_ );
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "BlockIdx = " << blockRowIdx << ", " <<blockColIdx << std::endl;
          statusOFS << procRow << ", " << procCol << ", " 
            << dest << std::endl;
#endif
          sizeSend[dest]++;
        } // for (i)
      } // for (j)

      // All-to-all exchange of size information
      MPI_Alltoall( 
          &sizeSend[0], 1, MPI_INT,
          &sizeRecv[0], 1, MPI_INT, grid_->comm );

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << std::endl << "sizeSend: " << sizeSend << std::endl;
      statusOFS << std::endl << "sizeRecv: " << sizeRecv << std::endl;
#endif



      // Reserve the space
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( ip == 0 ){
          displsSend[ip] = 0;
        }
        else{
          displsSend[ip] = displsSend[ip-1] + sizeSend[ip-1];
        }

        if( ip == 0 ){
          displsRecv[ip] = 0;
        }
        else{
          displsRecv[ip] = displsRecv[ip-1] + sizeRecv[ip-1];
        }
      }

      Int sizeSendTotal = displsSend[mpisize-1] + sizeSend[mpisize-1];
      Int sizeRecvTotal = displsRecv[mpisize-1] + sizeRecv[mpisize-1];

      rowSend.resize( sizeSendTotal );
      colSend.resize( sizeSendTotal );
      valSend.resize( sizeSendTotal );

      rowRecv.resize( sizeRecvTotal );
      colRecv.resize( sizeRecvTotal );
      valRecv.resize( sizeRecvTotal );

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "displsSend = " << displsSend << std::endl;
      statusOFS << "displsRecv = " << displsRecv << std::endl;
#endif

      // Put (row, col) to the sending buffer
      std::vector<Int>   cntSize( mpisize, 0 );

      rowPtr = A.rowindLocal.Data();
      colPtr = A.colptrLocal.Data();

      for( Int j = 0; j < numColLocal; j++ ){
        Int col         = perm( firstCol + j );
        Int blockColIdx = BlockIdx( col, super_ );
        Int procCol     = PCOL( blockColIdx, grid_ );
        for( Int i = colPtr[j] - 1; i < colPtr[j+1] - 1; i++ ){
          Int row         = perm( *(rowPtr++) - 1 );
          Int blockRowIdx = BlockIdx( row, super_ );
          Int procRow     = PROW( blockRowIdx, grid_ );
          Int dest = PNUM( procRow, procCol, grid_ );
          rowSend[displsSend[dest] + cntSize[dest]] = row;
          colSend[displsSend[dest] + cntSize[dest]] = col;
          cntSize[dest]++;
        } // for (i)
      } // for (j)


      // Check sizes match
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( cntSize[ip] != sizeSend[ip] )
          throw std::runtime_error( "Sizes of the sending information do not match." );
      }

      // Alltoallv to exchange information
      mpi::Alltoallv( 
          &rowSend[0], &sizeSend[0], &displsSend[0],
          &rowRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );
      mpi::Alltoallv( 
          &colSend[0], &sizeSend[0], &displsSend[0],
          &colRecv[0], &sizeRecv[0], &displsRecv[0],
          grid_->comm );

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "Alltoallv communication of nonzero indices finished." << std::endl;
#endif


#if ( _DEBUGlevel_ >= 1 )
      for( Int ip = 0; ip < mpisize; ip++ ){
        statusOFS << "rowSend[" << ip << "] = " << rowSend[ip] << std::endl;
        statusOFS << "rowRecv[" << ip << "] = " << rowRecv[ip] << std::endl;
        statusOFS << "colSend[" << ip << "] = " << colSend[ip] << std::endl;
        statusOFS << "colRecv[" << ip << "] = " << colRecv[ip] << std::endl;
      }
#endif

      // For each (row, col), fill the nonzero values to valRecv locally.
      for( Int g = 0; g < sizeRecvTotal; g++ ){
        Int row = rowRecv[g];
        Int col = colRecv[g];

        Int blockRowIdx = BlockIdx( row, super_ );
        Int blockColIdx = BlockIdx( col, super_ );

        // Search for the nzval
        bool isFound = false;

        if( blockColIdx <= blockRowIdx ){
          // Data on the L side

          std::vector<LBlock>&  Lcol = this->L( LBj( blockColIdx, grid_ ) );

          for( Int ib = 0; ib < Lcol.size(); ib++ ){
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << "blockRowIdx = " << blockRowIdx << ", Lcol[ib].blockIdx = " << Lcol[ib].blockIdx << ", blockColIdx = " << blockColIdx << std::endl;
#endif
            if( Lcol[ib].blockIdx == blockRowIdx ){
              IntNumVec& rows = Lcol[ib].rows;
              for( int iloc = 0; iloc < Lcol[ib].numRow; iloc++ ){
                if( rows[iloc] == row ){
                  Int jloc = col - FirstBlockCol( blockColIdx, super_ );
                  valRecv[g] = Lcol[ib].nzval( iloc, jloc );
                  isFound = true;
                  break;
                } // found the corresponding row
              }
            }
            if( isFound == true ) break;  
          } // for (ib)
        } 
        else{
          // Data on the U side

          std::vector<UBlock>&  Urow = this->U( LBi( blockRowIdx, grid_ ) );

          for( Int jb = 0; jb < Urow.size(); jb++ ){
            if( Urow[jb].blockIdx == blockColIdx ){
              IntNumVec& cols = Urow[jb].cols;
              for( int jloc = 0; jloc < Urow[jb].numCol; jloc++ ){
                if( cols[jloc] == col ){
                  Int iloc = row - FirstBlockRow( blockRowIdx, super_ );
                  valRecv[g] = Urow[jb].nzval( iloc, jloc );
                  isFound = true;
                  break;
                } // found the corresponding col
              }
            }
            if( isFound == true ) break;  
          } // for (jb)
        } // if( blockColIdx <= blockRowIdx ) 

        // Did not find the corresponding row, set the value to zero.
        if( isFound == false ){
          statusOFS << "In the permutated order, (" << row << ", " << col <<
            ") is not found in PMatrix." << std::endl;
          valRecv[g] = SCALAR_ZERO;
        }

      } // for (g)


      // Feed back valRecv to valSend through Alltoallv. NOTE: for the
      // values, the roles of "send" and "recv" are swapped.
      mpi::Alltoallv( 
          &valRecv[0], &sizeRecv[0], &displsRecv[0],
          &valSend[0], &sizeSend[0], &displsSend[0],
          grid_->comm );

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "Alltoallv communication of nonzero values finished." << std::endl;
#endif

      // Put the nonzero values from valSend to the matrix B.
      B.size = A.size;
      B.nnz  = A.nnz;
      B.nnzLocal = A.nnzLocal;
      B.colptrLocal = A.colptrLocal;
      B.rowindLocal = A.rowindLocal;
      B.nzvalLocal.Resize( B.nnzLocal );
      SetValue( B.nzvalLocal, SCALAR_ZERO );
      // Make sure that the communicator of A and B are the same.
      // FIXME Find a better way to compare the communicators
      //			if( grid_->comm != A.comm ){
      //				throw std::runtime_error( "The DistSparseMatrix providing the pattern has a different communicator from PMatrix." );
      //			}
      B.comm = grid_->comm;

      for( Int i = 0; i < mpisize; i++ )
        cntSize[i] = 0;

      rowPtr = B.rowindLocal.Data();
      colPtr = B.colptrLocal.Data();
      Scalar* valPtr = B.nzvalLocal.Data();

      for( Int j = 0; j < numColLocal; j++ ){
        Int col         = perm( firstCol + j );
        Int blockColIdx = BlockIdx( col, super_ );
        Int procCol     = PCOL( blockColIdx, grid_ );
        for( Int i = colPtr[j] - 1; i < colPtr[j+1] - 1; i++ ){
          Int row         = perm( *(rowPtr++) - 1 );
          Int blockRowIdx = BlockIdx( row, super_ );
          Int procRow     = PROW( blockRowIdx, grid_ );
          Int dest = PNUM( procRow, procCol, grid_ );
          *(valPtr++) = valSend[displsSend[dest] + cntSize[dest]];
          cntSize[dest]++;
        } // for (i)
      } // for (j)

      // Check sizes match
      for( Int ip = 0; ip < mpisize; ip++ ){
        if( cntSize[ip] != sizeSend[ip] )
          throw std::runtime_error( "Sizes of the sending information do not match." );
      }


#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    }     // -----  end of method PMatrix::PMatrixToDistSparseMatrix2  ----- 




    Int PMatrix::NnzLocal	(  )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::NnzLocal");
#endif
      Int numSuper = this->NumSuper();
      Int nnzLocal = 0;
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          std::vector<LBlock>& Lcol = this->L( LBj( ksup, grid_ ) );
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            nnzLocal += Lcol[ib].numRow * Lcol[ib].numCol;
          }
        } // if I own the column of ksup
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
          std::vector<UBlock>& Urow = this->U( LBi( ksup, grid_ ) );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            nnzLocal += Urow[jb].numRow * Urow[jb].numCol;
          }
        } // if I own the row of ksup
      }

#ifndef _RELEASE_
      PopCallStack();
#endif

      return nnzLocal;
    } 		// -----  end of method PMatrix::NnzLocal  ----- 


    LongInt PMatrix::Nnz	(  )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::Nnz");
#endif
      LongInt nnzLocal = LongInt( this->NnzLocal() );
      LongInt nnz;

      MPI_Allreduce( &nnzLocal, &nnz, 1, MPI_LONG_LONG, MPI_SUM, grid_->comm );

#ifndef _RELEASE_
      PopCallStack();
#endif

      return nnz;
    } 		// -----  end of method PMatrix::Nnz  ----- 

    void PMatrix::GetNegativeInertia	( Real& inertia )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::GetNegativeInertia");
#endif
      Int numSuper = this->NumSuper(); 

      Real inertiaLocal = 0.0;
      inertia          = 0.0;

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // I own the diagonal block	
        if( MYROW( grid_ ) == PROW( ksup, grid_ ) &&
            MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
          LBlock& LB = this->L( LBj( ksup, grid_ ) )[0];
          for( Int i = 0; i < LB.numRow; i++ ){
            if( LB.nzval(i, i).real() < 0 )
              inertiaLocal++;
          }
        }
      }

      // All processors own diag
      mpi::Allreduce( &inertiaLocal, &inertia, 1, MPI_SUM, grid_->comm );

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrix::GetNegativeInertia  ----- 


  } // namespace PEXSI

/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin and Lin Lin

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
/// @file ngchol_interf_impl.hpp
/// @brief Implementation of interface with NGCHOL.
/// @date 2014-07-08 Original version
#ifndef _PEXSI_NGCHOL_INTERF_IMPL_HPP_
#define _PEXSI_NGCHOL_INTERF_IMPL_HPP_

// Interface with NGCHOL
#include "ngchol.hpp"
#include "ngchol/sp_blas.hpp"

// Interface with PSelInv
#include "pexsi/pselinv.hpp"

// Interface with sparse matrix (CSC format)
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/environment.hpp"
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/NumMat.hpp"
#include "pexsi/NumVec.hpp"

// Interface with LAPACK
#include "pexsi/lapack.hpp"

namespace PEXSI{

template<typename T> void NGCHOLMatrixToSuperNode( 
    LIBCHOLESKY::SupernodalMatrix<T>& SMat,
    SuperNodeType& super ){
#ifndef _RELEASE_
	PushCallStack("NGCHOLMatrixToSuperNode");
#endif
  Int n = SMat.Size();

  // perm
  const LIBCHOLESKY::IntNumVec& SPerm = SMat.GetColPerm();
  super.perm.Resize( SPerm.m() );
  for( Int i = 0; i < SPerm.m(); i++ ){
    super.perm[i] = SPerm[i];
  }
  
  // permInv
  super.permInv.Resize( n );
  for( Int i = 0; i < n; i++ ){
    super.permInv[i] = i;
  }
  std::sort( super.permInv.Data(), super.permInv.Data() + n,
      IndexComp<IntNumVec&>(super.perm) );

  LIBCHOLESKY::IntNumVec& XSuper = SMat.GetSupernodalPartition();
  Int numSuper = XSuper.m() - 1;

  // superPtr
  IntNumVec& superPtr = super.superPtr;
  superPtr.Resize( numSuper + 1 );
  for( Int i = 0; i < numSuper + 1; i++ ){
    superPtr[i] = XSuper[i] - 1;
  }

  // superIdx
  IntNumVec& superIdx = super.superIdx;
  superIdx.Resize( n );
  const LIBCHOLESKY::IntNumVec& superMember = SMat.GetSupMembership();
  for( Int i = 0; i < n; i++ ){
    superIdx(i) = superMember[i] - 1;
  }

  // etree
  LIBCHOLESKY::ETree& etree = SMat.GetETree();
  LIBCHOLESKY::IntNumVec etreeVec(n);
  for( Int i = 0; i < n; i++ ){
    etreeVec[i] = i+1;
  }
  etreeVec = etree.ToPostOrder( etreeVec );
  super.etree.Resize(n);
  for( Int i = 0; i < n; i++ ){
    super.etree[i] = etreeVec[i]-1;
  }

#ifndef _RELEASE_
	PopCallStack();
#endif
}  // -----  end of NGCHOLMatrixToSuperNode ----- 



template<typename T> void NGCHOLMatrixToPMatrix( 
    LIBCHOLESKY::SupernodalMatrix<T>& SMat,
    PMatrix<T>& PMat ){
#ifndef _RELEASE_
	PushCallStack("NGCHOLMatrixToPMatrix");
#endif
  // This routine assumes that the g, supernode and options of PMatrix
  // has been set outside this routine.
  
  Int mpirank, mpisize;
  const GridType *g = PMat.Grid();

  // FIXME Check PMatrix and SupernodalMatrix has the same communicator
  MPI_Comm comm = g->comm;
  MPI_Comm_size(comm, &mpisize); 
  MPI_Comm_rank(comm, &mpirank);

  Int nprow = g->numProcRow, npcol = g->numProcCol;

  Int n = SMat.Size();
  PMat.ColBlockIdx().clear();
  PMat.RowBlockIdx().clear();
  PMat.ColBlockIdx().resize( PMat.NumLocalBlockCol() );
  PMat.RowBlockIdx().resize( PMat.NumLocalBlockRow() );

  const IntNumVec& superIdx = PMat.SuperNode()->superIdx;
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "superIdx = " << superIdx << std::endl;
#endif


  // for loop over all supernodes
  //
  //   if the current processor owns the supernode (the ownership is block
  //   cyclic)
  //     serialize the information
  //     broadcast the information to all processors
  //   elseif 
  //     receive the information from the processor owning the
  //     supernode, and save the information in a deserialized buffer
  //   endif
  //
  //   (Local operation from now)
  //   if the current processor owns the correct processor column
  //     for loop over the local blocks
  //       for loop over each row
  //         if the current row number belong to the current processor
  //           add the row index to a vector for LBuf
  //         endif
  //       endfor
  //     endfor
  //
  //     Allocate the sign of LBuf for saving the nonzero values
  //     Convert the format of rowind
  //
  //     for loop over the local blocks
  //       for loop over each row
  //         if the current row number belong to the current processor
  //           Append the nonzero values to nzval
  //         endif
  //       endfor
  //     endfor
  //   endif
  //
  //   discard the temporary information in the buffer for all
  //   processors
  //
  // endfor
  //
  // Perform SendToCrossDiagonal to fill the U part.
  //
  // Remark:
  //   1.  The communiation is performed in a supernode-by-supernode
  //   way.  This may not be very fast.  A first improvement is to
  //   broadcast all local supernodes at once and then do local
  //   post-processing.

  Int numSuper = PMat.NumSuper();
  LIBCHOLESKY::Icomm snodeIcomm;
  std::vector<char> buffer;
  for( Int iSuper = 0; iSuper < numSuper; iSuper++ ){
    LIBCHOLESKY::SuperNode<T> snode;

    if( mpirank == ( iSuper % mpisize ) ){
      // Get the local supernode
      Int iSuperLocal = iSuper / mpisize;
      LIBCHOLESKY::SuperNode<T>& snodetmp = 
        SMat.GetLocalSupernode(iSuperLocal);
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "iSuper = " << iSuper << ", iSuperLocal = " <<
        iSuperLocal << ", id = " << snodetmp.Id() << ", size = " << 
        snodetmp.Size() << ", #Block = " << snodetmp.NZBlockCnt() <<
        std::endl;
    statusOFS << "snode (before bcast) = " << snodetmp << std::endl;
#endif
      // Serialize the information in the current supernode
//      std::stringstream sstm;
//      serialize( snodetmp.Id(), sstm, NO_MASK );
//      serialize( snodetmp.Size(), sstm, NO_MASK );
//      serialize( snodetmp.NZBlockCnt(), sstm, NO_MASK );
//      for( Int blkidx = 0; blkidx < snodetmp.NZBlockCnt(); blkidx++){
//        LIBCHOLESKY::NZBlockDesc & nzblk_desc =
//          snodetmp.GetNZBlockDesc( blkidx );
//      } // for (blkidx)

      LIBCHOLESKY::Serialize( snodeIcomm, snodetmp );
      Int msgSize = snodeIcomm.size();
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "msgSize = " << msgSize << std::endl;
#endif

      // Communicate the supernode
      MPI_Bcast( &msgSize, 1, MPI_INT, mpirank, comm );
      MPI_Bcast( snodeIcomm.front(), msgSize, MPI_CHAR, mpirank, comm );
      // Copy the data from buffer to snode
      LIBCHOLESKY::Deserialize( snodeIcomm.front(), snode );
    } // if owning the supernode
    else{
      // Receive the supernode
      Int rootRank = ( iSuper % mpisize );
      Int msgSize;
      MPI_Bcast( &msgSize, 1, MPI_INT, rootRank, comm );
      buffer.resize(msgSize);
      MPI_Bcast( &buffer[0], msgSize, MPI_CHAR, rootRank, comm );
      LIBCHOLESKY::Deserialize( &buffer[0], snode );
    } // if not owning the supernode but is in the receiving column

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "All communication is finished." << std::endl;
#endif

    // Local operation from now
    if( ( mpirank % npcol ) == ( iSuper % npcol ) ){
      Int jb = iSuper / npcol;
      std::vector<LBlock<T> >& Lcol = PMat.L(jb);
      std::set<Int> blkSet;
      Int superSize = snode.Size();

      // Count the number of blocks in the supernode belonging to this
      // processor
      for( Int blkidx = 0; blkidx < snode.NZBlockCnt(); blkidx++ ){
        LIBCHOLESKY::NZBlockDesc desc = snode.GetNZBlockDesc( blkidx );
        Int nrows = snode.NRows(blkidx);
        Int firstRow = desc.GIndex - 1;
        Int lastRow = firstRow + nrows - 1;
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "firstRow = " << firstRow << ", lastRow = " << lastRow << std::endl;
#endif
        for( Int i = superIdx(firstRow); i <= superIdx(lastRow); i++ ){
          if( ( mpirank % nprow ) == ( i % nprow ) ){
            blkSet.insert( i );
          } // if the current processor is in the right processor row
        }
      } // for ( blkidx )
      
      Int numBlkLocal = blkSet.size();
      std::vector<Int> blkVec;
      blkVec.insert( blkVec.end(), blkSet.begin(), blkSet.end() );
      Lcol.resize( blkVec.size() );

      // Allocate the nonzero rows and nzvals 
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Lcol.size = " << Lcol.size() << std::endl;
      statusOFS << "blkSet.size = " << blkSet.size() << std::endl;
      statusOFS << "blkVec = " << blkVec << std::endl;
#endif

      std::vector<std::vector<Int> > rowsBlk( Lcol.size() );
      std::vector<std::vector<T> > nzvalBlk( Lcol.size() );

      for( Int blkidx = 0; blkidx < snode.NZBlockCnt(); blkidx++ ){
        LIBCHOLESKY::NZBlockDesc desc = snode.GetNZBlockDesc( blkidx );
        Int nrows = snode.NRows(blkidx);
        Int firstRow = desc.GIndex - 1;
        Int lastRow = firstRow + nrows - 1;
        T* nzval = snode.GetNZval( desc.Offset );
        std::vector<Int>::iterator vi;
        Int pos;
        for( Int i = firstRow; i <= lastRow; i++ ){
          vi = lower_bound( blkVec.begin(), blkVec.end(), superIdx(i) );
          pos = vi - blkVec.begin();
          rowsBlk[pos].push_back(i);
          nzvalBlk[pos].insert(
              nzvalBlk[pos].end(),
              nzval, nzval + superSize ); 
          nzval += superSize;
        }
      } // for ( blkidx )
      
      // Save the information to Lcol
      for ( Int iblk = 0; iblk < Lcol.size(); iblk++ ){
        std::vector<Int>& rows = rowsBlk[iblk];
        std::vector<T>& nzval = nzvalBlk[iblk];
        LBlock<T>& LB = Lcol[iblk];
        LB.blockIdx = GBi( blkVec[iblk], g );
        LB.numRow = rows.size();
        LB.numCol = superSize;
        LB.rows = IntNumVec( LB.numRow, true, &rows[0] );
        if( LB.numRow * LB.numCol != nzval.size() ){
          std::ostringstream msg;
          msg << "message size does not match for the blockIdx " << LB.blockIdx << std::endl
            << "LB.numRow * LB.numCol = " << LB.numRow * LB.numCol << std::endl
            << "nzval.size            = " << nzval.size() << std::endl;
          throw std::runtime_error( msg.str().c_str() );
        }
        // Convert the row major format to column major format
        Transpose( NumMat<T>( LB.numCol, LB.numRow, true,
            &nzval[0] ), LB.nzval );
      }
    } // if the current processor is in the right processor column

    // Set the MPI Barrier
    MPI_Barrier( comm );
  }




#ifndef _RELEASE_
	PopCallStack();
#endif

}  // -----  end of NGCHOLMatrixToPMatrix ----- 


}


#endif //_PEXSI_NGCHOL_INTERF_IMPL_HPP_


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
/// @file pselinv_asym.hpp
/// @brief Main file for parallel selected inversion on asymmetric matrices.
/// @date 2013-08-05
#ifndef _PEXSI_PSELINV_ASYM_HPP_
#define _PEXSI_PSELINV_ASYM_HPP_

// *********************************************************************
//  Common utilities
// *********************************************************************

#include  "environment.hpp"
#include	"NumVec.hpp"
#include	"NumMat.hpp" 
#include  "sparse_matrix.hpp"

#include  "superlu_dist_interf.hpp"
#include  "mpi_interf.hpp"
#include	"utility.hpp"
#include	"blas.hpp"
#include	"lapack.hpp"
#include <set>

#include	"pselinv.hpp"


namespace PEXSI{




  /**********************************************************************
   * Main data structure in PSelInv: PMatrix
   **********************************************************************/

  /// @class PMatrixAsym
  ///
  /// @brief PMatrixAsym contains the main data structure and the
  /// computational routine for the parallel selected inversion.  
  /// 
  /// **NOTE** The following is a bit obsolete.
  ///
  /// Procedure for Selected Inversion 
  /// --------------------------------
  ///
  /// After factorizing a SuperLUMatrix luMat (See SuperLUMatrix for
  /// information on how to perform factorization), perform the following
  /// steps for parallel selected inversion.
  /// 
  /// - Conversion from SuperLU_DIST.
  ///   
  ///   Symbolic information
  ///
  ///       SuperNodeType super; 
  ///       PMatrix PMloc;
  ///       luMat.SymbolicToSuperNode( super );  
  ///   
  ///   Numerical information, both L and U.
  ///
  ///       luMat.LUstructToPMatrix( PMloc ); 
  ///
  /// - Preparation.
  ///
  ///   Construct the communication pattern for SelInv.
  ///
  ///       PMloc.ConstructCommunicationPattern(); 
  ///       or PMloc.ConstructCommunicationPattern_P2p(); 
  ///       or PMloc.ConstructCommunicationPattern_Collectives(); 
  ///   
  ///   Numerical preparation so that SelInv only involves Gemm.
  ///
  ///       PMloc.PreSelInv();  
  ///
  /// - Selected inversion.
  ///
  ///       PMloc.SelInv();
  ///       or PMloc.SelInv_P2p();
  ///       or PMloc.SelInv_Collectives();
  ///
  /// - Postprocessing.
  ///
  ///   Get the information in DistSparseMatrix format 
  ///
  ///       DistSparseMatrix<Scalar> Ainv;
  ///       PMloc.PMatrixToDistSparseMatrix( Ainv );  
  ///
  /// Note
  /// ----
  ///
  /// - All major operations of PMatrix, including the selected inversion
  /// are defined directly as the member function of PMatrix.
  ///
  /// - In the current version of PMatrix, square grid is assumed.  This
  /// assumption is only used when sending the information to
  /// cross-diagonal blocks, i.e. from L(isup, ksup) to U(ksup, isup).
  /// This assumption can be relaxed later.

  template<typename T>
  class PMatrixAsym: public PMatrix<T>{

    protected:
      // *********************************************************************
      // Variables
      // *********************************************************************
      // Data variables

      std::vector<std::vector<UBlock<T> > > Ucol_;
      std::vector<std::vector<LBlock<T> > > Lrow_;

      // Communication variables
      // This is the tag used for mpi communication for selinv

      enum{
        SELINV_TAG_U_SIZE,
        SELINV_TAG_U_CONTENT,
        SELINV_TAG_L_SIZE,
        SELINV_TAG_L_CONTENT,
        SELINV_TAG_UCOL_SIZE,
        SELINV_TAG_UCOL_CONTENT,
        SELINV_TAG_LROW_SIZE,
        SELINV_TAG_LROW_CONTENT,
        SELINV_TAG_L_REDUCE,
        SELINV_TAG_U_REDUCE,
        SELINV_TAG_D_SIZE,
        SELINV_TAG_D_CONTENT,
        SELINV_TAG_D_REDUCE,
        SELINV_TAG_COUNT
      };




      struct SuperNodeBufferTypeAsym:public PMatrix<T>::SuperNodeBufferType {
        NumMat<T>    UUpdateBuf;
        std::vector<Int>  ColLocalPtr;
        std::vector<Int>  BlockIdxLocalU;
        std::vector<char> SstrLrowSend;
        std::vector<char> SstrUcolSend;
        std::vector<char> SstrLrowRecv;
        std::vector<char> SstrUcolRecv;
        Int               SizeSstrLrowSend;
        Int               SizeSstrUcolSend;
        Int               SizeSstrLrowRecv;
        Int               SizeSstrUcolRecv;
      };

      /// @brief SelInvIntra_P2p
      inline void SelInvIntra_P2p(Int lidx);

      /// @brief SelInv_lookup_indexes
      inline void SelInv_lookup_indexes(SuperNodeBufferTypeAsym & snode,
                                        std::vector<LBlock<T> > & LcolRecv,
                                        std::vector<LBlock<T> > & LrowRecv,
                                        std::vector<UBlock<T> > & UcolRecv,
                                        std::vector<UBlock<T> > & UrowRecv,
                                        NumMat<T> & AinvBuf,
                                        NumMat<T> & LBuf,
                                        NumMat<T> & UBuf);

      /// @brief UnpackData
      inline void UnpackData( SuperNodeBufferTypeAsym & snode,
                              std::vector<LBlock<T> > & LcolRecv,
                              std::vector<LBlock<T> > & LrowRecv,
                              std::vector<UBlock<T> > & UcolRecv,
                              std::vector<UBlock<T> > & UrowRecv
                            );

      /// @brief ComputeDiagUpdate
      inline void ComputeDiagUpdate(SuperNodeBufferTypeAsym & snode);

      /// @brief SendRecvCD_UpdateU
      inline void SendRecvCD(
              std::vector<SuperNodeBufferTypeAsym > & arrSuperNodes,
              Int stepSuper
                            );

    public:
      // *********************************************************************
      // Public member functions 
      // *********************************************************************

      PMatrixAsym():PMatrix<T>() {}

      PMatrixAsym( const GridType* g, const SuperNodeType* s, const PEXSI::SuperLUOptions * o );

      virtual ~PMatrixAsym() {  }

      void Setup( const GridType* g, const SuperNodeType* s, const PEXSI::SuperLUOptions * o );

      /// @brief NumBlockL returns the number of nonzero L blocks for the
      /// local block column jLocal.
      //Int NumBlockL( Int jLocal ) const { return L_[jLocal].size(); }

      /// @brief NumBlockU returns the number of nonzero U blocks for the
      /// local block row iLocal.
      //Int NumBlockU( Int iLocal ) const { return U_[iLocal].size(); }


      /// @brief Lrow returns the vector of nonzero L blocks for the local
      /// block row iLocal.
      std::vector<LBlock<T> >& Lrow( Int iLocal ) { return Lrow_[iLocal]; } 	

      /// @brief Ucol returns the vector of nonzero U blocks for the local
      /// block col jLocal.
      std::vector<UBlock<T> >& Ucol( Int jLocal ) { return Ucol_[jLocal]; }


      /// @brief ConstructCommunicationPattern constructs the communication
      /// pattern to be used later in the selected inversion stage.
      /// The supernodal elimination tree is used to add an additional level of parallelism between supernodes.
      /// [ConstructCommunicationPattern_P2p](@ref PEXSI::PMatrix::ConstructCommunicationPattern_P2p) is called by default.
      virtual void ConstructCommunicationPattern( );


      /// @brief ConstructCommunicationPattern_P2p constructs the communication
      /// pattern to be used later in the selected inversion stage.
      /// The supernodal elimination tree is used to add an additional level of parallelism between supernodes.
      void ConstructCommunicationPattern_P2p( );


      /// @brief PreSelInv prepares the structure in L_ and U_ so that
      /// SelInv only involves matrix-matrix multiplication.
      ///
      /// @todo Move documentation to a more proper place and update the
      /// information.
      ///
      /// Procedure
      /// ---------
      /// PreSelInv performs
      ///
      /// - Compute the inverse of the diagonal blocks
      ///
      ///   L_{kk} <- (L_{kk} U_{kk})^{-1}
      ///
      /// - Update the lower triangular L blocks
      ///
      ///   L_{ik} <- L_{ik} L_{kk}^{-1}
      ///
      /// - Update the upper triangular U blocks which saves redundant
      /// information as in L
      ///
      ///   U_{kj} <- L_{ik}
      ///
      /// Note
      /// ----
      ///
      /// PreSelInv assumes that
      /// PEXSI::PMatrix::ConstructCommunicationPattern has been executed.
      virtual void PreSelInv( );

      /// @brief SelInv is the main function for the selected inversion.
      ///
      /// @todo Move documentation to a more proper place and update the
      /// information.
      ///
      /// Procedure
      /// ---------
      ///
      /// PSelInv is a right-looking based parallel selected inversion
      /// subroutine for sparse symmetric matrices.  Static pivoting is
      /// used in this version.
      ///
      /// Although the matrix is symmetric, the key idea of the current
      /// implementation of PSelInv is that the upper-triangular matrix is
      /// saved (in the form of UBlock).  Such redundant information is
      /// effective for reducing the complexity for designing the
      /// communication pattern.  
      ///
      /// At each supernode ksup, the lower triangular part Ainv(isup, ksup)
      /// (isup > ksup) are first updated.  The blocks in the processor row
      /// of ksup first sends the nonzero blocks in U(ksup, jsup) (which is
      /// L(isup, ksup)^T) to the Schur complements Ainv(isup, jsup).  At
      /// the same time the blocks in the processor column of ksup sends the
      /// nonzero blocks (only nonzero row indices) to the Schur complement
      /// Ainv(isup, jsup).  Then 
      ///
      /// sum_{jsup} Ainv(isup, jsup) U^{T}(ksup, jsup)
      ///
      /// is performed.  In this procedure, only processors with
      /// isRecvFromAbove[ksup] == true && isRecvFromLeft[ksup] == true
      /// participate in the computation.
      ///
      ///
      /// The result is reduced to the processor column ksup within the same
      /// processor row.  The diagonal block Ainv(ksup, ksup) is simply updated
      /// by a reduce procedure within the column processor group of ksup.
      ///
      /// Then we update the Ainv(ksup, isup) blocks, simply via the update
      /// from the cross diagonal processors.
      ///
      /// <b> NOTE </b>: The cross diagonal processor is only well here
      /// defined for square grids.  For a P x P square grid, (ip,
      /// jp) is the cross diagonal processor of (jp, ip) if ip != jp.  The
      /// current version of SelInv only works for square processor grids.
      ///
      ///
      /// Communication pattern
      /// ---------------------
      ///
      /// The communication is controlled by 3 sending varaibles and 3
      /// receiving variables. The first dimension of all the sending and
      /// receiving variables are numSuper.  The information contains
      /// redundancy since not all processors have access to all the
      /// supernodes.  However, this increases the readability of the output
      /// significantly and only increases a small amount of memory cost for
      /// indexing.  This set of sending / receiving mechanism avoids the
      /// double indexing of the supernodes and can scale to matrices of
      /// large size.
      ///
      /// - isSendToBelow:  
      ///
      ///   Dimension: numSuper x numProcRow
      ///
      ///   Role     : At supernode ksup, if isSendToBelow(ksup, ip) == true, send
      ///   all local blocks {U(ksup, jsup) | jsup > ksup} to the processor row ip.
      ///
      /// - isRecvFromAbove:
      ///
      ///   Dimension: numSuper
      ///
      ///   Role     : 
      ///
      ///     * At supernode ksup, if isRecvFromAbove(ksup) == true,
      ///       receive blocks from the processor owning the block row of ksup
      ///       within the same column processor group.
      ///
      ///     * If isRecvFromAbove(ksup) == true && isRecvFromLeft(ksup) ==
      ///     true, the ucrrent processor participate in updating Ainv(isup,
      ///     ksup).
      ///
      ///
      /// - isSendToRight:
      ///
      ///   Dimension: numSuper x numProcCol
      ///
      ///   Role     : At supernode ksup, if isSendToRight(ksup, jp) == true, send
      ///   all local blocks (mainly the nonzero row indicies, without the
      ///   values to save the communication cost) {L(isup, ksup) | isup >
      ///   ksup} to the processor column jp.
      ///
      /// - isRecvFromLeft:
      ///   
      ///   Dimension: numSuper
      ///
      ///   Role     : 
      ///
      ///     * At supernode ksup, if isRecvFromLeft(ksup) == true, receive
      ///     blocks from the processor owning the block column of ksup
      ///     within the same row processor group.
      ///
      ///     * If isRecvFromAbove(ksup) == true && isRecvFromLeft(ksup) ==
      ///     true, the ucrrent processor participate in updating Ainv(isup,
      ///     ksup).
      ///
      /// - isSendToCrossDiagonal:
      ///
      ///   Dimension: numSuper
      ///
      ///   Role     : At supernode ksup, if isSendToCrossDiagonal(ksup) ==
      ///   true, send all local blocks {(isup, ksup) | isup > ksup} to the
      ///   cross-diagonal processor.  <b> NOTE </b>: This requires a square
      ///   processor grid.
      ///
      /// - isRecvCrossDiagonal:
      ///
      ///   Dimension: numSuper
      ///
      ///   Role     : At supernode ksup, if isRecvFromCrossDiagonal(ksup) ==
      ///   true, receive from the cross-diagonal processor.  <b> NOTE </b>:
      ///   This requires a square processor grid.
      ///   
      ///
      ///
      ///
      virtual void SelInv( );
      /// @brief Point-to-point version of the selected inversion.
      void SelInv_P2p( );


      /// @brief GetDiagonal extracts the diagonal elements of the PMatrix.
      ///
      /// 1) diag is permuted back to the natural order
      ///
      /// 2) diag is shared by all processors in grid_->comm through a
      /// Allreduce procedure.
      //void GetDiagonal( NumVec<T>& diag );
      //void GetColumn	( Int colIdx,  NumVec<T>& col );

//      virtual void DumpSuperNodes(Int count){
//        Int first_snode = max(0,this->NumSuper() -count );
//        //dump the last supernodes
//        for(Int I = first_snode; I<this->NumSuper();++I){
//          statusOFS<<"****** "<<I<<" *******"<<std::endl;
//  
//          //I own blocks of that supernode
//          if(MYCOL(this->grid_) == PCOL(I,this->grid_)){
//            std::vector<LBlock<T> >&  Lcol = this->L( LBj( I, this->grid_ ) );
//            for(Int bidx = 0; bidx < Lcol.size(); ++bidx){
//              LBlock<T> & block = Lcol[bidx];
//              statusOFS<<block.blockIdx<<std::endl;
//              statusOFS<<block.nzval<<std::endl;
//            }
//
//
//            std::vector<UBlock<T> >&  Urow = this->U( LBi( I, this->grid_ ) );
//            for(Int bidx = 0; bidx < Urow.size(); ++bidx){
//              UBlock<T> & block = Urow[bidx];
//              statusOFS<<block.blockIdx<<std::endl;
//              //NumMat<T> tNzval;
//              //Transpose(block.nzval, tNzval);
//              //statusOFS<<tNzval<<std::endl;
//              statusOFS<<block.nzval<<std::endl;
//            }
//
//
//
//
//          }
//        }
//      }
//

  };



} // namespace PEXSI


#include "pselinv_asym_impl.hpp"

#endif //_PEXSI_PSELINV_ASYM_HPP_

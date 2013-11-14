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
/// @file pselinv.hpp
/// @brief Main file for parallel selected inversion.
/// @date 2013-08-05
#ifndef _PSELINV_HPP_
#define _PSELINV_HPP_

// *********************************************************************
//  Common utilities
// *********************************************************************

#include  "environment.hpp"
#include	"numvec_impl.hpp"
#include	"nummat_impl.hpp" 
#include  "sparse_matrix_impl.hpp"
#include  "superlu_dist_interf.hpp"
#include  "mpi_interf.hpp"
#include	"utility.hpp"
#include	"blas.hpp"
#include	"lapack.hpp"
#include <set>

namespace PEXSI{

  struct SuperLUOptions;

  typedef std::vector<bool> bitMask;
  typedef std::map<bitMask , std::vector<Int> > bitMaskSet;

#ifdef SANITY_CHECK
  struct SelInvError{
    Real Value;
    Int ksup;
    Int ib;
    Int i;
    Int j;

    SelInvError(Real val):Value(val), ksup(-1),ib(-1),i(-1),j(-1){}
    SelInvError(Real val, Int pksup, Int pib, Int pi, Int pj):Value(val), ksup(pksup),ib(pib),i(pi),j(pj){}
    inline void Set(Real val, Int pksup, Int pib, Int pi, Int pj) {Value = val;ksup=pksup;ib=pib;i=pi;j=pj;}

    std::ostream& print(std::ostream &o) const
    {
      return o << "Value = " << Value << " (ksup="<<ksup<<", ib="<<ib<<", i="<<i<<", j="<<j<<")";
    }

  };


  inline std::ostream& operator << (std::ostream &o,const SelInvError &a){
    return a.print(o);
  }

  struct SelInvErrors{
    SelInvError MaxRelError;
    SelInvError CorrAbsError;
    SelInvError MaxAbsError;


    SelInvError MaxNwiseRelError;
    SelInvError CorrNwiseAbsError;
    SelInvError MaxNwiseAbsError;

    SelInvError MaxRwiseRelError;
    SelInvError CorrRwiseAbsError;
    SelInvError MaxRwiseAbsError;

    SelInvError MaxCwiseRelError ;
    SelInvError CorrCwiseAbsError;
    SelInvError MaxCwiseAbsError;



    // Member functions to setup the default value
    SelInvErrors(): MaxRelError(0),CorrAbsError(0),MaxAbsError(0),
    MaxNwiseRelError(0), CorrNwiseAbsError(0), MaxNwiseAbsError(0),
    MaxRwiseRelError(0), CorrRwiseAbsError(0), MaxRwiseAbsError(0),
    MaxCwiseRelError (0), CorrCwiseAbsError(0), MaxCwiseAbsError(0) {}

    std::ostream & print(std::ostream & output) const
    {
      output <<std::endl<< "Element-wise errors:"<<std::endl;
      output << "Max relative error = " << MaxRelError << std::endl;
      output << "Corresp. absolute error = " << CorrAbsError << std::endl;
      output << "Max absolute = " << MaxAbsError << std::endl;

      output <<std::endl<< "Norm-wise errors:"<<std::endl;
      output << "Max relative error = " << MaxNwiseRelError << std::endl;
      output << "Corresp. absolute error = " << CorrNwiseAbsError << std::endl;
      output << "Max absolute = " << MaxNwiseAbsError << std::endl;

      output <<std::endl<< "Row-wise errors:"<<std::endl;
      output << "Max relative error = " << MaxRwiseRelError << std::endl;
      output << "Corresp. absolute error = " << CorrRwiseAbsError << std::endl;
      output << "Max absolute = " << MaxRwiseAbsError << std::endl;

      output <<std::endl<< "Column-wise errors:"<<std::endl;
      output << "Max relative error = " << MaxCwiseRelError << std::endl;
      output << "Corresp. absolute error = " << CorrCwiseAbsError << std::endl;
      output << "Max absolute = " << MaxCwiseAbsError << std::endl;
      return output;
    }

  };


  inline std::ostream& operator << (std::ostream &o,const SelInvErrors &a){
    return a.print(o);
  }
#endif



  /**********************************************************************
   * Basic PSelInv data structure
   **********************************************************************/

  /// @struct GridType
  ///
  /// @brief GridType is the PSelInv way of defining the grid.  
  ///
  /// GridType should be consistent with the grid used by SuperLU.
  ///
  /// NOTE: It is your responsibility to make sure that the SuperLUGrid
  /// and GridType used for SelInv are the same.
  struct GridType{
    // Data
    MPI_Comm    comm;
    MPI_Comm    rowComm;
    MPI_Comm    colComm;
    Int         mpirank;
    Int         mpisize; 
    Int         numProcRow;
    Int         numProcCol;

    // Member function
    GridType( MPI_Comm Bcomm, int nprow, int npcol );
    ~GridType();
  };

  /// @struct SuperNodeType
  ///
  /// @brief SuperNodeType describes mapping between supernode and column, the
  /// permutation information, and potentially the elimination tree (not
  /// implemented here).
  /// 
  /// superIdx[i] is the supernode index to which column i belongs. 
  /// This is the same as supno[i] in SuperLU.
  ///
  /// superPtr[s] is the leading column of the s-th supernode (as in
  /// colptr).  This is the same as xsup[s] in SuperLU.
  ///
  ///	e.g.   superIdx  0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
  ///	       superPtr  0 1 2 4 7 12
  ///
  /// This is allocated during symbolic factorization SYMBFACT.
  ///
  /// perm is the permutation vector.  Symmetric permutation is assumed.
  /// perm is the same as ScalePermstruct -> perm_c.
  ///
  /// permInv is the inverse of the permutation vector.
  ///
  struct SuperNodeType{
    IntNumVec   perm;              
    IntNumVec   permInv;
    IntNumVec   superIdx;
    IntNumVec   superPtr;
    IntNumVec   etree;
  };


  /// @struct LBlock
  ///
  /// @brief LBlock stores a nonzero block in the lower triangular part or
  /// the diagonal part in PSelInv.
  struct LBlock{
    // Variables
    /// @brief Block index (supernodal index)
    Int               blockIdx;

    /// @brief Number of nonzero rows. 
    Int               numRow;

    /// @brief Number of nonzero columns. 
    Int               numCol;

    /// @brief Dimension numRow * 1, index (0-based) for the number of nonzero rows.
    IntNumVec         rows;

    /// @brief Dimension numRow * numCol, nonzero elements.
    NumMat<Scalar>    nzval;

    // Member functions;
    LBlock() {blockIdx = -1; numRow = 0; numCol =0;}
    ~LBlock() {}
    LBlock& operator = (const LBlock& LB) {
      blockIdx    = LB.blockIdx;
      numRow      = LB.numRow;
      numCol      = LB.numCol;
      rows        = LB.rows;
      nzval       = LB.nzval;
      return *this;
    }
    friend std::ostream& operator<<(std::ostream& out, const LBlock& vec) // output
    {
      out << "(" << vec.blockIdx << ", " << vec.numRow << ", " << vec.numCol <<std::endl<< "rows " << vec.rows <<std::endl<< "nzval " <<std::endl<< vec.nzval << ")";
      return out;
    }


  };

  /// @struct UBlock
  ///
  /// @brief UBlock stores a nonzero block in the upper triangular part in PSelInv. 
  /// 
  /// In particular, the current version of PSelInv is for sparse
  /// symmetric matrices.  All UBlocks, labeled as U(i,j), i<j save the
  /// redundant information as saved in L(j,i). The purpose of having
  /// UBlocks is to facilitate the communication.  
  ///
  /// @see PMatrix::SelInv
  struct UBlock{
    // Variables
    /// @brief Block index (supernodal index)
    Int               blockIdx;

    /// @brief Number of nonzero rows. 
    Int               numRow;

    /// @brief Number of nonzero columns. 
    Int               numCol;

    /// @brief Dimension numRow * 1, index (0-based) for the number of nonzero rows.
    IntNumVec         cols;

    /// @brief Dimension numRow * numCol, nonzero elements.
    NumMat<Scalar>    nzval;

    // Member functions;
    UBlock() {blockIdx = -1; numRow = 0; numCol =0;}
    ~UBlock() {}
    UBlock& operator = (const UBlock& UB) {
      blockIdx    = UB.blockIdx;
      numRow      = UB.numRow;
      numCol      = UB.numCol;
      cols        = UB.cols;
      nzval       = UB.nzval;
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const UBlock& vec) // output
    {
      out << "(" << vec.blockIdx << ", " << vec.numRow << ", " << vec.numCol <<std::endl<< "cols " << vec.cols <<std::endl<< "nzval " <<std::endl<< vec.nzval << ")";
      return out;
    }


  };

  // *********************************************************************
  // SuperLU style utility functions
  // 
  // The SuperLU style macros are defined here as inline functions 
  // so that the code is more portable.
  // *********************************************************************

  /// @brief MYPROC returns the current processor rank.
  inline Int MYPROC( const GridType* g )
  { return g->mpirank; }

  /// @brief MYROW returns my processor row
  inline Int MYROW( const GridType* g )
  { return g->mpirank / g->numProcCol; }

  /// @brief MYCOL returns my processor column
  inline Int MYCOL( const GridType* g )
  { return g->mpirank % g->numProcCol; }

  /// @brief PROW returns the processor row that the bnum-th block
  /// (supernode) belongs to.
  inline Int PROW( Int bnum, const GridType* g ) 
  { return bnum % g->numProcRow; }

  /// @brief PCOL returns the processor column that the bnum-th block
  /// (supernode) belongs to.
  inline Int PCOL( Int bnum, const GridType* g ) 
  { return bnum % g->numProcCol; }

  /// @brief PNUM returns the processor rank that the bnum-th block
  /// (supernode) belongs to.
  inline Int PNUM( Int i, Int j, const GridType* g )
  { return  (i%g->numProcRow) * g->numProcCol + j%g->numProcCol; }

  /// @brief LBi returns the local block number on the processor at
  /// processor row PROW( bnum, g ).
  inline Int LBi( Int bnum, const GridType* g )
  { return bnum / g->numProcRow; }

  /// @brief LBj returns the local block number on the processor at
  /// processor column PCOL( bnum, g ).
  inline Int LBj( Int bnum, const GridType* g)
  { return bnum / g->numProcCol; }

  /// @brief GBi returns the global block number from a local block number
  /// in the row direction.
  inline Int GBi( Int iLocal, const GridType* g )
  { return iLocal * g->numProcRow + MYROW( g ); }

  /// @brief GBj returns the global block number from a local block number
  /// in the column direction.
  inline Int GBj( Int jLocal, const GridType* g )
  { return jLocal * g->numProcCol + MYCOL( g ); }

  /// @brief CEILING is used for computing the storage space for local
  /// number of blocks.
  inline Int CEILING( Int a, Int b )
  { return (a%b) ? ( a/b + 1 ) : ( a/b ); }

  /// @brief BlockIdx returns the block index of a column i.
  inline Int BlockIdx( Int i, const SuperNodeType *s )
  { return s->superIdx[i]; }

  /// @brief FirstBlockCol returns the first column of a block
  /// bnum.
  inline Int FirstBlockCol( Int bnum, const SuperNodeType *s )
  { return s->superPtr[bnum]; }	


  /// @brief FirstBlockRow returns the first column of a block
  /// bnum. Note: the functionality of FirstBlockRow is exactly the same
  /// as in FirstBlockCol.
  inline Int FirstBlockRow( Int bnum, const SuperNodeType *s )
  { return s->superPtr[bnum]; } 


  /// @brief SuperSize returns the size of the block bnum.
  inline Int SuperSize( Int bnum, const SuperNodeType *s )
  { return s->superPtr[bnum+1] - s->superPtr[bnum]; } 

  /// @brief NumSuper returns the total number of supernodes.
  inline Int NumSuper( const SuperNodeType *s )
  { return s->superPtr.m() - 1; }

  /// @brief NumCol returns the total number of columns for a supernodal
  /// partiiton.
  inline Int NumCol( const SuperNodeType *s )
  { return s->superIdx.m(); }


  // *********************************************************************
  // Serialize / Deserialize
  // *********************************************************************

  // L part

  /// @namespace LBlockMask 
  ///
  /// @brief LBlockMask allows one to compress the selected data in
  /// LBlock used for communication.
  ///
  /// Example
  /// -------
  /// std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
  ///
  /// assumes all information is to be communicated. 
  ///
  ///	mask[LBlockMask::NZVAL] = 0;
  ///
  ///	then nzval will not be serialized / deserialized.
  ///
  namespace LBlockMask{
    enum {
      BLOCKIDX,
      NUMROW,
      NUMCOL,
      ROWS,
      NZVAL,
      TOTAL_NUMBER
    };
  }

  Int inline serialize(LBlock& val, std::ostream& os, const std::vector<Int>& mask){
    Int i = 0;
    if(mask[i]==1) serialize(val.blockIdx, os, mask); i++;
    if(mask[i]==1) serialize(val.numRow,  os, mask); i++;
    if(mask[i]==1) serialize(val.numCol,  os, mask); i++;
    if(mask[i]==1) serialize(val.rows, os, mask);   i++;
    if(mask[i]==1) serialize(val.nzval, os, mask);  i++;
    return 0;
  }

  Int inline deserialize(LBlock& val, std::istream& is, const std::vector<Int>& mask){
    Int i = 0;
    if(mask[i]==1) deserialize(val.blockIdx, is, mask); i++;
    if(mask[i]==1) deserialize(val.numRow,  is, mask); i++;
    if(mask[i]==1) deserialize(val.numCol,  is, mask); i++;
    if(mask[i]==1) deserialize(val.rows,   is, mask); i++;
    if(mask[i]==1) deserialize(val.nzval,  is, mask); i++; 
    return 0;
  }

  // U part
  /// @namespace UBlockMask 
  ///
  /// @brief UBlockMask allows one to compress the selected data in
  /// UBlock used for communication.
  ///
  /// Example
  /// -------
  ///
  /// std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
  ///
  /// assumes all information is to be communicated. 
  ///
  ///	mask[UBlockMask::NZVAL] = 0;
  ///
  ///	then nzval will not be serialized / deserialized.
  namespace UBlockMask{
    enum {
      BLOCKIDX,
      NUMROW,
      NUMCOL,
      COLS,
      NZVAL,
      TOTAL_NUMBER
    };
  }

  Int inline serialize(UBlock& val, std::ostream& os, const std::vector<Int>& mask){
    Int i = 0;
    if(mask[i]==1) serialize(val.blockIdx, os, mask); i++;
    if(mask[i]==1) serialize(val.numRow,  os, mask); i++;
    if(mask[i]==1) serialize(val.numCol,  os, mask); i++;
    if(mask[i]==1) serialize(val.cols, os, mask);   i++;
    if(mask[i]==1) serialize(val.nzval, os, mask);  i++;
    return 0;
  }

  Int inline deserialize(UBlock& val, std::istream& is, const std::vector<Int>& mask){
    Int i = 0;
    if(mask[i]==1) deserialize(val.blockIdx, is, mask); i++;
    if(mask[i]==1) deserialize(val.numRow,  is, mask); i++;
    if(mask[i]==1) deserialize(val.numCol,  is, mask); i++;
    if(mask[i]==1) deserialize(val.cols,   is, mask); i++;
    if(mask[i]==1) deserialize(val.nzval,  is, mask); i++; 
    return 0;
  }


  /**********************************************************************
   * Main data structure in PSelInv: PMatrix
   **********************************************************************/

  /// @class PMatrix
  ///
  /// @brief PMatrix contains the main data structure and the
  /// computational routine for the parallel selected inversion.  
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
  ///   
  ///   Numerical preparation so that SelInv only involves Gemm.
  ///
  ///       PMloc.PreSelInv();  
  ///
  /// - Selected inversion.
  ///
  ///       PMloc.SelInv();
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

  class PMatrix{

    private:
      // *********************************************************************
      // Variables
      // *********************************************************************
      // Data variables

      const GridType*       grid_;

      const SuperNodeType*  super_;

      const SuperLUOptions * options_;

      std::vector<std::vector<LBlock> > L_;
      std::vector<std::vector<UBlock> > U_;

      std::vector<std::vector<Int> > workingSet_;


      // Communication variables
      NumMat<bool>                       isSendToBelow_;
      NumMat<bool>                       isSendToRight_;
      NumVec<bool>                       isSendToDiagonal_;
      NumMat<bool>                       isSendToCrossDiagonal_;

      NumVec<bool>                       isRecvFromAbove_;
      NumMat<bool>                       isRecvFromBelow_;
      NumVec<bool>                       isRecvFromLeft_;
      NumMat<bool>                       isRecvFromCrossDiagonal_;

      //Communicators for the Bcast variant

      IntNumVec maxCommSizes_;

      NumVec<Int>                       countSendToBelow_;
      bitMaskSet                        maskSendToBelow_;
      std::vector<MPI_Comm>             commSendToBelow_;
      std::vector<MPI_Comm*>            commSendToBelowPtr_;
      std::vector<Int>                  commSendToBelowRoot_;
      std::vector<bitMask *>         commSendToBelowMaskPtr_;
      std::vector<bitMask>         commSendToBelowMask_;

      NumVec<Int>                       countRecvFromBelow_;
      bitMaskSet                        maskRecvFromBelow_;
      std::vector<MPI_Comm>             commRecvFromBelow_;
      std::vector<MPI_Comm*>            commRecvFromBelowPtr_;
      std::vector<Int>                  commRecvFromBelowRoot_;

      std::vector<bitMask *>         commRecvFromBelowMaskPtr_;
      std::vector<bitMask>         commRecvFromBelowMask_;
      

      NumVec<Int>                       countSendToRight_;
      bitMaskSet                        maskSendToRight_;
      std::vector<MPI_Comm>             commSendToRight_;
      std::vector<MPI_Comm*>            commSendToRightPtr_;
      std::vector<Int>                  commSendToRightRoot_;
      std::vector<bitMask *>         commSendToRightMaskPtr_;
      std::vector<bitMask>         commSendToRightMask_;

      //NumVec<Int>                       countCrossDiag_;

      // This is the tag used for mpi communication for selinv
      enum{
        SELINV_TAG_U_SIZE,
        SELINV_TAG_U_CONTENT,
        SELINV_TAG_L_SIZE,
        SELINV_TAG_L_CONTENT,
        SELINV_TAG_L_REDUCE,
        SELINV_TAG_D_SIZE,
        SELINV_TAG_D_CONTENT,
        SELINV_TAG_D_REDUCE,
        SELINV_TAG_COUNT
      };





      struct SuperNodeBufferType{
        NumMat<Scalar>    LUpdateBuf;
        NumMat<Scalar>    DiagBuf;
        std::vector<Int>  RowLocalPtr;
        std::vector<Int>  BlockIdxLocal;
        std::vector<char> SstrLcolSend;
        std::vector<char> SstrUrowSend;
        std::vector<char> SstrLcolRecv;
        std::vector<char> SstrUrowRecv;
        Int               SizeSstrLcolSend;
        Int               SizeSstrUrowSend;
        Int               SizeSstrLcolRecv;
        Int               SizeSstrUrowRecv;
        Int               Index;
        Int               isReady;

        SuperNodeBufferType():
          SizeSstrLcolSend(0),
          SizeSstrUrowSend(0),
          SizeSstrLcolRecv(0),
          SizeSstrUrowRecv(0),
          Index(0), 
          isReady(0){}

        SuperNodeBufferType(Int &pIndex) :
          SizeSstrLcolSend(0),
          SizeSstrUrowSend(0),
          SizeSstrLcolRecv(0),
          SizeSstrUrowRecv(0),
          Index(pIndex),
          isReady(0) {}

      };

      inline void PMatrix::SelInvIntra_Collectives(Int lidx);

      inline void PMatrix::SelInvIntra_P2p(Int lidx);

      inline void PMatrix::SelInv_lookup_indexes(const Int ksup, std::vector<LBlock> & LcolRecv, std::vector<UBlock> & UrowRecv, NumMat<Scalar> & AinvBuf,NumMat<Scalar> & UBuf,NumMat<Scalar> & LUpdateBuf);

      inline void PMatrix::SelInv_lookup_indexes(SuperNodeBufferType & snode, std::vector<LBlock> & LcolRecv, std::vector<UBlock> & UrowRecv, NumMat<Scalar> & AinvBuf,NumMat<Scalar> & UBuf);

      inline void PMatrix::GetWorkSet(std::vector<Int> & snodeEtree, std::vector<std::vector<Int> > & WSet);

      inline void PMatrix::UnpackData(SuperNodeBufferType & snode, std::vector<LBlock> & LcolRecv, std::vector<UBlock> & UrowRecv);

      inline void PMatrix::ComputeDiagUpdate(SuperNodeBufferType & snode);

      inline void PMatrix::SendRecvCD_UpdateU(std::vector<SuperNodeBufferType> & arrSuperNodes, Int stepSuper);

      void PMatrix::getMaxCommunicatorSizes();














    public:
      // *********************************************************************
      // Public member functions 
      // *********************************************************************

      PMatrix( const GridType* g, const SuperNodeType* s, const PEXSI::SuperLUOptions * o );

      ~PMatrix();

      Int NumCol() const { return super_ -> superIdx.m(); }

      Int NumSuper() const { return super_ ->superPtr.m() - 1; }

      /// @brief NumLocalBlockCol returns the total number of block columns.
      Int NumLocalBlockCol() const { return CEILING( this->NumSuper(), grid_->numProcCol ); }

      /// @brief NumLocalBlockRow returns the total number of block rows.
      Int NumLocalBlockRow() const { return CEILING( this->NumSuper(), grid_->numProcRow); }

      /// @brief NumBlockL returns the number of nonzero L blocks for the
      /// local block column jLocal.
      Int NumBlockL( Int jLocal ) const { return L_[jLocal].size(); }

      /// @brief NumBlockU returns the number of nonzero U blocks for the
      /// local block row iLocal.
      Int NumBlockU( Int iLocal ) const { return U_[iLocal].size(); }

      /// @brief Grid returns the GridType structure of the current PMatrix.
      const GridType* Grid() const { return grid_; }

      /// @brief SuperNode returns the supernodal partition of the current
      /// PMatrix.
      const SuperNodeType* SuperNode() const { return super_; }	

      /// @brief L returns the vector of nonzero L blocks for the local
      /// block column jLocal.
      std::vector<LBlock>& L( Int jLocal ) { return L_[jLocal]; } 	

      /// @brief U returns the vector of nonzero U blocks for the local
      /// block row iLocal.
      std::vector<UBlock>& U( Int iLocal ) { return U_[iLocal]; }

      /// @brief WorkingSet returns the ordered list of supernodes which could
      /// be done in parallel.
      std::vector<std::vector<int> >& WorkingSet( ) { return workingSet_; } 	

      /// @brief CountSendToRight returns the number of processors 
      /// to the right of current processor with which it has to communicate
      Int CountSendToRight(Int ksup) {  Int count= std::count (isSendToRight_.VecData(ksup), isSendToRight_.VecData(ksup) + grid_->numProcCol, true); return (isSendToRight_(MYCOL(grid_),ksup)?count-1:count); }

      /// @brief CountRecvFromBelow returns the number of processors 
      /// below the current processor from which it receives data
      Int CountRecvFromBelow(Int ksup) {  Int count= std::count (isRecvFromBelow_.VecData(ksup), isRecvFromBelow_.VecData(ksup) + grid_->numProcRow, true); return (isRecvFromBelow_(MYROW(grid_),ksup)?count-1:count); }

      /// @brief CountSendToCrossDiagonal returns the number of cross diagonal
      /// processors with which current processor has to communicate
      Int CountSendToCrossDiagonal(Int ksup) {  Int count= std::count (isSendToCrossDiagonal_.VecData(ksup), isSendToCrossDiagonal_.VecData(ksup) + grid_->numProcCol, true);  return ((isSendToCrossDiagonal_(MYCOL(grid_),ksup) && MYROW(grid_)==PROW(ksup,grid_))?count-1:count); }

      /// @brief CountRecvFromCrossDiagonal returns the number of cross diagonal
      /// processors with which current processor has to communicate
      Int CountRecvFromCrossDiagonal(Int ksup) {  Int count= std::count (isRecvFromCrossDiagonal_.VecData(ksup), isRecvFromCrossDiagonal_.VecData(ksup) + grid_->numProcRow, true);  return ((isRecvFromCrossDiagonal_(MYROW(grid_),ksup) && MYCOL(grid_)==PCOL(ksup,grid_))?count-1:count); }




      /// @brief GetEtree computes the supernodal elimination tree
      /// to be used later in the pipelined selected inversion stage.
      void GetEtree(std::vector<Int> & etree_supno );

      /// @brief ConstructCommunicationPattern_P2p constructs the communication
      /// pattern to be used later in the selected inversion stage.
      /// The supernodal elimination tree is used to schedule the pipelined supernodes.
      void ConstructCommunicationPattern_P2p( );

      /// @brief ConstructCommunicationPattern_Collectives constructs the communication
      /// pattern to be used later in the selected inversion stage with the Bcast variant.
      /// The supernodal elimination tree is used to schedule the pipelined supernodes.
      void ConstructCommunicationPattern_Collectives( );

      /// @brief DestructCommunicators_Collectives frees the MPI communicators allocated
      /// by CreateCommunicators_Collectives.
      void DestructCommunicators_Collectives( );

      void PMatrix::ConstructCommunicators_Collectives(Int lidx);

      /// @brief ConstructCommunicationPattern constructs the communication
      /// pattern to be used later in the selected inversion stage.
      void ConstructCommunicationPattern( );


      /// @brief PreSelInv prepares the structure in L_ and U_ so that
      /// SelInv only involves matrix-matrix multiplication.
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
      void PreSelInv( );

      /// @brief SelInv is the main function for the selected inversion.
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
      /// Future work
      /// -----------
      ///
      /// - Pipelining:
      ///
      ///   Approximately 50% of speedup.
      ///   
      ///   Reference: S. Li and J. Demmel, SuperLU_DIST : A Scalable
      ///   Distributed-Memory Sparse Direct Solver for Unsymmetric Linear
      ///   Systems, ACM TOMS, 2003
      ///
      /// - Look-ahead and static scheduling:
      ///
      ///   Approximately 2.5~3.5 times speedup.
      ///
      ///   Reference: I. Yamazaki and S. Li, New Scheduling Strategies and Hybrid Programming for
      ///   a Parallel Right-looking Sparse LU Factorization Algorithm on
      ///   Multicore Cluster Systems, IPDPS 2012
      ///
      ///
      void SelInv( );
      void SelInv_Collectives( );
      void SelInv_P2p( );
      void SelInv_Hybrid(Int threshold);


      /// @brief GetDiagonal extracts the diagonal elements of the PMatrix.
      ///
      /// 1) diag is permuted back to the natural order
      ///
      /// 2) diag is shared by all processors in grid_->comm through a
      /// Allreduce procedure.
      void GetDiagonal( NumVec<Scalar>& diag );
      void GetColumn	( Int colIdx,  NumVec<Scalar>& col );

#ifdef SANITY_CHECK
      void CompareDiagonal	( PMatrix & Ref, SelInvErrors & errors);
      void CompareOffDiagonal	( PMatrix & Ref,SelInvErrors & errors);
#endif

      /// @brief PMatrixToDistSparseMatrix converts the PMatrix into a
      /// distributed compressed sparse column matrix format.
      /// The DistSparseMatrix follows the natural order.
      ///
      /// @param[out] A Output sparse matrix.
      void PMatrixToDistSparseMatrix( DistSparseMatrix<Scalar>& A );

      /// @brief PMatrixToDistSparseMatrix converts the PMatrix into a
      /// distributed compressed sparse column matrix format B, which has
      /// the same sparsity pattern as the DistSparseMatrix A.
      ///
      /// The DistSparseMatrix follows the natural order.
      ///
      /// @param[in]  A Input sparse matrix to provide the sparsity pattern.
      ///
      /// @param[out] B Output sparse matrix.
      void PMatrixToDistSparseMatrix( 
          const DistSparseMatrix<Scalar>& A,
          DistSparseMatrix<Scalar>& B	);


      /// @brief PMatrixToDistSparseMatrix2 is a more efficient version
      /// which performs the same job as PMatrixToDistSparseMatrix(A,B)
      /// especially when A contains much less nonzero elements than the
      /// current PMatrix.
      ///
      /// The DistSparseMatrix follows the natural order.
      ///
      /// @param[in]  A Input sparse matrix to provide the sparsity pattern.
      ///
      /// @param[out] B Output sparse matrix.
      void PMatrixToDistSparseMatrix2( 
          const DistSparseMatrix<Scalar>& A,
          DistSparseMatrix<Scalar>& B );


      /// @brief NnzLocal computes the number of nonzero elements (L and U)
      /// saved locally.
      Int  NnzLocal();

      /// @brief Nnz computes the total number of nonzero elements in the
      /// PMatrix.
      LongInt  Nnz();

      /// @brief GetNegativeInertia computes the negative inertia of a
      /// PMatrix.  This can be used to estimate e.g. the number of
      /// eigenvalues of a matrix below a certain threshold.
      void GetNegativeInertia	( Real& inertia );

  };



} // namespace PEXSI

#endif // _PSELINV_HPP_

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
/// @file pselinv_asym_impl.hpp
/// @brief Implementation of the parallel SelInv.
/// @date 2013-08-05
#ifndef _PEXSI_PSELINV_ASYM_IMPL_HPP_
#define _PEXSI_PSELINV_ASYM_IMPL_HPP_

#include "timer.h"
#include "superlu_dist_interf.hpp"

namespace PEXSI{
  template<typename T>
    PMatrixAsym<T>::PMatrixAsym ( 
        const GridType* g, 
        const SuperNodeType* s, 
        const PEXSI::SuperLUOptions * o 
        )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::PMatrix");
#endif

      this->Setup( g, s, o );

#ifndef _RELEASE_
      PopCallStack();
#endif
      return ;
    } 		// -----  end of method PMatrixAsym::PMatrixAsym  ----- 


  template<typename T>
    void PMatrixAsym<T>::Setup( 
        const GridType* g, 
        const SuperNodeType* s, 
        const PEXSI::SuperLUOptions * o 
        ) {
#ifndef _RELEASE_
      PushCallStack("PMatrix::Setup");
#endif

      PMatrix<T>::Setup(g,s,o);



      Lrow_.clear();
      Ucol_.clear();

      Lrow_.resize( this->NumLocalBlockRow() );
      Ucol_.resize( this->NumLocalBlockRow() );

#ifndef _RELEASE_
      PopCallStack();
#endif
      return ;
    } 		// -----  end of method PMatrixAsym::Setup   ----- 


  ///////////// Utility functions ///////////////////
  template<typename T>
    inline  void PMatrixAsym<T>::SelInv_lookup_indexes(
        SuperNodeBufferTypeAsym & snode, 
        std::vector<LBlock<T> > & LcolRecv, 
        std::vector<LBlock<T> > & LrowRecv, 
        std::vector<UBlock<T> > & UcolRecv, 
        std::vector<UBlock<T> > & UrowRecv, /*useless so far*/ 
        NumMat<T> & AinvBuf,
        NumMat<T> & LBuf,
        NumMat<T> & UBuf )
    {
      TIMER_START(Compute_Sinv_LT_Lookup_Indexes);

      TIMER_START(Build_colptr_rowptr);
      // rowPtrL[ib] gives the row index in snode.LUpdateBuf for the first
      // nonzero row in LcolRecv[ib]. The total number of rows in
      // snode.LUpdateBuf is given by rowPtr[end]-1
      std::vector<Int> rowPtrL(LcolRecv.size() + 1);
      // colPtrL[jb] gives the column index in LBuf for the first
      // nonzero column in UrowRecv[jb]. The total number of rows in
      // LBuf is given by colPtr[end]-1
      std::vector<Int> colPtrL(LrowRecv.size() + 1);

      rowPtrL[0] = 0;
      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
        rowPtrL[ib+1] = rowPtrL[ib] + LcolRecv[ib].numRow;
      }
      colPtrL[0] = 0;
      for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
        colPtrL[jb+1] = colPtrL[jb] + LrowRecv[jb].numCol;
      }


      statusOFS<<"UrowRecv blockIdx: "<<std::endl;
      for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
        statusOFS<<UrowRecv[jb].blockIdx<<" ";
      }
      statusOFS<<std::endl;

      statusOFS<<"UcolRecv blockIdx: "<<std::endl;
      for( Int jb = 0; jb < UcolRecv.size(); jb++ ){
        statusOFS<<UcolRecv[jb].blockIdx<<" ";
      }
      statusOFS<<std::endl;

      statusOFS<<"LcolRecv blockIdx: "<<std::endl;
      for( Int jb = 0; jb < LcolRecv.size(); jb++ ){
        statusOFS<<LcolRecv[jb].blockIdx<<" ";
      }
      statusOFS<<std::endl;


      statusOFS<<"LrowRecv blockIdx: "<<std::endl;
      for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
        statusOFS<<LrowRecv[jb].blockIdx<<" ";
      }
      statusOFS<<std::endl;








      std::vector<Int> colPtrU(UcolRecv.size() + 1);
//      for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
//        UBlock<T>& UB = UrowRecv[jb];
//        //Find the block in LrowRecv and put it at the same index
//        Int indexL = 0;
//        for( Int ib = 0; ib < LrowRecv.size(); ib++ ){
//          LBlock<T>& LB = LrowRecv[ib];
//          if(LB.blockIdx==UB.blockIdx){
//            indexL = ib;
//            break;
//          }
//        }
//        
//        colPtrU[jb] = colPtrL[indexL];
//      }
//      colPtrU[UrowRecv.size()] = colPtrL[UrowRecv.size()];

      colPtrU[0] = 0;
      for( Int jb = 0; jb < UcolRecv.size(); jb++ ){
        colPtrU[jb+1] = colPtrU[jb] + UcolRecv[jb].numCol;
      }


      Int numRowAinvBuf = *rowPtrL.rbegin();
      Int numColAinvBuf = *colPtrL.rbegin();
      TIMER_STOP(Build_colptr_rowptr);

      TIMER_START(Allocate_lookup);
      // Allocate for the computational storage
      AinvBuf.Resize( numRowAinvBuf, numColAinvBuf );
      LBuf.Resize( SuperSize( snode.Index, this->super_ ), numColAinvBuf );
      UBuf.Resize( SuperSize( snode.Index, this->super_ ), numColAinvBuf );
      TIMER_STOP(Allocate_lookup);

      TIMER_START(Fill_LBuf);
      // Fill LBuf first. Make the transpose later in the Gemm phase.
      for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
        LBlock<T>& LB = LrowRecv[jb];
        if( LB.numRow != SuperSize(snode.Index, this->super_) ){
#ifdef USE_ABORT
          abort();
#endif
          throw std::logic_error( "The size of LB is not right.  Something is seriously wrong." );
        }
        lapack::Lacpy( 'A', LB.numRow, LB.numCol, LB.nzval.Data(),
            LB.numRow, LBuf.VecData( colPtrL[jb] ), SuperSize( snode.Index, this->super_ ) );


        statusOFS<<"LB = "<<std::endl<<LB<<std::endl;
        statusOFS<<"LBuf = "<<std::endl<<LBuf<<std::endl;

      }
      TIMER_STOP(Fill_LBuf);

      TIMER_START(Fill_UBuf);
      // Fill UBuf first. 
      SetValue(UBuf, ZERO<T>());
      for( Int jb = 0; jb < UcolRecv.size(); jb++ ){
        UBlock<T>& UB = UcolRecv[jb];
        if( UB.numRow != SuperSize(snode.Index, this->super_) ){
#ifdef USE_ABORT
          abort();
#endif
          throw std::logic_error( "The size of UB is not right.  Something is seriously wrong." );
        }


        lapack::Lacpy( 'A', UB.numRow, UB.numCol, UB.nzval.Data(),
            UB.numRow, UBuf.VecData( colPtrU[jb] ), SuperSize( snode.Index, this->super_ ) );

        statusOFS<<"UB = "<<std::endl<<UB<<std::endl;
        statusOFS<<"UBuf = "<<std::endl<<UBuf<<std::endl;

      }
      TIMER_STOP(Fill_UBuf);

      // Calculate the relative indices for (isup, jsup)
      // Fill AinvBuf with the information in L or U block.
      TIMER_START(JB_Loop);


#ifdef SORT
    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      LBlock<T>& LrowB = LrowRecv[jb];
      Int jsup = LrowB.blockIdx;

      Int SinvColsSta = FirstBlockCol( jsup, this->super_ );

      // Column relative indicies
      std::vector<Int> relCols( LrowB.numCol );
      for( Int j = 0; j < LrowB.numCol; j++ ){
        relCols[j] = LrowB.rows[j] - SinvColsSta;
      }


      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
        LBlock<T>& LB = LcolRecv[ib];
        Int isup = LB.blockIdx;
        Int SinvRowsSta = FirstBlockCol( isup, this->super_ );
        T* nzvalAinv = &AinvBuf( rowPtrL[ib], colPtrL[jb] );
        Int     ldAinv    = numRowAinvBuf;

        // Pin down the corresponding block in the part of Sinv.
        if( isup >= jsup ){


          std::vector<LBlock<T> >&  LcolSinv = this->L( LBj(jsup, this->grid_ ) );
          bool isBlockFound = false;
          TIMER_START(PARSING_ROW_BLOCKIDX);
          for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( LcolSinv[ibSinv].blockIdx == isup ){
              LBlock<T>& SinvB = LcolSinv[ibSinv];

              // Row relative indices
              std::vector<Int> relRows( LB.numRow );
              Int* rowsLBPtr    = LB.rows.Data();
              Int* rowsSinvBPtr = SinvB.rows.Data();

              TIMER_START(STDFIND_ROW);
              Int * pos =&rowsSinvBPtr[0];
              Int * last =&rowsSinvBPtr[SinvB.numRow];
              for( Int i = 0; i < LB.numRow; i++ ){
                //                pos = std::find(pos, &rowsSinvBPtr[SinvB.numRow-1], rowsLBPtr[i]);
//                pos = std::find(pos/*rowsSinvBPtr*/, last, rowsLBPtr[LB.rowsPerm[i]]);

                pos = std::upper_bound(pos/*rowsSinvBPtr*/, last, rowsLBPtr[i])-1;
                if(pos != last){
                  relRows[i] =  (Int)(pos - rowsSinvBPtr);
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


              // Column relative indicies
              std::vector<Int> relCols( LrowB.numCol );
              Int SinvColsSta = FirstBlockCol( jsup, this->super_ );
              for( Int j = 0; j < LrowB.numCol; j++ ){
                relCols[j] = LrowB.rows[j] - SinvColsSta;
              }




              TIMER_START(Copy_Sinv_to_Ainv);
              // Transfer the values from Sinv to AinvBlock
              T* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < LrowB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
//                  statusOFS<< "nzvalAinv["<<i<<","<<j<<"] = nzvalSinv["<<relRows[i]<<","<<relCols[j]<<"] = "<< nzvalSinv[relRows[i] + relCols[j] * ldSinv]<< std::endl;
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }
              TIMER_STOP(Copy_Sinv_to_Ainv);




              isBlockFound = true;
              break;
            }
          } // for (ibSinv )
          TIMER_STOP(PARSING_ROW_BLOCKIDX);
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
          Int SinvRowsSta = FirstBlockCol( isup, this->super_ );
          for( Int i = 0; i < LB.numRow; i++ ){
            relRows[i] = LB.rows[i] - SinvRowsSta;
          }
          std::vector<UBlock<T> >&   UrowSinv = this->U( LBi( isup, this->grid_ ) );
          bool isBlockFound = false;
          TIMER_START(PARSING_COL_BLOCKIDX);
          for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
            // Found the (isup, jsup) block in Sinv
            if( UrowSinv[jbSinv].blockIdx == jsup ){
              UBlock<T>& SinvB = UrowSinv[jbSinv];

              // Column relative indices
              std::vector<Int> relCols( LrowB.numCol );
              Int* rowsLrowBPtr    = LrowB.rows.Data();
              Int* colsSinvBPtr = SinvB.cols.Data();

              TIMER_START(STDFIND_COL);
              Int * pos =&colsSinvBPtr[0];
              Int * last =&colsSinvBPtr[SinvB.numCol];
              for( Int j = 0; j < LrowB.numCol; j++ ){
                //rowsLrowB is sorted
                pos = std::upper_bound(pos, last, rowsLrowBPtr[j])-1;
                if(pos !=last){
                  relCols[j] = (Int)(pos - colsSinvBPtr);
                }
                else{
                  std::ostringstream msg;
                  msg << "Col " << rowsLrowBPtr[j] <<
                    " in UB cannot find the corresponding row in SinvB" << std::endl
                    << "LrowB.rows    = " << LrowB.rows << std::endl
                    << "UinvB.cols = " << SinvB.cols << std::endl;
                  throw std::runtime_error( msg.str().c_str() );
                }
              }
              TIMER_STOP(STDFIND_COL);


              TIMER_START(Copy_Sinv_to_Ainv);
              // Transfer the values from Sinv to AinvBlock
              T* nzvalSinv = SinvB.nzval.Data();
              Int     ldSinv    = SinvB.numRow;
              for( Int j = 0; j < LrowB.numCol; j++ ){
                for( Int i = 0; i < LB.numRow; i++ ){
//                  statusOFS<< "nzvalAinv["<<i<<","<<j<<"] = nzvalSinv["<<relRows[i]<<","<<relCols[j]<<"] = "<< nzvalSinv[relRows[i] + relCols[j] * ldSinv]<< std::endl;
                  nzvalAinv[i+j*ldAinv] =
                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                }
              }
              TIMER_STOP(Copy_Sinv_to_Ainv);

              isBlockFound = true;
              break;
            }
          } // for (jbSinv)
          TIMER_STOP(PARSING_COL_BLOCKIDX);
          if( isBlockFound == false ){
            std::ostringstream msg;
            msg << "Block(" << isup << ", " << jsup
              << ") did not find a matching block in Sinv." << std::endl;
            throw std::runtime_error( msg.str().c_str() );
          }
        } // if (isup, jsup) is in U

      } // for( ib )
    } // for ( jb )


      /*
      //    for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
      //
      //      UBlock& UB = UrowRecv[jb];
      //      Int jsup = UB.blockIdx;
      //      Int SinvColsSta = FirstBlockCol( jsup, this->super_ );
      //
      //      // Column relative indicies
      //      std::vector<Int> relCols( UB.numCol );
      //      for( Int j = 0; j < UB.numCol; j++ ){
      //        relCols[j] = UB.cols[j] - SinvColsSta;
      //      }
      //
      //
      //
      //
      //      for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
      //        LBlock& LB = LcolRecv[ib];
      //        Int isup = LB.blockIdx;
      //        Int SinvRowsSta = FirstBlockCol( isup, this->super_ );
      //        Scalar* nzvalAinv = &AinvBuf( rowPtr[ib], colPtr[jb] );
      //        Int     ldAinv    = numRowAinvBuf;
      //
      //        // Pin down the corresponding block in the part of Sinv.
      //        if( isup >= jsup ){
      //          std::vector<LBlock>&  LcolSinv = this->L( LBj(jsup, this->grid_ ) );
      //          bool isBlockFound = false;
      //          TIMER_START(PARSING_ROW_BLOCKIDX);
      //          for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
      //            // Found the (isup, jsup) block in Sinv
      //            if( LcolSinv[ibSinv].blockIdx == isup ){
      //              LBlock& SinvB = LcolSinv[ibSinv];
      //
      //              // Row relative indices
      //              std::vector<Int> relRows( LB.numRow );
      //              Int* rowsLBPtr    = LB.rows.Data();
      //              Int* rowsSinvBPtr = SinvB.rows.Data();
      //
      //              TIMER_START(STDFIND_ROW);
      //              Int * pos =&rowsSinvBPtr[0];
      //              Int * last =&rowsSinvBPtr[SinvB.numRow];
      //              for( Int i = 0; i < LB.numRow; i++ ){
      //                //                pos = std::find(pos, &rowsSinvBPtr[SinvB.numRow-1], rowsLBPtr[i]);
      //                pos = std::find(rowsSinvBPtr, last, rowsLBPtr[i]);
      //                if(pos != last){
      //                  relRows[i] = (Int)(pos - rowsSinvBPtr);
      //                }
      //                else{
      //                  std::ostringstream msg;
      //                  msg << "Row " << rowsLBPtr[i] << 
      //                    " in LB cannot find the corresponding row in SinvB" << std::endl
      //                    << "LB.rows    = " << LB.rows << std::endl
      //                    << "SinvB.rows = " << SinvB.rows << std::endl;
      //                  #ifdef USE_ABORT
      abort();
#endif
throw std::runtime_error( msg.str().c_str() );
      //                }
      //              }
      //              TIMER_STOP(STDFIND_ROW);
      //
      //              TIMER_START(Copy_Sinv_to_Ainv);
      //              // Transfer the values from Sinv to AinvBlock
      //              Scalar* nzvalSinv = SinvB.nzval.Data();
      //              Int     ldSinv    = SinvB.numRow;
      //              for( Int j = 0; j < UB.numCol; j++ ){
      //                for( Int i = 0; i < LB.numRow; i++ ){
      //                  nzvalAinv[i+j*ldAinv] =
      //                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
      //                }
      //              }
      //              TIMER_STOP(Copy_Sinv_to_Ainv);
      //
      //              isBlockFound = true;
      //              break;
      //            }	
      //          } // for (ibSinv )
      //          TIMER_STOP(PARSING_ROW_BLOCKIDX);
      //          if( isBlockFound == false ){
      //            std::ostringstream msg;
      //            msg << "Block(" << isup << ", " << jsup 
      //              << ") did not find a matching block in Sinv." << std::endl;
      //            #ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( msg.str().c_str() );
      //          }
      //        } // if (isup, jsup) is in L
      //        else{
      //          // Row relative indices
      //          std::vector<Int> relRows( LB.numRow );
      //          Int SinvRowsSta = FirstBlockCol( isup, this->super_ );
      //          for( Int i = 0; i < LB.numRow; i++ ){
      //            relRows[i] = LB.rows[i] - SinvRowsSta;
      //          }
      //          std::vector<UBlock>&   UrowSinv = this->U( LBi( isup, this->grid_ ) );
      //          bool isBlockFound = false;
      //          TIMER_START(PARSING_COL_BLOCKIDX);
      //          for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
      //            // Found the (isup, jsup) block in Sinv
      //            if( UrowSinv[jbSinv].blockIdx == jsup ){
      //              UBlock& SinvB = UrowSinv[jbSinv];
      //
      //
      //
      //              // Column relative indices
      //              std::vector<Int> relCols( UB.numCol );
      //              Int* colsUBPtr    = UB.cols.Data();
      //              Int* colsSinvBPtr = SinvB.cols.Data();
      //              TIMER_START(STDFIND_COL);
      //              Int * pos =&colsSinvBPtr[0];
      //              Int * last =&colsSinvBPtr[SinvB.numCol];
      //              for( Int j = 0; j < UB.numCol; j++ ){
      //                //colsUB is sorted
      //                pos = std::find(colsSinvBPtr, last, colsUBPtr[j]);
      //                if(pos !=last){
      //                  relCols[j] = (Int)(pos - colsSinvBPtr);
      //                }
      //                else{
      //                  std::ostringstream msg;
      //                  msg << "Col " << colsUBPtr[j] << 
      //                    " in UB cannot find the corresponding row in SinvB" << std::endl
      //                    << "UB.cols    = " << UB.cols << std::endl
      //                    << "UinvB.cols = " << SinvB.cols << std::endl;
      //                  #ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( msg.str().c_str() );
      //                }
      //              }
      //              TIMER_STOP(STDFIND_COL);
      //
      //
      //              TIMER_START(Copy_Sinv_to_Ainv);
      //              // Transfer the values from Sinv to AinvBlock
      //              Scalar* nzvalSinv = SinvB.nzval.Data();
      //              Int     ldSinv    = SinvB.numRow;
      //              for( Int j = 0; j < UB.numCol; j++ ){
      //                for( Int i = 0; i < LB.numRow; i++ ){
      //                  nzvalAinv[i+j*ldAinv] =
      //                    nzvalSinv[relRows[i] + relCols[j] * ldSinv];
      //                }
      //              }
      //              TIMER_STOP(Copy_Sinv_to_Ainv);
      //
      //              isBlockFound = true;
      //              break;
      //            }
      //          } // for (jbSinv)
      //          TIMER_STOP(PARSING_COL_BLOCKIDX);
      //          if( isBlockFound == false ){
      //            std::ostringstream msg;
      //            msg << "Block(" << isup << ", " << jsup 
      //              << ") did not find a matching block in Sinv." << std::endl;
      //            #ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( msg.str().c_str() );
      //          }
      //        } // if (isup, jsup) is in U
      //
      //      } // for( ib )
      //    } // for ( jb )
      */
#else
        for( Int jb = 0; jb < LrowRecv.size(); jb++ ){
          for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
            LBlock<T>& LB = LcolRecv[ib];
            LBlock<T>& LrowB = LrowRecv[jb];
            Int isup = LB.blockIdx;
            Int jsup = LrowB.blockIdx;
            T* nzvalAinv = &AinvBuf( rowPtrL[ib], colPtrL[jb] );
            Int     ldAinv    = AinvBuf.m();

            // Pin down the corresponding block in the part of Sinv.
            if( isup >= jsup ){
              std::vector<LBlock<T> >&  LcolSinv = this->L( LBj(jsup, this->grid_ ) );
              bool isBlockFound = false;
              TIMER_START(PARSING_ROW_BLOCKIDX);
              for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
                // Found the (isup, jsup) block in Sinv
                if( LcolSinv[ibSinv].blockIdx == isup ){
                  LBlock<T> & SinvB = LcolSinv[ibSinv];

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
#ifdef USE_ABORT
                      statusOFS<<msg.str();
                      abort();
#endif
                      throw std::runtime_error( msg.str().c_str() );
                    }
                  }

                  // Column relative indicies
                  std::vector<Int> relCols( LrowB.numCol );
                  Int SinvColsSta = FirstBlockCol( jsup, this->super_ );
                  for( Int j = 0; j < LrowB.numCol; j++ ){
                    relCols[j] = LrowB.rows[j] - SinvColsSta;
                  }

                  // Transfer the values from Sinv to AinvBlock
                  T* nzvalSinv = SinvB.nzval.Data();
                  Int     ldSinv    = SinvB.numRow;
                  for( Int j = 0; j < LrowB.numCol; j++ ){
                    for( Int i = 0; i < LB.numRow; i++ ){
                      nzvalAinv[i+j*ldAinv] =
                        nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                    }
                  }

                  isBlockFound = true;
                  break;
                }	
              } // for (ibSinv )
              TIMER_STOP(PARSING_ROW_BLOCKIDX);
              if( isBlockFound == false ){
                std::ostringstream msg;
                msg << "Block(" << isup << ", " << jsup 
                  << ") did not find a matching block in Sinv." << std::endl;
#ifdef USE_ABORT
                      statusOFS<<msg.str();
                abort();
#endif
                throw std::runtime_error( msg.str().c_str() );
              }
            } // if (isup, jsup) is in L
            else{
              //std::vector<LBlock<T> >& LrowSinv = this->Lrow( LBi( isup, this->grid_ ) );
              std::vector<UBlock<T> >& UrowSinv = this->U( LBi( isup, this->grid_ ) );
              bool isBlockFound = false;
              TIMER_START(PARSING_COL_BLOCKIDX);
              for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
                // Found the (isup, jsup) block in Sinv
                if( UrowSinv[jbSinv].blockIdx == jsup ){
                  UBlock<T> & SinvB = UrowSinv[jbSinv];

                  // Row relative indices
                  std::vector<Int> relRows( LB.numRow );
                  Int SinvRowsSta = FirstBlockCol( isup, this->super_ );
                  for( Int i = 0; i < LB.numRow; i++ ){
                    relRows[i] = LB.rows[i] - SinvRowsSta;
                  }

                  // Column relative indices
                  std::vector<Int> relCols( LrowB.numCol );
                  Int* rowsLrowBPtr    = LrowB.rows.Data();
                  Int* colsSinvBPtr = SinvB.cols.Data();
                  for( Int j = 0; j < LrowB.numCol; j++ ){
                    bool isColFound = false;
                    for( Int j1 = 0; j1 < SinvB.numCol; j1++ ){
                      if( rowsLrowBPtr[j] == colsSinvBPtr[j1] ){
                        isColFound = true;
                        relCols[j] = j1;
                        break;
                      }
                    }
                    if( isColFound == false ){
                      std::ostringstream msg;
                      msg << "Rows " << rowsLrowBPtr[j] << 
                        " in LrowB cannot find the corresponding row in SinvB" << std::endl
                        << "LrowB.rows    = " << LrowB.rows << std::endl
                        << "UinvB.cols = " << SinvB.cols << std::endl;
#ifdef USE_ABORT
                      statusOFS<<msg.str();
                      abort();
#endif
                      throw std::runtime_error( msg.str().c_str() );
                    }
                  }


                  // Transfer the values from Sinv to AinvBlock
                  T* nzvalSinv = SinvB.nzval.Data();
                  Int     ldSinv    = SinvB.numRow;
                  for( Int j = 0; j < LrowB.numCol; j++ ){
                    for( Int i = 0; i < LB.numRow; i++ ){
                      nzvalAinv[i+j*ldAinv] =
                        nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                    }
                  }

                  isBlockFound = true;
                  break;
                }
              } // for (jbSinv)
              TIMER_STOP(PARSING_COL_BLOCKIDX);
              if( !isBlockFound ){
                std::ostringstream msg;
                msg << "Block(" << isup << ", " << jsup 
                  << ") did not find a matching block in Sinv." << std::endl;
#ifdef USE_ABORT
                statusOFS<<msg.str();
                abort();
#endif
                throw std::runtime_error( msg.str().c_str() );
              }
            } // if (isup, jsup) is in U

          } // for( ib )
        } // for ( jb )
#endif
      TIMER_STOP(JB_Loop);

      TIMER_STOP(Compute_Sinv_LT_Lookup_Indexes);
    } // End of method PMatrixAsym::Selinv_lookup_indexes


  template<typename T>
    inline void PMatrixAsym<T>::SendRecvCD(
        std::vector<SuperNodeBufferTypeAsym > & arrSuperNodes, 
        Int stepSuper)
    {

      TIMER_START(Send_CD);
      //compute the number of requests
      Int sendCount = 0;
      Int recvCount = 0;
      Int sendOffset[stepSuper];
      Int recvOffset[stepSuper];
      Int recvIdxL=0;
      Int recvIdxU=0;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
        sendOffset[supidx]=sendCount;
        recvOffset[supidx]=recvCount;
        sendCount+= this->CountSendToCrossDiagonal(snode.Index);
        recvCount+= this->CountRecvFromCrossDiagonal(snode.Index);
      }

      //Buffers for L
      std::vector<MPI_Request > arrMpiReqsSendLCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeSendLCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsRecvLCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeRecvLCD(recvCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrLcolSendCD(sendCount);
      std::vector<int > arrSstrLcolSizeSendCD(sendCount);
      std::vector<std::vector<char> > arrSstrLrowRecvCD(recvCount);
      std::vector<int > arrSstrLrowSizeRecvCD(recvCount);

      //Buffers for U
      std::vector<MPI_Request > arrMpiReqsSendUCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeSendUCD(recvCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsRecvUCD(sendCount, MPI_REQUEST_NULL );
      std::vector<MPI_Request > arrMpiReqsSizeRecvUCD(sendCount, MPI_REQUEST_NULL );
      std::vector<std::vector<char> > arrSstrUrowSendCD(recvCount);
      std::vector<int > arrSstrUrowSizeSendCD(recvCount);
      std::vector<std::vector<char> > arrSstrUcolRecvCD(sendCount);
      std::vector<int > arrSstrUcolSizeRecvCD(sendCount);



      //Do Isend for size and content of L and Irecv for sizes of U 
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

        TIMER_START(Send_L_Recv_Size_U_CrossDiag);

        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) 
            && this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode.Index ) ){

          Int sendIdxL = 0;
          Int recvIdxU = 0;
          for(Int dstCol = 0; dstCol<this->grid_->numProcCol; dstCol++){
            if(this->isSendToCrossDiagonal_(dstCol,snode.Index) ){
              Int dest = PNUM(PROW(snode.Index,this->grid_),dstCol,this->grid_);

              if( MYPROC( this->grid_ ) != dest	){
                //Send size and content of L
                MPI_Request & mpiReqSizeSendL = 
                  arrMpiReqsSizeSendLCD[sendOffset[supidx]+sendIdxL];
                MPI_Request & mpiReqSendL = 
                  arrMpiReqsSendLCD[sendOffset[supidx]+sendIdxL];
                std::vector<char> & sstrLcolSend = 
                  arrSstrLcolSendCD[sendOffset[supidx]+sendIdxL];
                Int & sstrSizeL = 
                  arrSstrLcolSizeSendCD[sendOffset[supidx]+sendIdxL];

                std::stringstream sstm;
                serialize( snode.RowLocalPtr, sstm, NO_MASK );
                serialize( snode.BlockIdxLocal, sstm, NO_MASK );
                serialize( snode.LUpdateBuf, sstm, NO_MASK );
                sstrLcolSend.resize( Size(sstm) );
                sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
                sstrSizeL = sstrLcolSend.size();

                MPI_Isend( &sstrSizeL, 1, MPI_INT, dest, 
                    IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE),
                    this->grid_->comm, &mpiReqSizeSendL );
                MPI_Isend( (void*)&sstrLcolSend[0], sstrSizeL, MPI_BYTE, dest,
                    IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT),
                    this->grid_->comm, &mpiReqSendL );
                sendIdxL++;

                //Recv for U size
                Int & sstrSizeU = 
                  arrSstrUcolSizeRecvCD[sendOffset[supidx]+recvIdxU];
                MPI_Request & mpiReqSizeRecvU = 
                  arrMpiReqsSizeRecvUCD[sendOffset[supidx]+recvIdxU];

                MPI_Irecv( &sstrSizeU, 1, MPI_INT, dest, 
                    IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE),
                    this->grid_->comm, &mpiReqSizeRecvU );
                recvIdxU++;
              }
            }
          }
        } // sender
        TIMER_STOP(Send_L_Recv_Size_U_CrossDiag);
      }


      //Do Irecv for sizes of L and Isend for size and content of U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
        //If I'm a receiver
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) 
            && this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode.Index ) ){
          Int recvIdxL=0;
          Int sendIdxU=0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,this->grid_),this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                //Recv size of L
                Int & sstrSizeL = 
                  arrSstrLrowSizeRecvCD[recvOffset[supidx]+recvIdxL];
                MPI_Request & mpiReqSizeRecvL = 
                  arrMpiReqsSizeRecvLCD[recvOffset[supidx]+recvIdxL];

                MPI_Irecv( &sstrSizeL, 1, MPI_INT, src, 
                    IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE),
                    this->grid_->comm, &mpiReqSizeRecvL );
                recvIdxL++;


                //Send size and content of U
                MPI_Request & mpiReqSizeSendU = 
                  arrMpiReqsSizeSendUCD[recvOffset[supidx]+sendIdxU];
                MPI_Request & mpiReqSendU = 
                  arrMpiReqsSendUCD[recvOffset[supidx]+sendIdxU];
                std::vector<char> & sstrUrowSend = 
                  arrSstrUrowSendCD[recvOffset[supidx]+sendIdxU];
                Int & sstrSizeU = 
                  arrSstrUrowSizeSendCD[recvOffset[supidx]+sendIdxU];

                std::stringstream sstm;
                serialize( snode.ColLocalPtr, sstm, NO_MASK );
                serialize( snode.BlockIdxLocalU, sstm, NO_MASK );
                serialize( snode.UUpdateBuf, sstm, NO_MASK );
                sstrUrowSend.resize( Size(sstm) );
                sstm.read( &sstrUrowSend[0], sstrUrowSend.size() );
                sstrSizeU = sstrUrowSend.size();

                MPI_Isend( &sstrSizeU, 1, MPI_INT, src, 
                    IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE), 
                    this->grid_->comm, &mpiReqSizeSendU );
                MPI_Isend( (void*)&sstrUrowSend[0], sstrSizeU, MPI_BYTE, src,
                    IDX_TO_TAG(supidx,SELINV_TAG_U_CONTENT),
                    this->grid_->comm, &mpiReqSendU );
                sendIdxU++;
              }
            }
          }
        }//end if I'm a receiver
      }

      //waitall sizes of L
      mpi::Waitall(arrMpiReqsSizeRecvLCD);

      //Allocate content and do Irecv for content of L
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
        //If I'm a receiver
        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) 
            && this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode.Index ) ){
          Int recvIdxL=0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,this->grid_),this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                Int & sstrSizeL = 
                  arrSstrLrowSizeRecvCD[recvOffset[supidx]+recvIdxL];
                std::vector<char> & sstrLrowRecv = 
                  arrSstrLrowRecvCD[recvOffset[supidx]+recvIdxL];
                MPI_Request & mpiReqRecvL = 
                  arrMpiReqsRecvLCD[recvOffset[supidx]+recvIdxL];
                sstrLrowRecv.resize( sstrSizeL);

                MPI_Irecv( (void*)&sstrLrowRecv[0], sstrSizeL, MPI_BYTE, src,
                    IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT),
                    this->grid_->comm, &mpiReqRecvL );
                recvIdxL++;
              }
            }
          }
        }//end if I'm a receiver
      }

      //waitall sizes of U
      mpi::Waitall(arrMpiReqsSizeRecvUCD);

      //Allocate content and do Irecv for content of U
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) 
            && this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode.Index ) ){
          Int recvIdxU = 0;
          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,snode.Index) ){
              Int src = PNUM(PROW(snode.Index,this->grid_),srcCol,this->grid_);
              if( MYPROC( this->grid_ ) != src ){
                Int & sstrSizeU = 
                  arrSstrUcolSizeRecvCD[sendOffset[supidx]+recvIdxU];
                std::vector<char> & sstrUcolRecv = 
                  arrSstrUcolRecvCD[sendOffset[supidx]+recvIdxU];
                MPI_Request & mpiReqRecvU = 
                  arrMpiReqsRecvUCD[sendOffset[supidx]+recvIdxU];
                sstrUcolRecv.resize( sstrSizeU );

                MPI_Irecv( (void*)&sstrUcolRecv[0], sstrSizeU, MPI_BYTE, src,
                    IDX_TO_TAG(supidx,SELINV_TAG_U_CONTENT),
                    this->grid_->comm, &mpiReqRecvU );
                recvIdxU++;
              }
            }
          }
        }//end if I'm a receiver
      }

      //waitall content of L and U
      mpi::Waitall(arrMpiReqsRecvLCD);
      mpi::Waitall(arrMpiReqsRecvUCD);







      //Do the work
      NumMat<T> Ltmp;
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) &&
            this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, snode.Index ) ){

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << " ["<<snode.Index<<"] "
            << "Update the upper triangular block" 
            << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode.Index<<"] "
            << "blockIdxLocal:" << snode.BlockIdxLocal
            << std::endl << std::endl; 
          statusOFS << std::endl << " ["<<snode.Index<<"] "
            << "rowLocalPtr:" << snode.RowLocalPtr
            << std::endl << std::endl; 
#endif

          std::vector<LBlock<T> >& Lrow = this->Lrow( LBi( snode.Index, this->grid_ ) );
          std::vector<Int> isBlockFound(Lrow.size(),false);
          Int recvIdxL=0;

          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,snode.Index) ){
              Int src = PNUM(srcRow,PCOL(snode.Index,this->grid_),this->grid_);
              TIMER_START(Recv_L_CrossDiag);

              //Declaring some pointers to avoid copy if data is local
              std::vector<Int> *pRowLocalPtr;
              std::vector<Int> *pBlockIdxLocal;
              NumMat<T> *pLUpdateBuf;

              std::vector<Int> rowLocalPtrRecv;
              std::vector<Int> blockIdxLocalRecv;
              NumMat<T> LUpdateBuf;

              if( MYPROC( this->grid_ ) != src ){
                std::stringstream sstm;
                Int & sstrSizeL = 
                  arrSstrLrowSizeRecvCD[recvOffset[supidx]+recvIdxL];
                std::vector<char> & sstrLrowRecv = 
                  arrSstrLrowRecvCD[recvOffset[supidx]+recvIdxL];

                sstm.write( &sstrLrowRecv[0], sstrSizeL );
                deserialize( rowLocalPtrRecv, sstm, NO_MASK );
                deserialize( blockIdxLocalRecv, sstm, NO_MASK );
                deserialize( LUpdateBuf, sstm, NO_MASK );
                pRowLocalPtr = &rowLocalPtrRecv;
                pBlockIdxLocal = &blockIdxLocalRecv;
                pLUpdateBuf = &LUpdateBuf;
                recvIdxL++;
              } // sender is not the same as receiver
              else{
                pRowLocalPtr = &snode.RowLocalPtr;
                pBlockIdxLocal = &snode.BlockIdxLocal;
                pLUpdateBuf = &snode.LUpdateBuf;
              } // sender is the same as receiver

              TIMER_STOP(Recv_L_CrossDiag);

#if ( _DEBUGlevel_ >= 1 )
              statusOFS <<" ["<<snode.Index<<"] P"<<MYPROC(this->grid_)<<" ("
                <<MYROW(this->grid_) <<","<<MYCOL(this->grid_)<<") <--- LBj("
                <<snode.Index<<") <--- P"<<src<<std::endl;
              statusOFS << std::endl << " ["<<snode.Index<<"] "
                << "rowLocalPtrRecv:" << pRowLocalPtr
                << std::endl << std::endl; 
              statusOFS << std::endl << " ["<<snode.Index<<"] "
                << "blockIdxLocalRecv:" << pBlockIdxLocal
                << std::endl << std::endl; 
#endif

              // Update Lrow
              for( Int ib = 0; ib < pBlockIdxLocal->size(); ib++ ){
                for( Int iib = 0; iib < Lrow.size(); iib++ ){
                  LBlock<T> &  LrowB = Lrow[iib];
                  if( LrowB.blockIdx == (*pBlockIdxLocal)[ib] ){
                    
                    Ltmp.Resize( LrowB.numCol, LrowB.numRow );
                    lapack::Lacpy( 'A', Ltmp.m(), Ltmp.n(), 
                        &(*pLUpdateBuf)( (*pRowLocalPtr)[ib], 0 ),
                        pLUpdateBuf->m(), Ltmp.Data(), Ltmp.m() );
                    isBlockFound[ib] = 1;
                    Transpose( Ltmp, LrowB.nzval );
                    break;
                  }
                }
              }
            }
          }

          for( Int ib = 0; ib < Lrow.size(); ib++ ){
            if( !isBlockFound[ib] ){
#ifdef USE_ABORT
              abort();
#endif
              throw std::logic_error( 
                  "LBlock cannot find its update. Something is seriously wrong."
                  );
            }
          }
        } // Done copying L


        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) &&
            this->isSendToCrossDiagonal_(this->grid_->numProcCol, snode.Index ) ){

          std::vector<UBlock<T> >& Ucol = this->Ucol( LBj( snode.Index, this->grid_ ) );
          std::vector<Int> isBlockFound(Ucol.size(),false);
          Int recvIdxU=0;

          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,snode.Index) ){
              Int src = PNUM(PROW(snode.Index,this->grid_),srcCol,this->grid_);
              TIMER_START(Recv_U_CrossDiag);

              //Declaring some pointers to avoid copy if data is local
              std::vector<Int> *pColLocalPtr;
              std::vector<Int> *pBlockIdxLocal;
              NumMat<T> *pUUpdateBuf;

              std::vector<Int> colLocalPtrRecv;
              std::vector<Int> blockIdxLocalRecv;
              NumMat<T> UUpdateBuf;

              if( MYPROC( this->grid_ ) != src ){
                std::stringstream sstm;
                Int & sstrSizeU = 
                  arrSstrUcolSizeRecvCD[sendOffset[supidx]+recvIdxU];
                std::vector<char> & sstrUcolRecv = 
                  arrSstrUcolRecvCD[sendOffset[supidx]+recvIdxU];

                sstm.write( &sstrUcolRecv[0], sstrSizeU );
                deserialize( colLocalPtrRecv, sstm, NO_MASK );
                deserialize( blockIdxLocalRecv, sstm, NO_MASK );
                deserialize( UUpdateBuf, sstm, NO_MASK );
                pColLocalPtr = &colLocalPtrRecv;
                pBlockIdxLocal = &blockIdxLocalRecv;
                pUUpdateBuf = &UUpdateBuf;
                recvIdxU++;
              } // sender is not the same as receiver
              else{
                pColLocalPtr = &snode.ColLocalPtr;
                pBlockIdxLocal = &snode.BlockIdxLocalU;
                pUUpdateBuf = &snode.UUpdateBuf;
              } // sender is the same as receiver

              TIMER_STOP(Recv_U_CrossDiag);

#if ( _DEBUGlevel_ >= 1 )
              statusOFS <<" ["<<snode.Index<<"] P"<<MYPROC(this->grid_)<<" ("
                <<MYROW(this->grid_) <<","<<MYCOL(this->grid_)<<") <--- LBi("
                <<snode.Index<<") <--- P"<<src<<std::endl;
              statusOFS << std::endl << " ["<<snode.Index<<"] "
                << "colLocalPtrRecv:" << pColLocalPtr
                << std::endl << std::endl; 
              statusOFS << std::endl << " ["<<snode.Index<<"] "
                << "blockIdxLocalRecv:" << pBlockIdxLocal
                << std::endl << std::endl; 
#endif

              // Update Ucol
              for( Int jb = 0; jb < pBlockIdxLocal->size(); jb++ ){
                for( Int jjb = 0; jjb < Ucol.size(); jjb++ ){
                  UBlock<T> &  UcolB = Ucol[jjb];
                  if( UcolB.blockIdx == (*pBlockIdxLocal)[jb] ){
                    lapack::Lacpy( 'A', UcolB.nzval.m(), UcolB.nzval.n(), 
                        &(*pUUpdateBuf)( 0, (*pColLocalPtr)[jb] ),
                        pUUpdateBuf->m(), UcolB.nzval.Data(), UcolB.nzval.m() );
                    isBlockFound[jb] = true;
                    break;
                  }
                }
              }
            }
          }

          for( Int jb = 0; jb < Ucol.size(); jb++ ){
            if( !isBlockFound[jb] ){
#ifdef USE_ABORT
              abort();
#endif
              throw std::logic_error( 
                  "UBlock cannot find its update. Something is seriously wrong."
                  );
            }
          } // end for (jb)
        } // Done copying U

      }

      TIMER_STOP(Send_CD);

      mpi::Waitall(arrMpiReqsSizeSendLCD);
      mpi::Waitall(arrMpiReqsSendLCD);

      mpi::Waitall(arrMpiReqsSizeSendUCD);
      mpi::Waitall(arrMpiReqsSendUCD);
    } // End of method PMatrixAsym::SendRecvCD 

  template<typename T>
    inline void PMatrixAsym<T>::UnpackData(
        SuperNodeBufferTypeAsym & snode, 
        std::vector<LBlock<T> > & LcolRecv, 
        std::vector<LBlock<T> > & LrowRecv,
        std::vector<UBlock<T> > & UcolRecv, 
        std::vector<UBlock<T> > & UrowRecv
        )
    {

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Unpack the received data for processors participate in Gemm. " << std::endl << std::endl; 
#endif
      // Lrow part
      if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
        std::stringstream sstm;
        sstm.write( &snode.SstrLrowRecv[0], snode.SstrLrowRecv.size() );
        std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
        Int numLBlock;
        deserialize( numLBlock, sstm, NO_MASK );
        LrowRecv.resize( numLBlock );
        for( Int ib = 0; ib < numLBlock; ib++ ){
          deserialize( LrowRecv[ib], sstm, mask );
        } 
      } // sender is not the same as receiver
      else{
        // U is obtained locally, just make a copy. Include everything
        // (there is no diagonal block)
        // Is it a copy ?  LL: YES. Maybe we should replace the copy by
        // something more efficient especially for mpisize == 1
        LrowRecv.resize(this->Lrow( LBi( snode.Index, this->grid_ ) ).size());
        std::copy(this->Lrow( LBi( snode.Index, this->grid_ ) ).begin(),this->Lrow( LBi( snode.Index, this->grid_ )).end(),LrowRecv.begin());
      } // sender is the same as receiver


      //L part
      if( MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
        std::stringstream     sstm;
        sstm.write( &snode.SstrLcolRecv[0], snode.SstrLcolRecv.size() );
        std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
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
        std::vector<LBlock<T> >& Lcol =  this->L( LBj( snode.Index, this->grid_ ) );
        Int startIdx = ( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) )?1:0;
        LcolRecv.resize( Lcol.size() - startIdx );
        std::copy(Lcol.begin()+startIdx,Lcol.end(),LcolRecv.begin());
      } // sender is the same as receiver

      // U part
      if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
        std::stringstream sstm;
        sstm.write( &snode.SstrUrowRecv[0], snode.SstrUrowRecv.size() );
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
        // Is it a copy ?  LL: YES. Maybe we should replace the copy by
        // something more efficient especially for mpisize == 1
        UrowRecv.resize(this->U( LBi( snode.Index, this->grid_ ) ).size());
        std::copy(this->U( LBi( snode.Index, this->grid_ ) ).begin(),this->U( LBi( snode.Index, this->grid_ )).end(),UrowRecv.begin());
      } // sender is the same as receiver


      // Ucol part
      if( MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
        std::stringstream     sstm;
        sstm.write( &snode.SstrUcolRecv[0], snode.SstrUcolRecv.size() );
        std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
        Int numUBlock;
        deserialize( numUBlock, sstm, NO_MASK );
        UcolRecv.resize( numUBlock );
        for( Int jb = 0; jb < numUBlock; jb++ ){
          deserialize( UcolRecv[jb], sstm, mask );
        }
      } // sender is not the same as receiver
      else{
        // Ucol is obtained locally, just make a copy. 
        // Do not include the diagonal block
        std::vector<UBlock<T> >& Ucol =  this->Ucol( LBj( snode.Index, this->grid_ ) );
        UcolRecv.resize( Ucol.size() );
        std::copy(Ucol.begin(),Ucol.end(),UcolRecv.begin());
      } // sender is the same as receiver
    } // End of method PMatrixAsym<T>::UnpackData

  template<typename T>
    inline void PMatrixAsym<T>::ComputeDiagUpdate(SuperNodeBufferTypeAsym & snode)
    {

      //--------- Computing  Diagonal block, all processors in the column
      //--------- are participating to all pipelined supernodes
      if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Updating the diagonal block" << std::endl << std::endl;
#endif
        std::vector<LBlock<T> >&  Lcol = this->L( LBj( snode.Index, this->grid_ ) );
        std::vector<UBlock<T> >&  Ucol = this->Ucol( LBj( snode.Index, this->grid_ ) );

        //Allocate DiagBuf even if Lcol.size() == 0
        snode.DiagBuf.Resize(SuperSize( snode.Index, this->super_ ),
            SuperSize( snode.Index, this->super_ ));
        SetValue(snode.DiagBuf, ZERO<T>());

        // Do I own the diagonal block ?
        Int startIb = (MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ))?1:0;
        for( Int ib = startIb; ib < Lcol.size(); ib++ ){
          LBlock<T> & LcolB = Lcol[ib];
          //find corresponding U block
          Int jb = 0;
          for(jb = 0; jb < Ucol.size(); ++jb){
            if(Ucol[jb].blockIdx == LcolB.blockIdx){
              break;
            }
          }
        
          assert(jb < Ucol.size());

          UBlock<T> & UcolB = Ucol[jb];

          //Compute U S-1 L
          blas::Gemm( 'N', 'N', snode.DiagBuf.m(), snode.DiagBuf.n(), 
              LcolB.numRow, MINUS_ONE<T>(),
              UcolB.nzval.Data(), UcolB.nzval.m(), 
              &snode.LUpdateBuf( snode.RowLocalPtr[ib-startIb], 0 ),
              snode.LUpdateBuf.m(), 
              ONE<T>(), snode.DiagBuf.Data(), snode.DiagBuf.m() );
        } 

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Updated the diagonal block" << std::endl << std::endl; 
#endif
      }
    } // End of method PMatrixAsym<T>::ComputeDiagUpdate 


  template<typename T>
    inline void PMatrixAsym<T>::SelInvIntra_P2p(Int lidx)
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
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendLToBelow;
      arrMpireqsSendLToBelow.resize( stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcRow, MPI_REQUEST_NULL ));
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendLToRight;
      arrMpireqsSendLToRight.resize(stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcCol, MPI_REQUEST_NULL ));

      std::vector<std::vector<MPI_Request> >  arrMpireqsSendUToBelow;
      arrMpireqsSendUToBelow.resize( stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcRow, MPI_REQUEST_NULL ));
      std::vector<std::vector<MPI_Request> >  arrMpireqsSendUToRight;
      arrMpireqsSendUToRight.resize(stepSuper, std::vector<MPI_Request>( 2 * this->grid_->numProcCol, MPI_REQUEST_NULL ));

      //This is required to reduce L
      std::vector<MPI_Request>  arrMpireqsSendLToLeft;
      arrMpireqsSendLToLeft.resize(stepSuper, MPI_REQUEST_NULL );

      //This is required to reduce U
      std::vector<MPI_Request>  arrMpireqsSendUToAbove;
      arrMpireqsSendUToAbove.resize(stepSuper, MPI_REQUEST_NULL );

      //This is required to reduce D
      std::vector<MPI_Request>  arrMpireqsSendToAbove;
      arrMpireqsSendToAbove.resize(stepSuper, MPI_REQUEST_NULL );

      //------------------------------------------------------------------------

      //This is required to receive the size and content of U/L
      std::vector<MPI_Request>   arrMpireqsRecvLSizeFromAny;
      arrMpireqsRecvLSizeFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);
      std::vector<MPI_Request>   arrMpireqsRecvLContentFromAny;
      arrMpireqsRecvLContentFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);

      std::vector<MPI_Request>   arrMpireqsRecvUSizeFromAny;
      arrMpireqsRecvUSizeFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);
      std::vector<MPI_Request>   arrMpireqsRecvUContentFromAny;
      arrMpireqsRecvUContentFromAny.resize(stepSuper*2 , MPI_REQUEST_NULL);

      //allocate the buffers for this supernode
      std::vector<SuperNodeBufferTypeAsym> arrSuperNodes(stepSuper);
      for (Int supidx=0; supidx<stepSuper; supidx++){ 
        arrSuperNodes[supidx].Index = superList[lidx][supidx];  
      }

      NumMat<T> AinvBuf, UBuf, LBuf;

      TIMER_STOP(AllocateBuffer);

#ifndef _RELEASE_
      PushCallStack("PMatrix::SelInv_P2p::UpdateLU");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Communication to the Schur complement." << std::endl << std::endl; 
#endif

      {
        // Senders
        for (Int supidx=0; supidx<stepSuper; supidx++){
          SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
          std::vector<MPI_Request> & mpireqsSendLToBelow = arrMpireqsSendLToBelow[supidx];
          std::vector<MPI_Request> & mpireqsSendLToRight = arrMpireqsSendLToRight[supidx];
          std::vector<MPI_Request> & mpireqsSendUToBelow = arrMpireqsSendUToBelow[supidx];
          std::vector<MPI_Request> & mpireqsSendUToRight = arrMpireqsSendUToRight[supidx];

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl <<  "["<<snode.Index<<"] "
            << "Communication for the Lrow part." << std::endl << std::endl; 
#endif
          // Communication for the Lrow part.
          if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
            // Pack the data in Lrow
            TIMER_START(Serialize_LrowL);
            std::stringstream sstm;
            std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
            std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi(snode.Index, this->grid_) );
            // All blocks are to be sent down.
            serialize( (Int)Lrow.size(), sstm, NO_MASK );
            for( Int ib = 0; ib < Lrow.size(); ib++ ){
              serialize( Lrow[ib], sstm, mask );
            }
            snode.SstrLrowSend.resize( Size( sstm ) );
            sstm.read( &snode.SstrLrowSend[0], snode.SstrLrowSend.size() );
            snode.SizeSstrLrowSend = snode.SstrLrowSend.size();
            TIMER_STOP(Serialize_LrowL);

            for( Int iProcRow = 0; iProcRow < this->grid_->numProcRow; iProcRow++ ){
              if( MYROW( this->grid_ ) != iProcRow &&
                  this->isSendToBelow_( iProcRow,snode.Index ) == true ){
                // Use Isend to send to multiple targets
                MPI_Isend( &snode.SizeSstrLrowSend, 1, MPI_INT,  
                    iProcRow, IDX_TO_TAG(supidx,SELINV_TAG_LROW_SIZE), this->grid_->colComm, &mpireqsSendLToBelow[2*iProcRow] );
                MPI_Isend( (void*)&snode.SstrLrowSend[0], snode.SizeSstrLrowSend, MPI_BYTE, 
                    iProcRow, IDX_TO_TAG(supidx,SELINV_TAG_LROW_CONTENT), 
                    this->grid_->colComm, &mpireqsSendLToBelow[2*iProcRow+1] );
#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending Lrow " << snode.SizeSstrLrowSend << " BYTES"<< std::endl <<  std::endl; 
#endif
              } // Send 
            } // for (iProcRow)
          } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<< "Communication for the L part." << std::endl << std::endl; 
#endif
          // Communication for the L (Lcol) part.
          if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){

            TIMER_START(Serialize_LrowL);
            // Pack the data in L 
            std::stringstream sstm;
            std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
            mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

            std::vector<LBlock<T> >&  Lcol = this->L( LBj(snode.Index, this->grid_) );
            // All blocks except for the diagonal block are to be sent right

            if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) )
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
            TIMER_STOP(Serialize_LrowL);

            for( Int iProcCol = 0; iProcCol < this->grid_->numProcCol ; iProcCol++ ){
              if( MYCOL( this->grid_ ) != iProcCol &&
                  this->isSendToRight_( iProcCol, snode.Index ) == true ){
                // Use Isend to send to multiple targets
                MPI_Isend( &snode.SizeSstrLcolSend, 1, MPI_INT,  
                    iProcCol, IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE), 
                    this->grid_->rowComm, &mpireqsSendLToRight[2*iProcCol] );
                MPI_Isend( (void*)&snode.SstrLcolSend[0], snode.SizeSstrLcolSend, MPI_BYTE, 
                    iProcCol, IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT), 
                    this->grid_->rowComm, &mpireqsSendLToRight[2*iProcCol+1] );
#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending L " << snode.SizeSstrLcolSend<< " BYTES"  << std::endl <<  std::endl; 
#endif
              } // Send 
            } // for (iProcCol)
          } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl <<  "["<<snode.Index<<"] "
            << "Communication for the U part." << std::endl << std::endl; 
#endif
          // Communication for the U (Urow) part.
          if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){

            TIMER_START(Serialize_UcolU);
            // Pack the data in U 
            std::stringstream sstm;
            std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
            mask[UBlockMask::NZVAL] = 0; // nzval is excluded 

            std::vector<UBlock<T> >&  Urow = this->U( LBi(snode.Index, this->grid_) );
            // All blocks except for the diagonal block are to be sent right

            if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) )
              serialize( (Int)Urow.size() - 1, sstm, NO_MASK );
            else
              serialize( (Int)Urow.size(), sstm, NO_MASK );

            for( Int jb = 0; jb < Urow.size(); jb++ ){
              if( Urow[jb].blockIdx > snode.Index ){
#if ( _DEBUGlevel_ >= 2 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Serializing Block index " << Urow[jb].blockIdx << std::endl;
#endif
                serialize( Urow[jb], sstm, mask );
              }
            }
            snode.SstrUrowSend.resize( Size( sstm ) );
            sstm.read( &snode.SstrUrowSend[0], snode.SstrUrowSend.size() );
            snode.SizeSstrUrowSend = snode.SstrUrowSend.size();
            TIMER_STOP(Serialize_UcolU);

            for( Int iProcRow = 0; iProcRow < this->grid_->numProcRow; iProcRow++ ){
              if( MYROW( this->grid_ ) != iProcRow &&
                  this->isSendToBelow_( iProcRow,snode.Index ) == true ){
                // Use Isend to send to multiple targets
                MPI_Isend( &snode.SizeSstrUrowSend, 1, MPI_INT,  
                    iProcRow, IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE), this->grid_->colComm, &mpireqsSendUToBelow[2*iProcRow] );
                MPI_Isend( (void*)&snode.SstrLrowSend[0], snode.SizeSstrUrowSend, MPI_BYTE, 
                    iProcRow, IDX_TO_TAG(supidx,SELINV_TAG_U_CONTENT), 
                    this->grid_->colComm, &mpireqsSendUToBelow[2*iProcRow+1] );
#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending U " << snode.SizeSstrUrowSend << " BYTES"<< std::endl <<  std::endl; 
#endif
              } // Send 
            } // for (iProcRow)
          } // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << std::endl << "["<<snode.Index<<"] "<< "Communication for the Ucol part." << std::endl << std::endl; 
#endif
          // Communication for the Ucol part.
          if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
            // Pack the data in Ucol
            TIMER_START(Serialize_UcolU);
            std::stringstream sstm;
            std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
            std::vector<UBlock<T> >&  Ucol = this->Ucol( LBj(snode.Index, this->grid_) );
            // All blocks are to be sent down.
            serialize( (Int)Ucol.size(), sstm, NO_MASK );
            for( Int jb = 0; jb < Ucol.size(); jb++ ){
              serialize( Ucol[jb], sstm, mask );
            }
            snode.SstrUcolSend.resize( Size( sstm ) );
            sstm.read( &snode.SstrUcolSend[0], snode.SstrUcolSend.size() );
            snode.SizeSstrUcolSend = snode.SstrUcolSend.size();
            TIMER_STOP(Serialize_UcolU);

            for( Int iProcCol = 0; iProcCol < this->grid_->numProcCol ; iProcCol++ ){
              if( MYCOL( this->grid_ ) != iProcCol &&
                  this->isSendToRight_( iProcCol, snode.Index ) == true ){
                // Use Isend to send to multiple targets
                MPI_Isend( &snode.SizeSstrUcolSend, 1, MPI_INT,  
                    iProcCol, IDX_TO_TAG(supidx,SELINV_TAG_UCOL_SIZE), 
                    this->grid_->rowComm, &mpireqsSendUToRight[2*iProcCol] );
                MPI_Isend( (void*)&snode.SstrUcolSend[0], snode.SizeSstrUcolSend, MPI_BYTE, 
                    iProcCol, IDX_TO_TAG(supidx,SELINV_TAG_UCOL_CONTENT), 
                    this->grid_->rowComm, &mpireqsSendUToRight[2*iProcCol+1] );
#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Sending Ucol " << snode.SizeSstrUcolSend<< " BYTES"  << std::endl <<  std::endl; 
#endif
              } // Send 
            } // for (iProcCol)
          } // if I am the sender



        } //Senders

        //TODO Ideally, we should not receive data in sequence but in any order with ksup packed with the data
        // Receivers (Size)
        for (Int supidx=0; supidx<stepSuper ; supidx++){
          SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
          MPI_Request * mpireqsRecvLFromAbove = &arrMpireqsRecvLSizeFromAny[supidx*2];
          MPI_Request * mpireqsRecvLFromLeft = &arrMpireqsRecvLSizeFromAny[supidx*2+1];
          MPI_Request * mpireqsRecvUFromAbove = &arrMpireqsRecvUSizeFromAny[supidx*2];
          MPI_Request * mpireqsRecvUFromLeft = &arrMpireqsRecvUSizeFromAny[supidx*2+1];

          // Receive the size first
          if( this->isRecvFromAbove_( snode.Index ) && 
              MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
            MPI_Irecv( &snode.SizeSstrLrowRecv, 1, MPI_INT, PROW( snode.Index, this->grid_ ), 
                IDX_TO_TAG(supidx,SELINV_TAG_LROW_SIZE),
                this->grid_->colComm, mpireqsRecvLFromAbove );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving Lrow size on tag " << IDX_TO_TAG(supidx,SELINV_TAG_LROW_SIZE)<< std::endl <<  std::endl; 
#endif
            MPI_Irecv( &snode.SizeSstrUrowRecv, 1, MPI_INT, PROW( snode.Index, this->grid_ ), 
                IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE),
                this->grid_->colComm, mpireqsRecvUFromAbove );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving U size on tag " << IDX_TO_TAG(supidx,SELINV_TAG_U_SIZE)<< std::endl <<  std::endl; 
#endif
          } // if I need to receive from up


          if( this->isRecvFromLeft_( snode.Index ) &&
              MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
            MPI_Irecv( &snode.SizeSstrLcolRecv, 1, MPI_INT, PCOL( snode.Index, this->grid_ ), 
                IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE),
                this->grid_->rowComm, mpireqsRecvLFromLeft );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving L size on tag " << IDX_TO_TAG(supidx,SELINV_TAG_L_SIZE)<< std::endl <<  std::endl; 
#endif
            MPI_Irecv( &snode.SizeSstrUcolRecv, 1, MPI_INT, PCOL( snode.Index, this->grid_ ), 
                IDX_TO_TAG(supidx,SELINV_TAG_UCOL_SIZE),
                this->grid_->rowComm, mpireqsRecvUFromLeft );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving Ucol size on tag " << IDX_TO_TAG(supidx,SELINV_TAG_UCOL_SIZE)<< std::endl <<  std::endl; 
#endif
          } // if I need to receive from left
        }

        //Wait to receive all the sizes for L
        TIMER_START(WaitSize_LrowL);
        mpi::Waitall(arrMpireqsRecvLSizeFromAny);
        TIMER_STOP(WaitSize_LrowL);

        // Receivers (Content)
        for (Int supidx=0; supidx<stepSuper ; supidx++){
          SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

          MPI_Request * mpireqsRecvFromAbove = &arrMpireqsRecvLContentFromAny[supidx*2];
          MPI_Request * mpireqsRecvFromLeft = &arrMpireqsRecvLContentFromAny[supidx*2+1];

          if( this->isRecvFromAbove_( snode.Index ) && 
              MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
            snode.SstrLrowRecv.resize( snode.SizeSstrLrowRecv );
            MPI_Irecv( &snode.SstrLrowRecv[0], snode.SizeSstrLrowRecv, MPI_BYTE, 
                PROW( snode.Index, this->grid_ ), IDX_TO_TAG(supidx,SELINV_TAG_LROW_CONTENT), 
                this->grid_->colComm, mpireqsRecvFromAbove );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving Lrow " << snode.SizeSstrLrowRecv << " BYTES"<< std::endl <<  std::endl; 
#endif
          } // if I need to receive from up

          if( this->isRecvFromLeft_( snode.Index ) &&
              MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
            snode.SstrLcolRecv.resize( snode.SizeSstrLcolRecv );
            MPI_Irecv( &snode.SstrLcolRecv[0], snode.SizeSstrLcolRecv, MPI_BYTE, 
                PCOL( snode.Index, this->grid_ ), IDX_TO_TAG(supidx,SELINV_TAG_L_CONTENT), 
                this->grid_->rowComm,
                mpireqsRecvFromLeft );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving L " << snode.SizeSstrLcolRecv << " BYTES"<< std::endl <<  std::endl; 
#endif
          } // if I need to receive from left
        }


        //Wait to receive all the sizes for U
        TIMER_START(WaitSize_UcolU);
        mpi::Waitall(arrMpireqsRecvUSizeFromAny);
        TIMER_STOP(WaitSize_UcolU);

        // Receivers (Content)
        for (Int supidx=0; supidx<stepSuper ; supidx++){
          SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

          MPI_Request * mpireqsRecvFromAbove = &arrMpireqsRecvUContentFromAny[supidx*2];
          MPI_Request * mpireqsRecvFromLeft = &arrMpireqsRecvUContentFromAny[supidx*2+1];

          if( this->isRecvFromAbove_( snode.Index ) && 
              MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
            snode.SstrUrowRecv.resize( snode.SizeSstrUrowRecv );
            MPI_Irecv( &snode.SstrUrowRecv[0], snode.SizeSstrUrowRecv, MPI_BYTE, 
                PROW( snode.Index, this->grid_ ), IDX_TO_TAG(supidx,SELINV_TAG_U_CONTENT), 
                this->grid_->colComm, mpireqsRecvFromAbove );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving U " << snode.SizeSstrUrowRecv << " BYTES"<< std::endl <<  std::endl; 
#endif
          } // if I need to receive from up

          if( this->isRecvFromLeft_( snode.Index ) &&
              MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
            snode.SstrUcolRecv.resize( snode.SizeSstrUcolRecv );
            MPI_Irecv( &snode.SstrUcolRecv[0], snode.SizeSstrUcolRecv, MPI_BYTE, 
                PCOL( snode.Index, this->grid_ ), IDX_TO_TAG(supidx,SELINV_TAG_UCOL_CONTENT), 
                this->grid_->rowComm,
                mpireqsRecvFromLeft );
#if ( _DEBUGlevel_ >= 1 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<  "Receiving Ucol " << snode.SizeSstrUcolRecv << " BYTES"<< std::endl <<  std::endl; 
#endif
          } // if I need to receive from left
        }









      }



      TIMER_START(Compute_Sinv_LU);
      {
        Int gemmProcessed = 0;
        Int gemmToDo = 0;
        //      Int toRecvGemm = 0;
        //copy the list of supernodes we need to process
        std::vector<Int> readySupidx;
        //find local things to do
        for(Int supidx = 0;supidx<stepSuper;supidx++){
          SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
          if( this->isRecvFromAbove_( snode.Index ) && this->isRecvFromLeft_( snode.Index )){
            gemmToDo+=2;
            if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
              snode.isReady+=2;
            }

            if(  MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
              snode.isReady+=2;
            }

            if(snode.isReady==4){
              readySupidx.push_back(supidx);
#if ( _DEBUGlevel_ >= 1 )
              statusOFS<<std::endl<<"Locally processing ["<<snode.Index<<"]"<<std::endl;
#endif
            }
          }
          else{
            if( (this->isRecvFromLeft_( snode.Index )  ) && MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) )
            {
              MPI_Request & mpireqsSendToLeft = arrMpireqsSendLToLeft[supidx];
              //Dummy 0-b send If I was a receiver, I need to send my data to proc in column of snode.Index
              MPI_Isend( NULL, 0, MPI_BYTE, PCOL(snode.Index,this->grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_L_REDUCE), this->grid_->rowComm, &mpireqsSendToLeft );

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYPROC(this->grid_)<<" has sent "<< 0 << " bytes to " << PNUM(MYROW(this->grid_),PCOL(snode.Index,this->grid_),this->grid_) << std::endl;
#endif
            }// if( isRecvFromLeft_( snode.Index ))

            if( (this->isRecvFromAbove_( snode.Index )  ) && MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) )
            {
              MPI_Request & mpireqsSendToAbove = arrMpireqsSendUToAbove[supidx];
              //Dummy 0-b send If I was a receiver, I need to send my data to proc in column of snode.Index
              MPI_Isend( NULL, 0, MPI_BYTE, PROW(snode.Index,this->grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_U_REDUCE), this->grid_->colComm, &mpireqsSendToAbove );

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYPROC(this->grid_)<<" has sent "<< 0 << " bytes to " << PNUM(PROW(snode.Index,this->grid_),MYCOL(this->grid_),this->grid_) << std::endl;
#endif
            }// if( isRecvFromAbove_( snode.Index ) )



          }
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
          do{


            int reqIndicesL[arrMpireqsRecvLContentFromAny.size()];
            int numRecv = 0; 

            //then process with the remote ones

            TIMER_START(WaitContent_LrowL);
#if defined(PROFILE)
            if(begin_SendULWaitContentFirst==0){
              begin_SendULWaitContentFirst=1;
              TIMER_START(WaitContent_LrowL_First);
            }
#endif

            numRecv = 0;
            MPI_Waitsome(2*stepSuper, &arrMpireqsRecvLContentFromAny[0], &numRecv, reqIndicesL, MPI_STATUSES_IGNORE);

            for(int i =0;i<numRecv;i++){
              reqidx = reqIndicesL[i];
              //I've received something
              if(reqidx!=MPI_UNDEFINED)
              { 

                supidx = reqidx/2;
                SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
                snode.isReady++;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS<<std::endl<<"Received data for ["<<snode.Index<<"] reqidx%2="<<reqidx%2<<" is ready ?"<<snode.isReady<<std::endl;
#endif
                //if we received both L and U, the supernode is ready
                if(snode.isReady==4){
                  readySupidx.push_back(supidx);

#if defined(PROFILE)
                  if(end_SendULWaitContentFirst==0){
                    TIMER_STOP(WaitContent_LrowL_First);
                    end_SendULWaitContentFirst=1;
                  }
#endif
                }
              }

            }//end for waitsome
            TIMER_STOP(WaitContent_LrowL);


            int reqIndicesU[arrMpireqsRecvUContentFromAny.size()];
            numRecv = 0; 

            //then process with the remote ones

            TIMER_START(WaitContent_UcolU);
            numRecv = 0;
            MPI_Waitsome(2*stepSuper, &arrMpireqsRecvUContentFromAny[0], &numRecv, reqIndicesU, MPI_STATUSES_IGNORE);

            for(int i =0;i<numRecv;i++){
              reqidx = reqIndicesU[i];
              //I've received something
              if(reqidx!=MPI_UNDEFINED)
              { 

                supidx = reqidx/2;
                SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
                snode.isReady++;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS<<std::endl<<"Received data for ["<<snode.Index<<"] reqidx%2="<<reqidx%2<<" is ready ?"<<snode.isReady<<std::endl;
#endif
                //if we received both L and U, the supernode is ready
                if(snode.isReady==4){
                  readySupidx.push_back(supidx);
                }
              }

            }//end for waitsome
            TIMER_STOP(WaitContent_UcolU);

          } while(readySupidx.size()==0);

          //If I have some work to do 
          if(readySupidx.size()>0)
          {
            supidx = readySupidx.back();
            readySupidx.pop_back();
            SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];


            // Only the processors received information participate in the Gemm 
            if( this->isRecvFromAbove_( snode.Index ) && this->isRecvFromLeft_( snode.Index ) ){

              std::vector<LBlock<T> > LcolRecv;
              std::vector<LBlock<T> > LrowRecv;

              std::vector<UBlock<T> > UcolRecv;
              std::vector<UBlock<T> > UrowRecv;

              // Save all the data to be updated for { L( isup, snode.Index ) | isup > snode.Index }.
              // The size will be updated in the Gemm phase and the reduce phase

              UnpackData(snode, LcolRecv, LrowRecv, UcolRecv, UrowRecv);

              //NumMat<T> AinvBuf, UBuf;
              SelInv_lookup_indexes(snode, LcolRecv, LrowRecv, UcolRecv, UrowRecv,AinvBuf,LBuf, UBuf);

              snode.LUpdateBuf.Resize( AinvBuf.m(), SuperSize( snode.Index, this->super_ ) );
              TIMER_START(Compute_Sinv_LT_GEMM);
              blas::Gemm( 'N', 'T', AinvBuf.m(), LBuf.m(), AinvBuf.n(), MINUS_ONE<T>(), 
                  AinvBuf.Data(), AinvBuf.m(), 
                  LBuf.Data(), LBuf.m(),
                  ZERO<T>(), snode.LUpdateBuf.Data(), snode.LUpdateBuf.m() ); 
              TIMER_STOP(Compute_Sinv_LT_GEMM);

              snode.UUpdateBuf.Resize( SuperSize( snode.Index, this->super_ ), AinvBuf.m() );
              TIMER_START(Compute_Sinv_U_GEMM);
              blas::Gemm( 'N', 'T', UBuf.m(), AinvBuf.m(), AinvBuf.n(), MINUS_ONE<T>(), 
                  UBuf.Data(), UBuf.m(),
                  AinvBuf.Data(), AinvBuf.m(), 
                  ZERO<T>(), snode.UUpdateBuf.Data(), snode.UUpdateBuf.m() ); 
              TIMER_STOP(Compute_Sinv_U_GEMM);

              statusOFS<<"LBuf of "<<snode.Index<<":"<<std::endl;
              statusOFS<<LBuf<<std::endl;

              statusOFS<<"UBuf of "<<snode.Index<<":"<<std::endl;
              statusOFS<<UBuf<<std::endl;

#if ( _DEBUGlevel_ >= 2 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<<  "snode.LUpdateBuf: " << snode.LUpdateBuf << std::endl;
              statusOFS << std::endl << "["<<snode.Index<<"] "<<  "snode.UUpdateBuf: " << snode.UUpdateBuf << std::endl;
#endif
            } // if Gemm is to be done locally


            //If I was a receiver, I need to send my data to proc in column of snode.Index
            if( this->isRecvFromAbove_( snode.Index )  ){
              if( this->isRecvFromLeft_( snode.Index ) && MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) )
              {
                MPI_Request & mpireqsSendToLeft = arrMpireqsSendLToLeft[supidx];

                MPI_Isend( snode.LUpdateBuf.Data(), snode.LUpdateBuf.ByteSize(), MPI_BYTE, PCOL(snode.Index,this->grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_L_REDUCE), this->grid_->rowComm, &mpireqsSendToLeft );

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYCOL(this->grid_)<<" has sent "<< snode.LUpdateBuf.ByteSize() << " bytes to " << PCOL(snode.Index,this->grid_) << std::endl;
#endif

              }//Sender
            }

            //If I was a receiver, I need to send my data to proc in column of snode.Index
            if( this->isRecvFromAbove_( snode.Index )  ){
              if( this->isRecvFromLeft_( snode.Index ) && MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) )
              {
                MPI_Request & mpireqsSendToAbove = arrMpireqsSendUToAbove[supidx];

                MPI_Isend( snode.UUpdateBuf.Data(), snode.UUpdateBuf.ByteSize(), MPI_BYTE, PROW(snode.Index,this->grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_U_REDUCE), this->grid_->colComm, &mpireqsSendToAbove );

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYROW(this->grid_)<<" has sent "<< snode.UUpdateBuf.ByteSize() << " bytes to " << PROW(snode.Index,this->grid_) << std::endl;
#endif

              }//Sender
            }



            gemmProcessed+=2;


#if ( _DEBUGlevel_ >= 1 )
            statusOFS<<std::endl<<"gemmProcessed ="<<gemmProcessed<<"/"<<gemmToDo<<std::endl;
#endif


          }
        }

      }
      TIMER_STOP(Compute_Sinv_LU);

      //Reduce Sinv L to the processors in PCOL(ksup,this->grid_)
      TIMER_START(Reduce_Sinv_L);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];


        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
            //determine the number of rows in LUpdateBufReduced
            Int numRowLUpdateBuf;
            std::vector<LBlock<T> >&  Lcol = this->L( LBj( snode.Index, this->grid_ ) );
            if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
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
                snode.LUpdateBuf.Resize( numRowLUpdateBuf,SuperSize( snode.Index, this->super_ ) );
                // Fill zero is important
                SetValue( snode.LUpdateBuf, ZERO<T>() );
              }
            }

#if ( _DEBUGlevel_ >= 2 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<   "LUpdateBuf Before Reduction: " <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif

            Int totCountRecv = 0;
          Int numRecv = this->CountSendToRight(snode.Index);
            NumMat<T>  LUpdateBufRecv(numRowLUpdateBuf,SuperSize( snode.Index, this->super_ ) );
            for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
              //Do the blocking recv
              MPI_Status stat;
              Int size = 0;
              TIMER_START(L_RECV);
              MPI_Recv(LUpdateBufRecv.Data(), LUpdateBufRecv.ByteSize(), MPI_BYTE, MPI_ANY_SOURCE,IDX_TO_TAG(supidx,SELINV_TAG_L_REDUCE), this->grid_->rowComm,&stat);
              TIMER_STOP(L_RECV);
              MPI_Get_count(&stat, MPI_BYTE, &size);
              //if the processor contributes
              if(size>0){

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYCOL(this->grid_)<<" has received "<< size << " bytes from " << stat.MPI_SOURCE << std::endl;
#endif
#if ( _DEBUGlevel_ >= 2 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<   "LUpdateBufRecv: " <<  LUpdateBufRecv << std::endl << std::endl; 
#endif
                //do the sum
                blas::Axpy(snode.LUpdateBuf.Size(), ONE<T>(), LUpdateBufRecv.Data(), 1, snode.LUpdateBuf.Data(), 1 );
              }
            } // for (iProcCol)
#if ( _DEBUGlevel_ >= 2 ) 
          statusOFS << std::endl << "["<<snode.Index<<"] "<<   "LUpdateBuf After Reduction: " <<  snode.LUpdateBuf << std::endl << std::endl; 
#endif
        } // Receiver
      }

      TIMER_STOP(Reduce_Sinv_L);


      mpi::Waitall( arrMpireqsSendLToLeft );
      //--------------------- End of reduce of LUpdateBuf-------------------------

      //Reduce U Sinv  to the processors in PROW(ksup,this->grid_)
      TIMER_START(Reduce_Sinv_U);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];


        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
            //determine the number of rows in UUpdateBufReduced
            Int numColUUpdateBuf;
            //FIXME U must be revised to store the same structure as L ?
            std::vector<UBlock<T> >&  Urow = this->U( LBi( snode.Index, this->grid_ ) );
            std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi( snode.Index, this->grid_ ) );
            //if( MYCOL( this->grid_ ) != PCOL( snode.Index, this->grid_ ) ){
              snode.ColLocalPtr.resize( Urow.size() + 1 );
              snode.BlockIdxLocalU.resize( Urow.size() );
              snode.ColLocalPtr[0] = 0;

              std::vector<Int> colPtrL(Lrow.size()+1);
                colPtrL[0] = 0;
                 for( Int ib = 0; ib < Lrow.size(); ib++ ){
                colPtrL[ib+1] = colPtrL[ib] + Lrow[ib].numCol;
                  }



              for( Int jb = 0; jb < Urow.size(); jb++ ){
                 Int indexL =0;
                 for( Int ib = 0; ib < Lrow.size(); ib++ ){
                    if(Lrow[ib].blockIdx == Urow[jb].blockIdx){
                      indexL = ib;
                      break;
                    }
                  }
                statusOFS<<jb<<" vs "<<indexL<<std::endl;

                snode.ColLocalPtr[jb] = colPtrL[indexL];
                //snode.ColLocalPtr[jb+1] = snode.ColLocalPtr[jb] + Urow[jb].numCol;
                snode.BlockIdxLocalU[jb] = Lrow[jb].blockIdx;
              }
              snode.ColLocalPtr.back()=colPtrL.back();

            statusOFS<<colPtrL<<std::endl;
            statusOFS<<snode.ColLocalPtr<<std::endl;
            statusOFS<<snode.BlockIdxLocalU<<std::endl;

            //} // I do not own the diagonal block
            numColUUpdateBuf = *snode.ColLocalPtr.rbegin();


            if( numColUUpdateBuf > 0 ){
              if( snode.UUpdateBuf.m() == 0 && snode.UUpdateBuf.n() == 0 ){
                snode.UUpdateBuf.Resize( SuperSize( snode.Index, this->super_ ), numColUUpdateBuf );
                // Fill zero is important
                SetValue( snode.UUpdateBuf, ZERO<T>() );
              }
            }

#if ( _DEBUGlevel_ >= 2 )
            statusOFS << std::endl << "["<<snode.Index<<"] "<<   "UUpdateBuf Before Reduction: " <<  snode.UUpdateBuf << std::endl << std::endl; 
#endif

            Int totCountRecv = 0;

          Int numRecv = this->CountSendToBelow(snode.Index);

            NumMat<T>  UUpdateBufRecv(SuperSize( snode.Index, this->super_ ),numColUUpdateBuf );
            for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
              //Do the blocking recv
              MPI_Status stat;
              Int size = 0;
              TIMER_START(U_RECV);
              MPI_Recv(UUpdateBufRecv.Data(), UUpdateBufRecv.ByteSize(), MPI_BYTE, MPI_ANY_SOURCE,IDX_TO_TAG(supidx,SELINV_TAG_U_REDUCE), this->grid_->colComm,&stat);
              TIMER_STOP(U_RECV);
              MPI_Get_count(&stat, MPI_BYTE, &size);
              //if the processor contributes
              if(size>0){

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYROW(this->grid_)<<" has received "<< size << " bytes from " << stat.MPI_SOURCE << std::endl;
#endif
#if ( _DEBUGlevel_ >= 2 )
                statusOFS << std::endl << "["<<snode.Index<<"] "<<   "UUpdateBufRecv: " <<  UUpdateBufRecv << std::endl << std::endl; 
#endif
                //do the sum
                blas::Axpy(snode.UUpdateBuf.Size(), ONE<T>(), UUpdateBufRecv.Data(), 1, snode.UUpdateBuf.Data(), 1 );
              }
            } // for (iProcRow)
#if ( _DEBUGlevel_ >= 2 ) 
          statusOFS << std::endl << "["<<snode.Index<<"] "<<   "UUpdateBuf After Reduction: " <<  snode.UUpdateBuf << std::endl << std::endl; 
#endif
        } // Receiver
      }

      TIMER_STOP(Reduce_Sinv_U);



      mpi::Waitall( arrMpireqsSendUToAbove );


      //--------------------- End of reduce of UUpdateBuf-------------------------

#ifndef _RELEASE_
      PushCallStack("PMatrixAsym::SelInv_P2p::UpdateD");
#endif

      TIMER_START(Update_Diagonal);
      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

        ComputeDiagUpdate(snode);

        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
          if( MYROW( this->grid_ ) != PROW( snode.Index, this->grid_ ) ){
            if(this->isSendToDiagonal_(snode.Index)){
              //send to above
              MPI_Request & mpireqsSendToAbove = arrMpireqsSendToAbove[supidx];
              MPI_Isend( snode.DiagBuf.Data(),  snode.DiagBuf.ByteSize(), MPI_BYTE,
                  PROW(snode.Index,this->grid_) ,IDX_TO_TAG(supidx,SELINV_TAG_D_REDUCE), this->grid_->colComm, &mpireqsSendToAbove );

#if ( _DEBUGlevel_ >= 1 )
              statusOFS << std::endl << "["<<snode.Index<<"] "<< " P"<<MYROW(this->grid_)<<" has sent "<< snode.DiagBuf.ByteSize() << " bytes of DiagBuf to " << PROW(snode.Index,this->grid_) << " isSendToDiagonal = "<< this->isSendToDiagonal_(snode.Index) <<  std::endl;
#endif
            }
          }
        }
      }

      TIMER_STOP(Update_Diagonal);



      TIMER_START(Reduce_Diagonal);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];
        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) ){
          if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) ){
            if(snode.DiagBuf.Size()==0){
              snode.DiagBuf.Resize( SuperSize( snode.Index, this->super_ ), SuperSize( snode.Index, this->super_ ));
              SetValue(snode.DiagBuf, ZERO<T>());
            }
            //receive from below
            Int totCountRecv = 0;
            Int numRecv = this->CountRecvFromBelow(snode.Index);
            NumMat<T>  DiagBufRecv(snode.DiagBuf.m(),snode.DiagBuf.n());

            for( Int countRecv = 0; countRecv < numRecv ; ++countRecv ){
              //Do the blocking recv
              MPI_Status stat;
              Int size = 0;
              TIMER_START(D_RECV);
              MPI_Recv(DiagBufRecv.Data(), DiagBufRecv.ByteSize(), MPI_BYTE, MPI_ANY_SOURCE,IDX_TO_TAG(supidx,SELINV_TAG_D_REDUCE), this->grid_->colComm,&stat);
              TIMER_STOP(D_RECV);
              MPI_Get_count(&stat, MPI_BYTE, &size);
              //if the processor contributes
              if(size>0){
                // Add DiagBufRecv to diagonal block.
                blas::Axpy(snode.DiagBuf.Size(), ONE<T>(), DiagBufRecv.Data(),
                    1, snode.DiagBuf.Data(), 1 );
              }
            }
            LBlock<T> &  LB = this->L( LBj( snode.Index, this->grid_ ) )[0];
            // Symmetrize LB
            blas::Axpy( LB.numRow * LB.numCol, ONE<T>(), snode.DiagBuf.Data(), 1, LB.nzval.Data(), 1 );
            //Symmetrize( LB.nzval );
          }

        } 
      }


      TIMER_STOP(Reduce_Diagonal);

#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("PMatrixAsym::SelInv_P2p::SendRecvCD");
#endif

      SendRecvCD(arrSuperNodes, stepSuper);

#ifndef _RELEASE_
      PopCallStack();
#endif

#ifndef _RELEASE_
      PushCallStack("PMatrixAsym::SelInv_P2p::UpdateLFinal");
#endif

      TIMER_START(Update_L);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Finish updating the L part by filling LUpdateBufReduced"
          << " back to L" << std::endl << std::endl; 
#endif

        if( MYCOL( this->grid_ ) == PCOL( snode.Index, this->grid_ ) 
            && snode.LUpdateBuf.m() > 0 ){
          std::vector<LBlock<T> >&  Lcol = this->L( LBj( snode.Index, this->grid_ ) );
          //Need to skip the diagonal block if present
          Int startBlock = (MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ))?1:0;
          for( Int ib = startBlock; ib < Lcol.size(); ib++ ){
            LBlock<T> & LB = Lcol[ib];
            lapack::Lacpy( 'A', LB.numRow, LB.numCol, 
                &snode.LUpdateBuf(snode.RowLocalPtr[ib-startBlock], 0),
                snode.LUpdateBuf.m(), LB.nzval.Data(), LB.numRow );
          }
        } // Finish updating L	
      } // for (snode.Index) : Main loop


      TIMER_STOP(Update_L);

#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("PMatrixAsym::SelInv_P2p::UpdateUFinal");
#endif

      TIMER_START(Update_U);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        SuperNodeBufferTypeAsym & snode = arrSuperNodes[supidx];

#if ( _DEBUGlevel_ >= 1 )
        statusOFS << std::endl << "["<<snode.Index<<"] "
          << "Finish updating the U part by filling UUpdateBufReduced"
          << " back to U" << std::endl << std::endl; 
#endif

        if( MYROW( this->grid_ ) == PROW( snode.Index, this->grid_ ) 
            && snode.UUpdateBuf.m() > 0 ){
          std::vector<UBlock<T> >&  Urow = this->U( LBi( snode.Index, this->grid_ ) );
          std::vector<LBlock<T> >&  Lrow = this->Lrow( LBi( snode.Index, this->grid_ ) );
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock<T> & UB = Urow[jb];

            //assert(UB.blockIdx == snode.BlockIdxLocalU[jb]);
            //Inndices follow L order... look at Lrow

            lapack::Lacpy( 'A', UB.numRow, UB.numCol, 
                &snode.UUpdateBuf( 0, snode.ColLocalPtr[jb] ),
                snode.UUpdateBuf.m(), UB.nzval.Data(), UB.numRow );
          }
        } // Finish updating U
      } // for (snode.Index) : Main loop


      TIMER_STOP(Update_U);

#ifndef _RELEASE_
      PopCallStack();
#endif


      TIMER_START(Barrier);
      mpi::Waitall(arrMpireqsRecvLContentFromAny);
      mpi::Waitall(arrMpireqsRecvUContentFromAny);
      //Sync for reduce L
      //      mpi::Waitall( arrMpireqsSendToLeft );
      //Sync for reduce D
      mpi::Waitall(arrMpireqsSendToAbove);

      for (Int supidx=0; supidx<stepSuper; supidx++){
        Int ksup = superList[lidx][supidx];
        std::vector<MPI_Request> & mpireqsSendLToRight = arrMpireqsSendLToRight[supidx];
        std::vector<MPI_Request> & mpireqsSendLToBelow = arrMpireqsSendLToBelow[supidx];
        std::vector<MPI_Request> & mpireqsSendUToRight = arrMpireqsSendUToRight[supidx];
        std::vector<MPI_Request> & mpireqsSendUToBelow = arrMpireqsSendUToBelow[supidx];


        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          MPI_Barrier(this->grid_->colComm);
        }

        mpi::Waitall( mpireqsSendLToRight );
        mpi::Waitall( mpireqsSendLToBelow );
        mpi::Waitall( mpireqsSendUToRight );
        mpi::Waitall( mpireqsSendUToBelow );

      }

      TIMER_STOP(Barrier);


#ifdef LIST_BARRIER
#ifndef ALL_BARRIER
      if (this->options_->maxPipelineDepth!=-1)
#endif
      {
        MPI_Barrier(this->grid_->comm);
      }
#endif

    }

  template<typename T> 
    void PMatrixAsym<T>::SelInv	(  )
    {
      this->SelInv_P2p	(  );
    } 		// -----  end of method PMatrixAsym::SelInv  ----- 

  template<typename T> 
    void PMatrixAsym<T>::SelInv_P2p	(  )
    {
      TIMER_START(SelInv_P2p);

#ifndef _RELEASE_
      PushCallStack("PMatrixAsym::SelInv_P2p");
#endif


      Int numSuper = this->NumSuper(); 

      // Main loop
      std::vector<std::vector<Int> > & superList = this->WorkingSet();
      Int numSteps = superList.size();

      for (Int lidx=0; lidx<numSteps ; lidx++){
        Int stepSuper = superList[lidx].size(); 

        this->SelInvIntra_P2p(lidx);

//        if(lidx==1){ return;};
      }

#ifndef _RELEASE_
      PopCallStack();
#endif

      TIMER_STOP(SelInv_P2p);

      return ;
    } 		// -----  end of method PMatrixAsym::SelInv_P2p  ----- 

  template<typename T> 
    void PMatrixAsym<T>::PreSelInv	(  )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrixAsym::PreSelInv");
#endif

      Int numSuper = this->NumSuper(); 

#ifndef _RELEASE_
      PushCallStack("L(i,k) <- L(i,k) * L(k,k)^{-1}");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "L(i,k) <- L(i,k) * L(k,k)^{-1}"
        << std::endl << std::endl; 
#endif




      for( Int ksup = 0; ksup < numSuper; ksup++ ){



        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          // Broadcast the diagonal L block
          NumMat<T> nzvalLDiag;
          std::vector<LBlock<T> >& Lcol = this->L( LBj( ksup, this->grid_ ) );
          if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
            nzvalLDiag = Lcol[0].nzval;
            if( nzvalLDiag.m() != SuperSize(ksup, this->super_) ||
                nzvalLDiag.n() != SuperSize(ksup, this->super_) ){
#ifdef USE_ABORT
              abort();
#endif
              throw std::runtime_error( 
                  "The size of the diagonal block of L is wrong." );
            }
          } // Owns the diagonal block
          else {
            nzvalLDiag.Resize(SuperSize(ksup, this->super_), SuperSize(ksup, this->super_));
          }
          MPI_Bcast( (void*)nzvalLDiag.Data(), nzvalLDiag.ByteSize(),
              MPI_BYTE, PROW( ksup, this->grid_ ), this->grid_->colComm );

          // Triangular solve
          for( Int ib = 0; ib < Lcol.size(); ib++ ){
            LBlock<T> & LB = Lcol[ib];
            if( LB.blockIdx > ksup ){
#if ( _DEBUGlevel_ >= 2 )
              // Check the correctness of the triangular solve 
              //for the first local column
//              if( LBj( ksup, this->grid_ ) == 0 ){
//                statusOFS << "Diag   L(" << ksup << ", " << ksup << "): "
//                  << nzvalLDiag << std::endl;
//                statusOFS << "Before solve L(" << LB.blockIdx << ", " << ksup 
//                  << "): " << LB.nzval << std::endl;
//              }
              NumMat<T> Ljk = LB.nzval;
#endif
              blas::Trsm( 'R', 'L', 'N', 'U', LB.numRow, LB.numCol,
                  ONE<T>(), nzvalLDiag.Data(), LB.numCol, 
                  LB.nzval.Data(), LB.numRow );
#if ( _DEBUGlevel_ >= 2 )
              NumMat<T> res = LB.nzval;
              //Compute L'(jk)U(kk) which should be Ljk
              blas::Trmm('R','L','N','U',res.m(),res.n(),ONE<T>(),nzvalLDiag.Data(),nzvalLDiag.m(),res.Data(),res.m());
              blas::Axpy(res.Size(), MINUS_ONE<T>(), Ljk.Data(), 1, res.Data(), 1 );
              double norm = lapack::Lange('F',res.m(),res.n(),res.Data(),res.m(),Ljk.Data());
             
                statusOFS << "After solve norm of residual of L(" << LB.blockIdx << ", " << ksup
                  << "): " << norm << std::endl;






              // Check the correctness of the triangular solve
              // for the first local column
//              if( LBj( ksup, this->grid_ ) == 0 ){
//                statusOFS << "After solve  L(" << LB.blockIdx << ", " << ksup 
//                  << "): " << LB.nzval << std::endl;
//              }
#endif
            }
          }
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      } // for (ksup)


#ifndef _RELEASE_
      PopCallStack();
#endif


#ifndef _RELEASE_
      PushCallStack("U(k,j) <- U(k,k)^{-1} * U(k,j)");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "U(k,j) <- U(k,k)^{-1} * U(k,j)" 
        << std::endl << std::endl; 
#endif
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          // Broadcast the diagonal L block
          NumMat<T> nzvalUDiag;
          std::vector<UBlock<T> >& Urow = this->U( LBi( ksup, this->grid_ ) );
          if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
            std::vector<LBlock<T> >& Lcol = this->L( LBj( ksup, this->grid_ ) );
            nzvalUDiag = Lcol[0].nzval;
            if( nzvalUDiag.m() != SuperSize(ksup, this->super_) ||
                nzvalUDiag.n() != SuperSize(ksup, this->super_) ){
#ifdef USE_ABORT
              abort();
#endif
              throw std::runtime_error( 
                  "The size of the diagonal block of U is wrong." );
            }
          } // Owns the diagonal block
          else {
            nzvalUDiag.Resize(SuperSize(ksup, this->super_), SuperSize(ksup, this->super_));
          }
          MPI_Bcast( (void*)nzvalUDiag.Data(), nzvalUDiag.ByteSize(),
              MPI_BYTE, PCOL( ksup, this->grid_ ), this->grid_->rowComm );

          // Triangular solve
          for( Int jb = 0; jb < Urow.size(); jb++ ){
            UBlock<T> & UB = Urow[jb];
            if( UB.blockIdx > ksup ){
#if ( _DEBUGlevel_ >= 2 )
              // Check the correctness of the triangular solve for the first local column
//              if( LBi( ksup, this->grid_ ) == 0 ){
//                statusOFS << "Diag U(" << ksup << ", " << ksup << "): " 
//                  << nzvalUDiag << std::endl;
//                statusOFS << "Before solve U(" << ksup << ", " << UB.blockIdx
//                  << "): " << UB.nzval << std::endl;
//              }
              NumMat<T> Ukj = UB.nzval;
#endif
              blas::Trsm( 'L', 'U', 'N', 'N', UB.numRow, UB.numCol, 
                  ONE<T>(), nzvalUDiag.Data(), UB.numRow,
                  UB.nzval.Data(), UB.numRow );
#if ( _DEBUGlevel_ >= 2 )
              NumMat<T> res = UB.nzval;
              //Compute U(kk) * U^(kj) which should be Ukj
              blas::Trmm('L','U','N','N',res.m(),res.n(),ONE<T>(),nzvalUDiag.Data(),nzvalUDiag.m(),res.Data(),res.m());
              blas::Axpy(res.Size(), MINUS_ONE<T>(), Ukj.Data(), 1, res.Data(), 1 );
              double norm = lapack::Lange('F',res.m(),res.n(),res.Data(),res.m(),Ukj.Data());
              statusOFS << "After solve, norm of residual of U(" << ksup << ", " << UB.blockIdx
                  << "): " << norm << std::endl;

 
              // Check the correctness of the triangular solve for the first local column
//              if( LBi( ksup, this->grid_ ) == 0 ){
//                statusOFS << "After solve  U(" << ksup << ", " << UB.blockIdx
//                  << "): " << UB.nzval << std::endl;
//              }
#endif
            }
          }
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for (ksup)


#ifndef _RELEASE_
      PopCallStack();
#endif




#ifndef _RELEASE_
      PushCallStack("L(i,i) <- [L(k,k) * U(k,k)]^{-1} ");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "L(i,i) <- [L(k,k) * U(k,k)]^{-1}" << std::endl 
        << std::endl; 
#endif

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) &&
            MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ )	){
          IntNumVec ipiv( SuperSize( ksup, this->super_ ) );
          // Note that the pivoting vector ipiv should follow the FORTRAN
          // notation by adding the +1
          for(Int i = 0; i < SuperSize( ksup, this->super_ ); i++){
            ipiv[i] = i + 1;
          }
          LBlock<T> & LB = (this->L( LBj( ksup, this->grid_ ) ))[0];
#if ( _DEBUGlevel_ >= 2 )
          // Check the correctness of the matrix inversion 
          // for the first local column
//          statusOFS << "Factorized A (" << ksup << ", " << ksup << "): "
//            << LB.nzval << std::endl;
          NumMat<T> Lkk = LB.nzval;
#endif
          lapack::Getri( SuperSize( ksup, this->super_ ), LB.nzval.Data(), 
              SuperSize( ksup, this->super_ ), ipiv.Data() );

#if ( _DEBUGlevel_ >= 2 )
              NumMat<T> res (LB.nzval.m(), LB.nzval.n());
              NumMat<T> Akk (LB.nzval.m(), LB.nzval.n());
              //rebuild A
              SetValue(Akk,ZERO<T>());
              //copy u into A
              for(Int i = 0; i<Akk.m();++i){
                for(Int j = i; j<Akk.n();++j){
                  Akk(i,j)= Lkk(i,j);
                }
              }
              blas::Trmm('L','L','N','U',Akk.m(),Akk.n(),ONE<T>(),Lkk.Data(),Lkk.m(),Akk.Data(),Akk.m());
            
//              statusOFS << "After inversion, original A(" << ksup << ", " << ksup
//                  << "): " << Akk << std::endl;

              //Compute U(kk) * U'(kj) which should be Ukj
              blas::Gemm('N','N',res.m(),res.n(),res.m(),ONE<T>(),LB.nzval.Data(),res.m(),Akk.Data(),res.m(),ZERO<T>(),res.Data(),res.m());
              for(Int i = 0; i<res.m();++i){
                res(i,i)-=ONE<T>();
              }

//              statusOFS << "After inversion, residual of A(" << ksup << ", " << ksup
//                  << "): " << res << std::endl;

              double norm = lapack::Lange('F',res.m(),res.n(),res.Data(),res.m(),Lkk.Data());
              statusOFS << "After inversion, norm of residual of A(" << ksup << ", " << ksup
                  << "): " << norm << std::endl;




          // Check the correctness of the matrix inversion 
          // for the first local column
//          statusOFS << "Inverted   A (" << ksup << ", " << ksup << "): " 
//            << LB.nzval << std::endl;
#endif
        } // if I need to invert the diagonal block
      } // for (ksup)


#ifndef _RELEASE_
      PopCallStack();
#endif



/*
      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        //Build an array of SuperNodeBufferTypeAsym
        std::vector<SuperNodeBufferTypeAsym> arrSuperNodes(1);
        //allocate the buffers for this supernode
        for (Int supidx=0; supidx<arrSuperNodes.size(); supidx++){ 
          arrSuperNodes[supidx].Index = ksup;  
        }

        SendRecvCD(arrSuperNodes, arrSuperNodes.size());
      }
*/







#ifndef _RELEASE_
      PushCallStack("Lrow(k,i) <- L(i,k)");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Lrow(k,i) <- L(i,k)" << std::endl << std::endl; 
#endif

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        Int ksupProcRow = PROW( ksup, this->grid_ );
        Int ksupProcCol = PCOL( ksup, this->grid_ );

        Int sendCount = this->CountSendToCrossDiagonal(ksup);
        Int recvCount = this->CountRecvFromCrossDiagonal(ksup);

        std::vector<MPI_Request > arrMpiReqsSend(sendCount, MPI_REQUEST_NULL );
        std::vector<MPI_Request > arrMpiReqsSizeSend(sendCount, MPI_REQUEST_NULL);
        std::vector<std::vector<char> > arrSstrLcolSend(sendCount);
        std::vector<Int > arrSstrLcolSizeSend(sendCount);

        std::vector<MPI_Request > arrMpiReqsRecv(recvCount, MPI_REQUEST_NULL );
        std::vector<MPI_Request > arrMpiReqsSizeRecv(recvCount, MPI_REQUEST_NULL);
        std::vector<std::vector<char> > arrSstrLcolRecv(recvCount);
        std::vector<Int > arrSstrLcolSizeRecv(recvCount);



        // Sender of L
        if( this->isSendToCrossDiagonal_(this->grid_->numProcCol,ksup) ){
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" should send to "
            << this->CountSendToCrossDiagonal(ksup)<<" processors"
            << std::endl;
#endif

          Int sendIdx = 0;
          for(Int dstCol = 0; dstCol<this->grid_->numProcCol; dstCol++){
            if(this->isSendToCrossDiagonal_(dstCol,ksup) ){
              Int dst = PNUM(PROW(ksup,this->grid_),dstCol,this->grid_);
              if(MYPROC(this->grid_)!= dst){
                // Pack L data
                std::stringstream sstm;
                std::vector<char> & sstrLcolSend = arrSstrLcolSend[sendIdx];
                Int & sstrSize = arrSstrLcolSizeSend[sendIdx];
                MPI_Request & mpiReqSend = arrMpiReqsSend[sendIdx];
                MPI_Request & mpiReqSizeSend = arrMpiReqsSizeSend[sendIdx];

                std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
                std::vector<LBlock<T> >&  Lcol = this->L( LBj(ksup, this->grid_) );
                // All blocks except for the diagonal block are to be sent right
                //TODO not true > this is a scatter operation ! Can we know the destination ?

                //Skip the diagonal block if necessary 
                Int startIdx = ( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ) ? 1:0;
//                Int startIdx = ( MYPROC(this->grid_) == PNUM(ksupProcRow,ksupProcCol, this->grid_) ) ? 1:0;
                Int count = 0;
                for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > ksup && 
                      (Lcol[ib].blockIdx % this->grid_->numProcCol) == dstCol  ){
                    count++;
                  }
                }

                serialize( (Int)count, sstm, NO_MASK );

                for( Int ib = startIdx; ib < Lcol.size(); ib++ ){
                  if( Lcol[ib].blockIdx > ksup &&  
                      (Lcol[ib].blockIdx % this->grid_->numProcCol) == dstCol  ){ 
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS << "["<<ksup<<"] SEND contains "<<Lcol[ib].blockIdx
                      << " which corresponds to "<<GBj(ib,this->grid_)
                      << std::endl;
#endif
                    serialize( Lcol[ib], sstm, mask );
                  }
                }

                sstrLcolSend.resize( Size(sstm) );
                sstm.read( &sstrLcolSend[0], sstrLcolSend.size() );
                sstrSize = sstrLcolSend.size();

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") ---> LBj("<<ksup<<")="
                  << LBj(ksup,this->grid_)<<" ---> P"<<dst<<" ("<<ksupProcRow
                  << ","<<dstCol<<")"<<std::endl;
#endif
                MPI_Isend( &sstrSize, 1, MPI_INT, dst, SELINV_TAG_D_SIZE,
                    this->grid_->comm, &mpiReqSizeSend );
                MPI_Isend( (void*)&sstrLcolSend[0], sstrSize, MPI_BYTE, dst,
                    SELINV_TAG_D_CONTENT, this->grid_->comm, &mpiReqSend );
                sendIdx++;
              }
            }
          }
        } // if I am a sender of L

        // Receiver of L
        if( this->isRecvFromCrossDiagonal_(this->grid_->numProcRow,ksup) ){

#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" should receive from "
            << this->CountRecvFromCrossDiagonal(ksup)<<" processors"
            << std::endl;
#endif

          std::vector<UBlock<T> >& Urow = this->U( LBi( ksup, this->grid_ ) );
          std::vector<LBlock<T> >& Lrow = this->Lrow( LBi( ksup, this->grid_ ) );
          std::vector<LBlock<T> >& Lcol = this->L( LBj( ksup, this->grid_ ) );


          statusOFS<<"Lrow : "<<std::endl;
          for(Int ib=0;ib<Lrow.size();++ib){statusOFS<<Lrow[ib].blockIdx<<" ";}
          statusOFS<<std::endl;

          statusOFS<<"Lcol : "<<std::endl;
          for(Int ib=0;ib<Lcol.size();++ib){statusOFS<<Lcol[ib].blockIdx<<" ";}
          statusOFS<<std::endl;

          std::vector<Int> isBlockFound(Lrow.size(),ZERO<Int>());
          //If I own the diagonal block, mark it as found
//          if(MYPROC(this->grid_)== PNUM(PROW(ksup,this->grid_),PCOL(ksup,this->grid_),this->grid_)){
//            isBlockFound[0] = true;
//          }

          Int recvIdx = 0;
          //receive size first
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,ksup) ){
              std::vector<LBlock<T> > LcolRecv;
              Int src = PNUM(srcRow,PCOL(ksup,this->grid_),this->grid_);
              if(MYPROC(this->grid_)!= src){
                MPI_Request & mpiReqSizeRecv = arrMpiReqsSizeRecv[recvIdx];
                Int & sstrSize = arrSstrLcolSizeRecv[recvIdx];

                MPI_Irecv( &sstrSize, 1, MPI_INT, src, SELINV_TAG_D_SIZE, 
                    this->grid_->comm, &mpiReqSizeRecv );
                recvIdx++;
              }
            }
          }

          mpi::Waitall(arrMpiReqsSizeRecv);

          //receive content
          recvIdx = 0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,ksup) ){
              Int src = PNUM(srcRow,PCOL(ksup,this->grid_),this->grid_);
              if(MYPROC(this->grid_)!= src){
                MPI_Request & mpiReqRecv = arrMpiReqsRecv[recvIdx];
                Int & sstrSize = arrSstrLcolSizeRecv[recvIdx];
                std::vector<char> & sstrLcolRecv = arrSstrLcolRecv[recvIdx];
                sstrLcolRecv.resize(sstrSize);

                MPI_Irecv( &sstrLcolRecv[0], sstrSize, MPI_BYTE, src, 
                    SELINV_TAG_D_CONTENT, this->grid_->comm, &mpiReqRecv );
                recvIdx++;
              }
            }
          }

          mpi::Waitall(arrMpiReqsRecv);


          //Process the content
          recvIdx = 0;
          for(Int srcRow = 0; srcRow<this->grid_->numProcRow; srcRow++){
            if(this->isRecvFromCrossDiagonal_(srcRow,ksup) ){
              std::vector<LBlock<T> > * pLcol;
              std::vector<LBlock<T> > LcolRecv;
              Int src = PNUM(srcRow,PCOL(ksup,this->grid_),this->grid_);
              if(MYPROC(this->grid_)!= src){

                Int & sstrSize = arrSstrLcolSizeRecv[recvIdx];
                std::vector<char> & sstrLcolRecv = arrSstrLcolRecv[recvIdx];
                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") <--- LBj("<<ksup<<") <--- P"
                  << src <<" ("<<srcRow<<","<<ksupProcCol<<")"
                  << std::endl;
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
                  statusOFS << "["<<ksup<<"] RECV contains "
                    << LcolRecv[ib].blockIdx<< " which corresponds to "
                    << ib * this->grid_->numProcRow + srcRow;
                  statusOFS << " L is on row "<< srcRow 
                    << " whereas Lrow is on col "
                    << LcolRecv[ib].blockIdx % this->grid_->numProcCol 
                    << std::endl;
#endif
                }
                recvIdx++;
                pLcol = &LcolRecv;
              } // sender is not the same as receiver
              else{
                // L is obtained locally, just make a copy. Do not include the diagonal block
                //              Int offset = ( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )?1:0;
                //              LcolRecv.resize( Lcol.size() - offset );
                //              for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
                //                LcolRecv[ib] = Lcol[ib+offset];
                //              }

                pLcol = &Lcol;
              } // sender is the same as receiver

              //Update Lrow
              //We can directly put LcolRecv in Lrow (stil need to transpose)

              //Skip the diagonal block if necessary 
              Int startIdx = ( MYPROC( this->grid_ ) == src && MYPROC(this->grid_) == PNUM(ksupProcRow,ksupProcCol, this->grid_) ) ? 1:0;
              for( Int ib = startIdx; ib < pLcol->size(); ib++ ){
                LBlock<T> & LB = (*pLcol)[ib];
                if( LB.blockIdx <= ksup ){
#ifdef USE_ABORT
                  abort();
#endif
                  throw std::logic_error( "LcolRecv contains the wrong blocks." );
                }

                //std::vector<LBlock<T> > &  LrowCD = this->Lrow( LBi( ksup, super_ ) );
                for( Int iib = 0; iib < Lrow.size(); iib++ ){
                  LBlock<T> &  LrowB = Lrow[ iib ];
                  if( LB.blockIdx == LrowB.blockIdx ){
                    // Compare size (LrowB is the transpose)
//                    if( LB.numRow != LrowB.numCol || LB.numCol != LrowB.numRow ){
//                      std::ostringstream msg;
//                      msg << "LB(" << LB.blockIdx << ", " << ksup 
//                        << ") and LrowB(" << LrowB.blockIdx << ", " << ksup 
//                        << ")	do not share the same size." << std::endl
//                        << "LB: " << LB.numRow << " x " << LB.numCol 
//                        << std::endl
//                        << "LrowB: " << LrowB.numRow << " x " << LrowB.numCol 
//                        << std::endl;
//#ifdef USE_ABORT
//                      abort();
//#endif
//                      throw std::runtime_error( msg.str().c_str() );
//                    }

                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    LrowB = LB;
                    //Store in "row major" format / i.e transpose is in col-major
                    Transpose(LB.nzval, LrowB.nzval);
                    LrowB.numCol = LB.numRow;
                    LrowB.numRow = LB.numCol;

#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<ksup<<"] USING LB "<<LB.blockIdx<< std::endl;
#endif
                    isBlockFound[iib] = 1;//true;
                    break;
                  } // if( LB.blockIdx == LrowB.blockIdx )
                } // for (iib)
              } // for (ib)
            }
          }

          for( Int ib = 0; ib < Lrow.size(); ib++ ){
            if( !isBlockFound[ib] ){
#ifdef USE_ABORT
              LBlock<T> & LB = Lrow[ib];
              statusOFS<<isBlockFound<<std::endl;
              abort();
#endif
              throw std::logic_error( 
                  "LBlock cannot find its update. Something is seriously wrong." );
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
      PushCallStack("Ucol(j,k) <- U(k,j)");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "Ucol(j,k) <- U(k,j)" << std::endl << std::endl;
#endif

      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        Int ksupProcRow = PROW( ksup, this->grid_ );
        Int ksupProcCol = PCOL( ksup, this->grid_ );

        Int recvCount = this->CountSendToCrossDiagonal(ksup);
        Int sendCount = this->CountRecvFromCrossDiagonal(ksup);

        std::vector<MPI_Request > arrMpiReqsSend(sendCount, MPI_REQUEST_NULL );
        std::vector<MPI_Request > arrMpiReqsSizeSend(sendCount, MPI_REQUEST_NULL);
        std::vector<std::vector<char> > arrSstrUrowSend(sendCount);
        std::vector<Int > arrSstrUrowSizeSend(sendCount);

        std::vector<MPI_Request > arrMpiReqsRecv(recvCount, MPI_REQUEST_NULL );
        std::vector<MPI_Request > arrMpiReqsSizeRecv(recvCount, MPI_REQUEST_NULL);
        std::vector<std::vector<char> > arrSstrUrowRecv(recvCount);
        std::vector<Int > arrSstrUrowSizeRecv(recvCount);

        // Sender for U is a receiver for L
        if( this->isRecvFromCrossDiagonal_(this->grid_->numProcRow,ksup) ){
#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" should send to "
            << this->CountRecvFromCrossDiagonal(ksup)<<" processors"
            << std::endl;
#endif

          Int sendIdx = 0;
          for(Int dstRow = 0; dstRow<this->grid_->numProcRow; dstRow++){
            if(this->isRecvFromCrossDiagonal_(dstRow,ksup) ){
              Int dst = PNUM(dstRow,PCOL(ksup,this->grid_),this->grid_);
              if(MYPROC(this->grid_)!= dst){
                // Pack U data
                std::stringstream sstm;
                std::vector<char> & sstrUrowSend = arrSstrUrowSend[sendIdx];
                Int & sstrSize = arrSstrUrowSizeSend[sendIdx];
                MPI_Request & mpiReqSend = arrMpiReqsSend[sendIdx];
                MPI_Request & mpiReqSizeSend = arrMpiReqsSizeSend[sendIdx];

                std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
                std::vector<UBlock<T> >&  Urow = this->U( LBi(ksup, this->grid_) );
                // All blocks except for the diagonal block are to be sent right
                //TODO not true > this is a scatter operation ! Can we know the destination ?

                Int count = 0;
                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  if( Urow[jb].blockIdx > ksup 
                      && (Urow[jb].blockIdx % this->grid_->numProcRow) == dstRow ){
                    count++;
                  }
                }

                serialize( (Int)count, sstm, NO_MASK );

                for( Int jb = 0; jb < Urow.size(); jb++ ){
                  if( Urow[jb].blockIdx > ksup 
                      && (Urow[jb].blockIdx % this->grid_->numProcRow) == dstRow ){ 
#if ( _DEBUGlevel_ >= 1 )
                    statusOFS << "["<<ksup<<"] SEND contains "<<Urow[jb].blockIdx
                      << " which corresponds to "<<GBi(jb,this->grid_)
                      << std::endl;
#endif
                    serialize( Urow[jb], sstm, mask );
                  }
                }

                sstrUrowSend.resize( Size(sstm) );
                sstm.read( &sstrUrowSend[0], sstrUrowSend.size() );
                sstrSize = sstrUrowSend.size();

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") ---> LBi("<<ksup<<")="
                  << LBi(ksup,this->grid_)<<" ---> P"<<dst<<" ("<<ksupProcCol
                  << ","<<dstRow<<")"<<std::endl;
#endif
                MPI_Isend( &sstrSize, 1, MPI_INT, dst, SELINV_TAG_D_SIZE, 
                    this->grid_->comm, &mpiReqSizeSend );
                MPI_Isend( (void*)&sstrUrowSend[0], sstrSize, MPI_BYTE, dst, 
                    SELINV_TAG_D_CONTENT, this->grid_->comm, &mpiReqSend );
                sendIdx++;
              }
            }
          }
        } // if I am a sender of U

        // Receiver for U is a sender for L
        if( this->isSendToCrossDiagonal_(this->grid_->numProcCol,ksup) ){


#if ( _DEBUGlevel_ >= 1 )
          statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" should receive from "
            << this->CountSendToCrossDiagonal(ksup)<<" processors"
            << std::endl;
#endif

          std::vector<UBlock<T> >& Ucol = this->Ucol( LBj( ksup, this->grid_ ) );
          std::vector<Int> isBlockFound(Ucol.size(),0);

          Int recvIdx = 0;
          //receive size first
          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,ksup) ){
              std::vector<UBlock<T> > UrowRecv;
              Int src = PNUM(PROW(ksup,this->grid_),srcCol,this->grid_);
              if(MYPROC(this->grid_)!= src){
                MPI_Request & mpiReqSizeRecv = arrMpiReqsSizeRecv[recvIdx];
                Int & sstrSize = arrSstrUrowSizeRecv[recvIdx];
                MPI_Irecv( &sstrSize, 1, MPI_INT, src, SELINV_TAG_D_SIZE,
                    this->grid_->comm, &mpiReqSizeRecv );
                recvIdx++;
              }
            }
          }

          mpi::Waitall(arrMpiReqsSizeRecv);

          //receive content
          recvIdx = 0;
          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,ksup) ){
              std::vector<UBlock<T> > UrowRecv;
              Int src = PNUM(PROW(ksup,this->grid_),srcCol,this->grid_);
              if(MYPROC(this->grid_)!= src){
                MPI_Request & mpiReqRecv = arrMpiReqsRecv[recvIdx];
                Int & sstrSize = arrSstrUrowSizeRecv[recvIdx];
                std::vector<char> & sstrUrowRecv = arrSstrUrowRecv[recvIdx];
                sstrUrowRecv.resize(sstrSize);
                MPI_Irecv( &sstrUrowRecv[0], sstrSize, MPI_BYTE, src, 
                    SELINV_TAG_D_CONTENT, this->grid_->comm, &mpiReqRecv );
                recvIdx++;
              }
            }
          }

          mpi::Waitall(arrMpiReqsRecv);

          //Process the content
          recvIdx = 0;
          for(Int srcCol = 0; srcCol<this->grid_->numProcCol; srcCol++){
            if(this->isSendToCrossDiagonal_(srcCol,ksup) ){
              std::vector<UBlock<T> > * pUrow;
              std::vector<UBlock<T> > UrowRecv;
              Int src = PNUM(PROW(ksup,this->grid_),srcCol,this->grid_);
              if(MYPROC(this->grid_)!= src){
                Int & sstrSize = arrSstrUrowSizeRecv[recvIdx];
                std::vector<char> & sstrUrowRecv = arrSstrUrowRecv[recvIdx];
                std::stringstream sstm;

#if ( _DEBUGlevel_ >= 1 )
                statusOFS << "["<<ksup<<"] P"<<MYPROC(this->grid_)<<" ("<<MYROW(this->grid_)
                  << ","<<MYCOL(this->grid_)<<") <--- LBi("<<ksup<<") <--- P"
                  << src<<" ("<<ksupProcCol<<","<<srcCol<<")"<<std::endl;
#endif
                sstm.write( &sstrUrowRecv[0], sstrSize );

                // Unpack U data.  
                Int numUBlock;
                std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
                deserialize( numUBlock, sstm, NO_MASK );
                UrowRecv.resize(numUBlock);
                for( Int jb = 0; jb < numUBlock; jb++ ){
                  deserialize( UrowRecv[jb], sstm, mask );
#if ( _DEBUGlevel_ >= 1 )
                  statusOFS << "["<<ksup<<"] RECV contains "
                    << UrowRecv[jb].blockIdx<< " which corresponds to "
                    << jb * this->grid_->numProcCol + srcCol;
                  statusOFS << " Urow is on col "<< srcCol 
                    << " whereas U is on row "
                    << UrowRecv[jb].blockIdx % this->grid_->numProcRow 
                    << std::endl;
#endif
                }
                pUrow = &UrowRecv;
                recvIdx++;
              } // sender is not the same as receiver
              else{
                // U is obtained locally, just make a copy. 
                std::vector<UBlock<T> >& Urow = this->U( LBi( ksup, this->grid_ ) );
                pUrow = &Urow;
                //              Int offset = ( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )?1:0;
                //              UrowRecv.resize( Urow.size() - offset );
                //              for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
                //                UrowRecv[jb] = Urow[jb+offset];
                //              }
              } // sender is the same as receiver

              //Update Ucol
              //We can directly put UrowRecv in Ucol
              for( Int jb = 0; jb < pUrow->size(); jb++ ){
                UBlock<T> & UB = (*pUrow)[jb];
                if( UB.blockIdx <= ksup ){
#ifdef USE_ABORT
                  abort();
#endif
                  throw std::logic_error( "UrowRecv contains the wrong blocks." );
                }

                //std::vector<UBlock<T> > &  UcolCD = this->Ucol( LBi( ksup, super_ ) );
                for( Int jjb = 0; jjb < Ucol.size(); jjb++ ){
                  UBlock<T> &  UcolB = Ucol[jjb];
                  if( UB.blockIdx == UcolB.blockIdx ){
                    // Compare size
//                    if( UB.numRow != UcolB.numRow || UB.numCol != UcolB.numCol ){
//                      std::ostringstream msg;
//                      msg << "UB(" << ksup << ", " << UB.blockIdx << ") and UcolB(" 
//                        << ksup << ", " << UcolB.blockIdx 
//                        << ")	do not share the same size." << std::endl
//                        << "UB: " << UB.numRow << " x " << UB.numCol 
//                        << std::endl
//                        << "UcolB: " << UcolB.numRow << " x " << UcolB.numCol 
//                        << std::endl;
//#ifdef USE_ABORT
//                      abort();
//#endif
//                      throw std::runtime_error( msg.str().c_str() );
//                    }

                    // Note that the order of the column indices of the U
                    // block may not follow the order of the row indices,
                    // overwrite the information in U.
                    UcolB = UB;

#if ( _DEBUGlevel_ >= 1 )
                    statusOFS<<"["<<ksup<<"] USING UB "<<UB.blockIdx<< std::endl;
#endif
                    isBlockFound[jjb] = 1;
                    break;
                  } // if( UB.blockIdx == UcolB.blockIdx )
                } // for (jjb)
              } // for (jb)
            }
          }

          for( Int jb = 0; jb < Ucol.size(); jb++ ){
            if( !isBlockFound[jb] ){
#ifdef USE_ABORT
              abort();
#endif
              throw std::logic_error( 
                  "UBlock cannot find its update. Something is seriously wrong." );
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
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method PMatrixAsym::PreSelInv  ----- 


  template<typename T>
    void PMatrixAsym<T>::ConstructCommunicationPattern	(  )
    {
      ConstructCommunicationPattern_P2p();
    } 		// -----  end of method PMatrix::ConstructCommunicationPattern  ----- 



  template<typename T>
    void PMatrixAsym<T>::ConstructCommunicationPattern_P2p	(  )
    {
#ifndef _RELEASE_
      PushCallStack("PMatrix::ConstructCommunicationPattern_P2p");
#endif


      Int numSuper = this->NumSuper();

      TIMER_START(Allocate);

#ifndef _RELEASE_
      PushCallStack( "Initialize the communication pattern" );
#endif
      this->isSendToBelow_.Resize(this->grid_->numProcRow, numSuper);
      this->isSendToRight_.Resize(this->grid_->numProcCol, numSuper);
      this->isSendToDiagonal_.Resize( numSuper );
      SetValue( this->isSendToBelow_, false );
      SetValue( this->isSendToRight_, false );
      SetValue( this->isSendToDiagonal_, false );

      this->isSendToCrossDiagonal_.Resize(this->grid_->numProcCol+1, numSuper );
      SetValue( this->isSendToCrossDiagonal_, false );
      this->isRecvFromCrossDiagonal_.Resize(this->grid_->numProcRow+1, numSuper );
      SetValue( this->isRecvFromCrossDiagonal_, false );

      this->isRecvFromAbove_.Resize( numSuper );
      this->isRecvFromLeft_.Resize( numSuper );
      this->isRecvFromBelow_.Resize( this->grid_->numProcRow, numSuper );
      SetValue( this->isRecvFromAbove_, false );
      SetValue( this->isRecvFromBelow_, false );
      SetValue( this->isRecvFromLeft_, false );
#ifndef _RELEASE_
      PopCallStack();
#endif

      TIMER_STOP(Allocate);

      // Skip the communication pattern for one single processor.
      // This does not work for now.
      //    if( this->grid_ -> mpisize == 1 ){
      //#if ( _DEBUGlevel_ >= 0 )
      //      statusOFS << "Skip the construction of communication pattern for "
      //        << "a single processor" << std::endl;
      //      statusOFS << "The values of all communication variables are set to be true, " 
      //        << "but no MPI process will be involved in the selected inversion phase."
      //        << std::endl;
      //#endif
      //      SetValue( this->isSendToBelow_, true );
      //      SetValue( this->isSendToRight_, true );
      //      SetValue( this->isSendToDiagonal_, true );
      //      SetValue( this->isSendToCrossDiagonal_, true );
      //      SetValue( this->isRecvFromCrossDiagonal_, true );
      //      SetValue( this->isRecvFromAbove_, true );
      //      SetValue( this->isRecvFromBelow_, true );
      //      SetValue( this->isRecvFromLeft_, true );
      //      return;
      //    }



      TIMER_START(GetEtree);
      std::vector<Int> snodeEtree(this->NumSuper());
      this->GetEtree(snodeEtree);
      TIMER_STOP(GetEtree);



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


      TIMER_START(Column_communication);


      for( Int ksup = 0; ksup < numSuper; ksup++ ){
        // All block columns perform independently
        std::vector<Int> tAllBlockRowIdx;
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          // Communication
          std::vector<Int> & colBlockIdx = this->ColBlockIdx(LBj(ksup, this->grid_));
          TIMER_START(Allgatherv_Column_communication);
          if( this->grid_ -> mpisize != 1 )
            mpi::Allgatherv( colBlockIdx, tAllBlockRowIdx, this->grid_->colComm );
          else
            tAllBlockRowIdx = colBlockIdx;

          TIMER_STOP(Allgatherv_Column_communication);

          localColBlockRowIdx[LBj( ksup, this->grid_ )].insert(
              tAllBlockRowIdx.begin(), tAllBlockRowIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
          statusOFS 
            << " Column block " << ksup 
            << " has the following nonzero block rows" << std::endl;
          for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, this->grid_ )].begin();
              si != localColBlockRowIdx[LBj( ksup, this->grid_ )].end();
              si++ ){
            statusOFS << *si << "  ";
          }
          statusOFS << std::endl; 
#endif

            std::vector< LBlock<T> > & Lcol = this->L( LBj(ksup, this->grid_ ) );
            statusOFS<<"Lcol of "<<ksup<<std::endl;
            for(Int ib=0;ib<Lcol.size();++ib){statusOFS<<Lcol[ib].blockIdx<<" ";}
            statusOFS<<std::endl;

//          if( MYROW( this->grid_ ) != PROW( ksup, this->grid_ ) ){
            std::vector< UBlock<T> > & Ucol = this->Ucol( LBj(ksup, this->grid_ ) );
            //Allocate Ucol
            for(Int ib = 0; ib < tAllBlockRowIdx.size(); ++ib){
              if(tAllBlockRowIdx[ib] > ksup && 
                   (tAllBlockRowIdx[ib] % this->grid_->numProcRow) == MYROW(this->grid_)  ){ 
                Ucol.push_back(UBlock<T>());
                Ucol.back().blockIdx = tAllBlockRowIdx[ib];
              }
            }
            statusOFS<<"Ucol of "<<ksup<<std::endl;
            for(Int ib=0;ib<Ucol.size();++ib){statusOFS<<Ucol[ib].blockIdx<<" ";}
            statusOFS<<std::endl;



            if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
              assert(Lcol.size() == Ucol.size()+1);
            }
            else{
              assert(Lcol.size() == Ucol.size());
            }


//          }

        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )

        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          // Communication
          TIMER_START(Allgatherv_Column_communication_row);
          if( this->grid_ -> mpisize != 1 )
          {
            //Broadcast size
            Int sizeLcol = tAllBlockRowIdx.size();
            MPI_Bcast(&sizeLcol,sizeof(sizeLcol),MPI_BYTE,PCOL(ksup, this->grid_),this->grid_->rowComm);
            if(sizeLcol>0){
              if(MYCOL(this->grid_) != PCOL(ksup, this->grid_)){
                tAllBlockRowIdx.resize(sizeLcol);
              }
              MPI_Bcast(&tAllBlockRowIdx[0],sizeLcol*sizeof(Int),MPI_BYTE,PCOL(ksup, this->grid_),this->grid_->rowComm);
            }
          }

          TIMER_STOP(Allgatherv_Column_communication_row);



//          gdb_lock();

          if( MYCOL( this->grid_ ) != PCOL( ksup, this->grid_ ) ){
            //          std::vector< LBlock<T> > & Lcol = this->L( LBj(ksup, this->grid_ ) );
            std::vector< UBlock<T> > & Urow = this->U( LBi(ksup, this->grid_ ) );
            //Allocate Lrow and extend Urow
            std::vector< LBlock<T> > & Lrow = this->Lrow( LBi(ksup, this->grid_ ) );
            for(Int ib = 0; ib < tAllBlockRowIdx.size(); ++ib){
              if(tAllBlockRowIdx[ib] > ksup && 
                   (tAllBlockRowIdx[ib] % this->grid_->numProcCol) == MYCOL(this->grid_)  ){ 

                Lrow.push_back(LBlock<T>());
                Lrow.back().blockIdx = tAllBlockRowIdx[ib];


                bool isFound = false;
                for(Int jb = 0; jb < Urow.size(); ++jb){
                  if(Urow[jb].blockIdx == tAllBlockRowIdx[ib]){
                    isFound = true;
                    break;
                  }
                }
                if(!isFound){
                  Urow.push_back(UBlock<T>());
                  Urow.back().blockIdx = tAllBlockRowIdx[ib];
                }

              }

            }

            statusOFS<<"Lrow of "<<ksup<<std::endl;
            for(Int ib=0;ib<Lrow.size();++ib){statusOFS<<Lrow[ib].blockIdx<<" ";}
            statusOFS<<std::endl;
            assert(Urow.size() == Lrow.size());
          }
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for(ksup)


#ifndef _RELEASE_
      PopCallStack();
#endif

      TIMER_STOP(Column_communication);

      TIMER_START(Row_communication);
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
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){

          // Communication
          std::vector<Int> tAllBlockColIdx;
          std::vector<Int> & rowBlockIdx = this->RowBlockIdx(LBi(ksup, this->grid_));
          TIMER_START(Allgatherv_Row_communication);
          if( this->grid_ -> mpisize != 1 )
            mpi::Allgatherv( rowBlockIdx, tAllBlockColIdx, this->grid_->rowComm );
          else
            tAllBlockColIdx = rowBlockIdx;

          TIMER_STOP(Allgatherv_Row_communication);

          localRowBlockColIdx[LBi( ksup, this->grid_ )].insert(
              tAllBlockColIdx.begin(), tAllBlockColIdx.end() );

#if ( _DEBUGlevel_ >= 1 )
          statusOFS 
            << " Row block " << ksup 
            << " has the following nonzero block columns" << std::endl;
          for( std::set<Int>::iterator si = localRowBlockColIdx[LBi( ksup, this->grid_ )].begin();
              si != localRowBlockColIdx[LBi( ksup, this->grid_ )].end();
              si++ ){
            statusOFS << *si << "  ";
          }
          statusOFS << std::endl; 
#endif

            std::vector< UBlock<T> > & Urow = this->U( LBi(ksup, this->grid_ ) );
            statusOFS<<"Urow of "<<ksup<<std::endl;
            for(Int ib=0;ib<Urow.size();++ib){statusOFS<<Urow[ib].blockIdx<<" ";}
            statusOFS<<std::endl;
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for(ksup)

#ifndef _RELEASE_
      PopCallStack();
#endif

      TIMER_STOP(Row_communication);

      TIMER_START(STB_RFA);
#ifndef _RELEASE_
      PushCallStack("SendToBelow / RecvFromAbove");
#endif
      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        // Loop over all the supernodes to the right of ksup


        Int jsup = snodeEtree[ksup];
        while(jsup<numSuper){
          Int jsupLocalBlockCol = LBj( jsup, this->grid_ );
          Int jsupProcCol = PCOL( jsup, this->grid_ );
          if( MYCOL( this->grid_ ) == jsupProcCol ){

            // SendToBelow / RecvFromAbove only if (ksup, jsup) is nonzero.
            if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ) {
              for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
                  si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
                Int isup = *si;
                Int isupProcRow = PROW( isup, this->grid_ );
                if( isup > ksup ){
                  if( MYROW( this->grid_ ) == isupProcRow ){
                    this->isRecvFromAbove_(ksup) = true;
                  }
                  if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
                    this->isSendToBelow_( isupProcRow, ksup ) = true;
                  }
                } // if( isup > ksup )
              } // for (si)
            } // if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 )

          } // if( MYCOL( this->grid_ ) == PCOL( jsup, this->grid_ ) )
          jsup = snodeEtree[jsup];

        } // for(jsup)
      } // for(ksup)

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "this->isSendToBelow:" << std::endl;
      for(int j = 0;j< this->isSendToBelow_.n();j++){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < this->isSendToBelow_.m();i++){
          statusOFS<< this->isSendToBelow_(i,j) << " ";
        }
        statusOFS<<std::endl;
      }

      statusOFS << std::endl << "this->isRecvFromAbove:" << std::endl;
      for(int j = 0;j< this->isRecvFromAbove_.m();j++){
        statusOFS << "["<<j<<"] "<< this->isRecvFromAbove_(j)<<std::endl;
      }
#endif

#ifndef _RELEASE_
      PopCallStack();
#endif


      TIMER_STOP(STB_RFA);








      TIMER_START(STR_RFL_RFB);


#ifndef _RELEASE_
      PushCallStack("SendToRight / RecvFromLeft");
#endif
      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        // Loop over all the supernodes below ksup

        Int isup = snodeEtree[ksup];
        while(isup<numSuper){
          Int isupLocalBlockRow = LBi( isup, this->grid_ );
          Int isupProcRow       = PROW( isup, this->grid_ );
          if( MYROW( this->grid_ ) == isupProcRow ){
            // SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.
            if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
              for( std::set<Int>::iterator si = localRowBlockColIdx[isupLocalBlockRow].begin();
                  si != localRowBlockColIdx[isupLocalBlockRow].end(); si++ ){
                Int jsup = *si;
                Int jsupProcCol = PCOL( jsup, this->grid_ );
                if( jsup > ksup ){

                  if( MYCOL( this->grid_ ) == jsupProcCol ){
                    this->isRecvFromLeft_(ksup) = true;
                  }
                  if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
                    this->isSendToRight_( jsupProcCol, ksup ) = true;
                  }
                }
              } // for (si)
            } // if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 )
          } // if( MYROW( this->grid_ ) == isupProcRow )


          if( MYCOL( this->grid_ ) == PCOL(ksup, this->grid_) ){

            if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){ 
              this->isRecvFromBelow_(isupProcRow,ksup) = true;
            }    
            else if (MYROW(this->grid_) == isupProcRow){
              this->isSendToDiagonal_(ksup)=true;
            }    
          } // if( MYCOL( this->grid_ ) == PCOL(ksup, this->grid_) )
          isup = snodeEtree[isup];

        } // for (isup)
      }	 // for (ksup)


#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "this->isSendToRight:" << std::endl;
      for(int j = 0;j< this->isSendToRight_.n();j++){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < this->isSendToRight_.m();i++){
          statusOFS<< this->isSendToRight_(i,j) << " ";
        }
        statusOFS<<std::endl;
      }

      statusOFS << std::endl << "this->isRecvFromLeft:" << std::endl;
      for(int j = 0;j< this->isRecvFromLeft_.m();j++){
        statusOFS << "["<<j<<"] "<< this->isRecvFromLeft_(j)<<std::endl;
      }

      statusOFS << std::endl << "this->isRecvFromBelow:" << std::endl;
      for(int j = 0;j< this->isRecvFromBelow_.n();j++){
        statusOFS << "["<<j<<"] ";
        for(int i =0; i < this->isRecvFromBelow_.m();i++){
          statusOFS<< this->isRecvFromBelow_(i,j) << " ";
        }
        statusOFS<<std::endl;
      }
#endif

#ifndef _RELEASE_
      PopCallStack();
#endif


      TIMER_STOP(STR_RFL_RFB);




      TIMER_START(STCD_RFCD);


#ifndef _RELEASE_
      PushCallStack("SendToCrossDiagonal / RecvFromCrossDiagonal");
#endif
      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) ){
          for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, this->grid_ )].begin();
              si != localColBlockRowIdx[LBj( ksup, this->grid_ )].end(); si++ ){
            Int isup = *si;
            Int isupProcRow = PROW( isup, this->grid_ );
            Int isupProcCol = PCOL( isup, this->grid_ );
            if( isup > ksup && MYROW( this->grid_ ) == isupProcRow ){
              this->isSendToCrossDiagonal_(this->grid_->numProcCol, ksup ) = true;
              this->isSendToCrossDiagonal_(isupProcCol, ksup ) = true;
            }
          } // for (si)
        } // if( MYCOL( this->grid_ ) == PCOL( ksup, this->grid_ ) )
      } // for (ksup)

      for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
        if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) ){
          for( std::set<Int>::iterator si = localRowBlockColIdx[ LBi(ksup, this->grid_) ].begin();
              si != localRowBlockColIdx[ LBi(ksup, this->grid_) ].end(); si++ ){
            Int jsup = *si;
            Int jsupProcCol = PCOL( jsup, this->grid_ );
            Int jsupProcRow = PROW( jsup, this->grid_ );
            if( jsup > ksup && MYCOL(this->grid_) == jsupProcCol ){
              this->isRecvFromCrossDiagonal_(this->grid_->numProcRow, ksup ) = true;
              this->isRecvFromCrossDiagonal_(jsupProcRow, ksup ) = true;
            }
          } // for (si)
        } // if( MYROW( this->grid_ ) == PROW( ksup, this->grid_ ) )
      } // for (ksup)
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "this->isSendToCrossDiagonal:" << std::endl;
      for(int j =0; j < this->isSendToCrossDiagonal_.n();j++){
        if(this->isSendToCrossDiagonal_(this->grid_->numProcCol,j)){
          statusOFS << "["<<j<<"] ";
          for(int i =0; i < this->isSendToCrossDiagonal_.m()-1;i++){
            if(this->isSendToCrossDiagonal_(i,j))
            {
              statusOFS<< PNUM(PROW(j,this->grid_),i,this->grid_)<<" ";
            }
          }
          statusOFS<<std::endl;
        }
      }

      statusOFS << std::endl << "this->isRecvFromCrossDiagonal:" << std::endl;
      for(int j =0; j < this->isRecvFromCrossDiagonal_.n();j++){
        if(this->isRecvFromCrossDiagonal_(this->grid_->numProcRow,j)){
          statusOFS << "["<<j<<"] ";
          for(int i =0; i < this->isRecvFromCrossDiagonal_.m()-1;i++){
            if(this->isRecvFromCrossDiagonal_(i,j))
            {
              statusOFS<< PNUM(i,PCOL(j,this->grid_),this->grid_)<<" ";
            }
          }
          statusOFS<<std::endl;
        }
      }


#endif

#ifndef _RELEASE_
      PopCallStack();
#endif

      TIMER_STOP(STCD_RFCD);

#ifndef _RELEASE_
      PopCallStack();
#endif

      //Build the list of supernodes based on the elimination tree from SuperLU
      this->GetWorkSet(snodeEtree,this->WorkingSet());

      return ;
    } 		// -----  end of method PMatrix::ConstructCommunicationPattern_P2p  ----- 





} // namespace PEXSI

#endif //_PEXSI_PSELINV_ASYM_IMPL_HPP_

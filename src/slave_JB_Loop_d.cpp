#include <cstdio>
#include "pexsi/slave_kernel.hpp"

extern "C" {
#include <slave.h>
}

#ifndef MAX_BLOCK_SIZE 
#define MAX_BLOCK_SIZE 256
#endif

extern "C" void  JB_Loop_d(JB_Loop_param_d_t* param){

    LBlock_d_t* LcolRecvTmp = param->LcolRecvTmp;
    int LcolRecvTmpSize = param->LcolRecvTmpSize;
    UBlock_d_t* UrowRecvTmp = param->UrowRecvTmp;
    int UrowRecvTmpSize = param->UrowRecvTmpSize;
    LBlocks_d_t* LTmp = param->LTmp;
    int LTmpSize = param->LTmpSize;
    UBlocks_d_t* UTmp = param->UTmp;
    int UTmpSize = param->UTmpSize;
    Matrix_d_t* AinvBufTmp = param->AinvBufTmp;
    int* rowPtr = param->rowPtr;
    int rowPtrSize = param->rowPtrSize;
    int* colPtr = param->colPtr;
    int colPtrSize = param->colPtrSize;
    int numProcCol = param->numProcCol;
    int numProcRow = param->numProcRow;
    int* superPtr = param->superPtr;
    int superPtrSize = param->superPtrSize;

    int local_i_start = (LcolRecvTmpSize * _ROW)/8;
    int local_i_end = (LcolRecvTmpSize * (_ROW+1))/8;

    int local_j_start = (UrowRecvTmpSize * _COL)/8;
    int local_j_end = (UrowRecvTmpSize * (_COL+1))/8;

    for( int jb = local_j_start; jb < local_j_end; jb++ ){
        for( int ib = local_i_start; ib < local_i_end; ib++ ){
            LBlock_d_t* LB = &LcolRecvTmp[ib];
            UBlock_d_t* UB = &UrowRecvTmp[jb];
            int isup = LB->blockIdx;
            int jsup = UB->blockIdx;
            double* nzvalAinv = &AinvBufTmp->val[rowPtr[ib]+colPtr[jb]*AinvBufTmp->ld];
            int     ldAinv    = AinvBufTmp->m;

            // Pin down the corresponding block in the part of Sinv.
            if( isup >= jsup ){
            // std::vector<LBlock<double> >&  LcolSinv = this->L( LBj(jsup, grid_ ) );
            // LBlocks_d_t LcolSinv = LTmp[LBj(jsup, grid_)];
            LBlocks_d_t LcolSinv = LTmp[jsup/numProcCol];
            bool isBlockFound = false;
            for( int ibSinv = 0; ibSinv < LcolSinv.len; ibSinv++ ){
                // Found the (isup, jsup) block in Sinv
                if( LcolSinv.lblocks[ibSinv].blockIdx == isup ){
                LBlock_d_t* SinvB = &LcolSinv.lblocks[ibSinv];

                // Row relative indices
                int* rowsLBPtr    = LB->rows;
                int* rowsSinvBPtr = SinvB->rows;

                // Column relative indicies
                // int SinvColsSta = FirstBlockCol( jsup, super_ );
                int SinvColsSta = superPtr[jsup];

                int relCols[MAX_BLOCK_SIZE];
                for( int j = 0; j < UB->numCol; j++ ){
                    relCols[j] = UB->cols[j] - SinvColsSta;
                }

                int relRows[MAX_BLOCK_SIZE];
                for( int i = 0; i < LB->numRow; i++ ){
                    bool isRowFound = false;
                    for( int i1 = 0; i1 < SinvB->numRow; i1++ ){
                    if( rowsLBPtr[i] == rowsSinvBPtr[i1] ){
                        isRowFound = true;
                        relRows[i] = i1;
                        break;
                    }
                    }
                    if( isRowFound == false ){
                        printf("Row %d, in LB cannot find the corresponding row in SinvB\n",rowsLBPtr[i]);
                    }
                }

                // Transfer the values from Sinv to AinvBlock
                double* nzvalSinv = SinvB->nzval;
                int     ldSinv    = SinvB->numRow;
                for( int j = 0; j < UB->numCol; j++ ){
                    for( int i = 0; i < LB->numRow; i++ ){
                    nzvalAinv[i+j*ldAinv] =
                        nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                    }
                }

                isBlockFound = true;
                break;
                }	
            } // for (ibSinv )
            if( isBlockFound == false ){
                printf("Block(%d, %d) did not find a matching block in Sinv.",isup,jsup);
            }
            } // if (isup, jsup) is in L
            else{
            // UBlocks_d_t UrowSinv = UTmp[LBi( isup, grid_)];
            UBlocks_d_t UrowSinv = UTmp[isup/numProcRow];
            bool isBlockFound = false;
            for( int jbSinv = 0; jbSinv < UrowSinv.len; jbSinv++ ){
                // Found the (isup, jsup) block in Sinv
                if( UrowSinv.ublocks[jbSinv].blockIdx == jsup ){
                UBlock_d_t* SinvB = &UrowSinv.ublocks[jbSinv];

                // Row relative indices
                // int SinvRowsSta = FirstBlockCol( isup, super_ );
                int SinvRowsSta = superPtr[isup];

                int* colsUBPtr    = UB->cols;
                int* colsSinvBPtr = SinvB->cols;

                int relRows[MAX_BLOCK_SIZE];
                for( int i = 0; i < LB->numRow; i++ ){
                    relRows[i] = LB->rows[i] - SinvRowsSta;
                }

                int relCols[MAX_BLOCK_SIZE];
                // Column relative indices
                for( int j = 0; j < UB->numCol; j++ ){
                    bool isColFound = false;
                    for( int j1 = 0; j1 < SinvB->numCol; j1++ ){
                    if( colsUBPtr[j] == colsSinvBPtr[j1] ){
                        isColFound = true;
                        relCols[j] = j1;
                        break;
                    }
                    }
                    if( isColFound == false ){
                        printf("Col %d, in UB cannot find the corresponding row in SinvB\n",colsUBPtr[j]);
                    }
                }

                // Transfer the values from Sinv to AinvBlock
                double* nzvalSinv = SinvB->nzval;
                int     ldSinv    = SinvB->numRow;
                for( int j = 0; j < UB->numCol; j++ ){
                    for( int i = 0; i < LB->numRow; i++ ){
                    nzvalAinv[i+j*ldAinv] =
                        nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                    }
                }

                isBlockFound = true;
                break;
                }
            } // for (jbSinv)
            if( isBlockFound == false ){
                printf("Block(%d, %d) did not find a matching block in Sinv.\n",isup,jsup);
            }
            } // if (isup, jsup) is in U

        } // for( ib )
    } // for ( jb )
}
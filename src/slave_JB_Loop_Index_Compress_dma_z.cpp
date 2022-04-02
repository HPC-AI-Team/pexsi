#include <cstdio>
#include <complex>
#include "pexsi/slave_kernel.hpp"

extern "C" {
#include <slave.h>
}

#ifndef MAX_BLOCK_SIZE 
#define MAX_BLOCK_SIZE 256
#endif

typedef struct{
    int segment_ptr[MAX_BLOCK_SIZE+1];
    int segment_offset[MAX_BLOCK_SIZE];
    int segment_count;
} indirect_index_segment_compress_t;

void indirect_index_segment_compress_init(indirect_index_segment_compress_t* segment_compress,const int *indirect_index, int len){
    if(len < 1){
        segment_compress->segment_count = 0;
        return;
    }

    int* segment_ptr = segment_compress->segment_ptr;
    int* segment_offset = segment_compress->segment_offset;

    // segment compress
    segment_ptr[0] = 0;
    segment_offset[0] = indirect_index[0];
    int segment_prev = indirect_index[0]; 
    int segment_count = 1;
    for(int i = 1; i < len; i++){
        int segment_cur = indirect_index[i] - i;
        if(segment_prev != segment_cur){
            segment_ptr[segment_count] = i;
            segment_offset[segment_count] = segment_cur;
            segment_count += 1;
            segment_prev = segment_cur;
        }
    }
    segment_ptr[segment_count] = len;

    // output
    segment_compress->segment_count = segment_count;
}

static void indirect_index_segment_compress_destroy(indirect_index_segment_compress_t* segment_compress){
}

extern "C" void  JB_Loop_index_compress_dma_z(JB_Loop_param_z_t* param){

    LBlock_z_t* LcolRecvTmp = param->LcolRecvTmp;
    int LcolRecvTmpSize = param->LcolRecvTmpSize;
    UBlock_z_t* UrowRecvTmp = param->UrowRecvTmp;
    int UrowRecvTmpSize = param->UrowRecvTmpSize;
    LBlocks_z_t* LTmp = param->LTmp;
    int LTmpSize = param->LTmpSize;
    UBlocks_z_t* UTmp = param->UTmp;
    int UTmpSize = param->UTmpSize;
    Matrix_z_t* AinvBufTmp = param->AinvBufTmp;
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
            LBlock_z_t* LB = &LcolRecvTmp[ib];
            UBlock_z_t* UB = &UrowRecvTmp[jb];
            int isup = LB->blockIdx;
            int jsup = UB->blockIdx;
            std::complex<double>* nzvalAinv = &AinvBufTmp->val[rowPtr[ib]+colPtr[jb]*AinvBufTmp->ld];
            int     ldAinv    = AinvBufTmp->m;

            // Pin down the corresponding block in the part of Sinv.
            if( isup >= jsup ){
            // std::vector<LBlock<std::complex<double>> >&  LcolSinv = this->L( LBj(jsup, grid_ ) );
            // LBlocks_z_t LcolSinv = LTmp[LBj(jsup, grid_)];
            LBlocks_z_t LcolSinv = LTmp[jsup/numProcCol];
            bool isBlockFound = false;
            for( int ibSinv = 0; ibSinv < LcolSinv.len; ibSinv++ ){
                // Found the (isup, jsup) block in Sinv
                if( LcolSinv.lblocks[ibSinv].blockIdx == isup ){
                LBlock_z_t* SinvB = &LcolSinv.lblocks[ibSinv];

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

                indirect_index_segment_compress_t segment_compress;
                indirect_index_segment_compress_init(&segment_compress,relRows,LB->numRow);

                std::complex<double> buffer[MAX_BLOCK_SIZE];

                // Transfer the values from Sinv to AinvBlock
                std::complex<double>* nzvalSinv = SinvB->nzval;
                int ldSinv    = SinvB->numRow;
                for( int j = 0; j < UB->numCol; j++ ){
                  std::complex<double>* nzvalAinv_j = &nzvalAinv[j*ldAinv];
                  std::complex<double>* nzvalSinv_j = &nzvalSinv[relCols[j] * ldSinv];
                  // for( int i = 0; i < LB.numRow; i++ ){
                  //   nzvalAinv_j[i] = nzvalSinv_j[relRows[i]];
                  // }
                  for(int ptr = 0; ptr < segment_compress.segment_count; ++ptr){
                      int i_start = segment_compress.segment_ptr[ptr];
                      int i_end = segment_compress.segment_ptr[ptr+1];
                      int offset = segment_compress.segment_offset[ptr];
                      std::complex<double>* nzvalSinv_j_offset = nzvalSinv_j + offset;
                      // #pragma omp simd
                      // for(int i = i_start; i < i_end; i++){
                      //     nzvalAinv_j[i] = nzvalSinv_j_offset[i];
                      // } 
                      athread_dma_get(buffer, nzvalSinv_j_offset + i_start, (i_end - i_start) * sizeof(std::complex<double>));
                      athread_dma_put(nzvalAinv_j + i_start, buffer, (i_end - i_start) * sizeof(std::complex<double>));
                    //   memcpy(nzvalAinv_j + i_start,nzvalSinv_j_offset + i_start,(i_end - i_start) * sizeof(std::complex<double>));
                  }
                }               

                // // Transfer the values from Sinv to AinvBlock
                // std::complex<double>* nzvalSinv = SinvB->nzval;
                // int     ldSinv    = SinvB->numRow;
                // for( int j = 0; j < UB->numCol; j++ ){
                //     for( int i = 0; i < LB->numRow; i++ ){
                //     nzvalAinv[i+j*ldAinv] =
                //         nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                //     }
                // }

                isBlockFound = true;
                break;
                }	
            } // for (ibSinv )
            if( isBlockFound == false ){
                printf("Block(%d, %d) did not find a matching block in Sinv.",isup,jsup);
            }
            } // if (isup, jsup) is in L
            else{
            // UBlocks_z_t UrowSinv = UTmp[LBi( isup, grid_)];
            UBlocks_z_t UrowSinv = UTmp[isup/numProcRow];
            bool isBlockFound = false;
            for( int jbSinv = 0; jbSinv < UrowSinv.len; jbSinv++ ){
                // Found the (isup, jsup) block in Sinv
                if( UrowSinv.ublocks[jbSinv].blockIdx == jsup ){
                UBlock_z_t* SinvB = &UrowSinv.ublocks[jbSinv];

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

                indirect_index_segment_compress_t segment_compress;
                indirect_index_segment_compress_init(&segment_compress,relRows,LB->numRow);
                
                std::complex<double> buffer[MAX_BLOCK_SIZE];

                // Transfer the values from Sinv to AinvBlock
                std::complex<double>* nzvalSinv = SinvB->nzval;
                int ldSinv    = SinvB->numRow;
                for( int j = 0; j < UB->numCol; j++ ){
                  std::complex<double>* nzvalAinv_j = &nzvalAinv[j*ldAinv];
                  std::complex<double>* nzvalSinv_j = &nzvalSinv[relCols[j] * ldSinv];
                  // for( int i = 0; i < LB.numRow; i++ ){
                  //   nzvalAinv_j[i] = nzvalSinv_j[relRows[i]];
                  // }
                  for(int ptr = 0; ptr < segment_compress.segment_count; ++ptr){
                      int i_start = segment_compress.segment_ptr[ptr];
                      int i_end = segment_compress.segment_ptr[ptr+1];
                      int offset = segment_compress.segment_offset[ptr];
                      std::complex<double>* nzvalSinv_j_offset = nzvalSinv_j + offset;
                      // #pragma omp simd
                      // for(int i = i_start; i < i_end; i++){
                      //     nzvalAinv_j[i] = nzvalSinv_j_offset[i];
                      // } 
                      // memcpy(nzvalAinv_j + i_start,nzvalSinv_j_offset + i_start,(i_end - i_start) * sizeof(std::complex<double>));
                      athread_dma_get(buffer, nzvalSinv_j_offset + i_start, (i_end - i_start) * sizeof(std::complex<double>));
                      athread_dma_put(nzvalAinv_j + i_start, buffer, (i_end - i_start) * sizeof(std::complex<double>));
                  }
                }

                indirect_index_segment_compress_destroy(&segment_compress);

                // // Transfer the values from Sinv to AinvBlock
                // std::complex<double>* nzvalSinv = SinvB->nzval;
                // int     ldSinv    = SinvB->numRow;
                // for( int j = 0; j < UB->numCol; j++ ){
                //     for( int i = 0; i < LB->numRow; i++ ){
                //     nzvalAinv[i+j*ldAinv] =
                //         nzvalSinv[relRows[i] + relCols[j] * ldSinv];
                //     }
                // }

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
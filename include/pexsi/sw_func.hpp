#pragma once
#include "pexsi/sw_def.hpp"
#include <complex>
#include <vector>

using PEXSI::LBlock;
using PEXSI::UBlock;
using PEXSI::NumMat;
using std::complex;
using std::vector;

void LBlock_copy_d(PEXSI::LBlock<double>& src,LBlock_d_t* dest){
    dest->blockIdx = src.blockIdx;
    dest->numRow = src.numRow;
    dest->numCol = src.numCol;
    dest->rows = src.rows.Data();
    dest->rows_size = src.rows.m();
    dest->nzval = src.nzval.Data();
    dest->nzval_m = src.nzval.m();
    dest->nzval_n = src.nzval.n();
    dest->nzval_bufsize = src.nzval.AllocatedSize();
}

void UBlock_copy_d(PEXSI::UBlock<double>& src,UBlock_d_t* dest){
    dest->blockIdx = src.blockIdx;
    dest->numRow = src.numRow;
    dest->numCol = src.numCol;
    dest->cols = src.cols.Data();
    dest->cols_size = src.cols.m();
    dest->nzval = src.nzval.Data();
    dest->nzval_m = src.nzval.m();
    dest->nzval_n = src.nzval.n();
    dest->nzval_bufsize = src.nzval.AllocatedSize();
}

void LBlock_vector_d_init(std::vector<PEXSI::LBlock<double>>& src,LBlock_d_t** dest, int* len){
    int size = src.size();
    *len = size;
    LBlock_d_t* lblock_ptr = (LBlock_d_t*)malloc(size * sizeof(LBlock_d_t));
    *dest = lblock_ptr;
    for(int i = 0; i < size; ++i){
        LBlock_copy_d(src[i],&lblock_ptr[i]);
    }
}

void UBlock_vector_d_init(std::vector<PEXSI::UBlock<double>>& src,UBlock_d_t** dest, int* len){
    int size = src.size();
    *len = size;
    UBlock_d_t* ublock_ptr = (UBlock_d_t*)malloc(size * sizeof(UBlock_d_t));
    *dest = ublock_ptr;
    for(int i = 0; i < size; ++i){
        UBlock_copy_d(src[i],&ublock_ptr[i]);
    }
}

void LBlock_vector_d_destory(LBlock_d_t* dest){
    free(dest);
}

void UBlock_vector_d_destory(UBlock_d_t* dest){
    free(dest);
}

void LBlock_matrix_d_init(vector<vector<LBlock<double>>>& src,LBlocks_d_t** dest, int* len){
    int size = src.size();
    *len = size;
    LBlocks_d_t* lblocks_ptr = (LBlocks_d_t*)malloc(size * sizeof(LBlocks_d_t));
    *dest = lblocks_ptr;
    for(int i = 0; i < size; ++i){
        LBlock_d_t* lblocks;
        int len;
        LBlock_vector_d_init(src[i], &lblocks, &len);
        lblocks_ptr[i].lblocks = lblocks;
        lblocks_ptr[i].len = len;
    }
}

void LBlock_matrix_d_destroy(LBlocks_d_t* dest, int len){
    for(int i = 0; i< len;++i){
        LBlock_vector_d_destory(dest[i].lblocks);
    }
}

void UBlock_matrix_d_init(vector<vector<UBlock<double>>>& src,UBlocks_d_t** dest, int* len){
    int size = src.size();
    *len = size;
    UBlocks_d_t* ublocks_ptr = (UBlocks_d_t*)malloc(size * sizeof(UBlocks_d_t));
    *dest = ublocks_ptr;
    for(int i = 0; i < size; ++i){
        UBlock_d_t* ublocks;
        int len;
        UBlock_vector_d_init(src[i], &ublocks, &len);
        ublocks_ptr[i].ublocks = ublocks;
        ublocks_ptr[i].len = len;
    }
}

void UBlock_matrix_d_destroy(UBlocks_d_t* dest, int len){
    for(int i = 0; i< len;++i){
        UBlock_vector_d_destory(dest[i].ublocks);
    }
}

void Matrix_d_init(NumMat<double>& src,Matrix_d_t* dest){
    dest->val = src.Data();
    dest->m = src.m();
    dest->n = src.n();
    dest->ld = src.m();
}

// z ----------------------------------------------------------------------------------------------

void LBlock_copy_z(LBlock<complex<double>>& src,LBlock_z_t* dest){
    dest->blockIdx = src.blockIdx;
    dest->numRow = src.numRow;
    dest->numCol = src.numCol;
    dest->rows = src.rows.Data();
    dest->rows_size = src.rows.m();
    dest->nzval = src.nzval.Data();
    dest->nzval_m = src.nzval.m();
    dest->nzval_n = src.nzval.n();
    dest->nzval_bufsize = src.nzval.AllocatedSize();
}

void UBlock_copy_z(UBlock<complex<double>>& src,UBlock_z_t* dest){
    dest->blockIdx = src.blockIdx;
    dest->numRow = src.numRow;
    dest->numCol = src.numCol;
    dest->cols = src.cols.Data();
    dest->cols_size = src.cols.m();
    dest->nzval = src.nzval.Data();
    dest->nzval_m = src.nzval.m();
    dest->nzval_n = src.nzval.n();
    dest->nzval_bufsize = src.nzval.AllocatedSize();
}

void LBlock_vector_z_init(vector<LBlock<complex<double>>>& src,LBlock_z_t** dest, int* len){
    int size = src.size();
    *len = size;
    LBlock_z_t* lblock_ptr = (LBlock_z_t*)malloc(size * sizeof(LBlock_z_t));
    *dest = lblock_ptr;
    for(int i = 0; i < size; ++i){
        LBlock_copy_z(src[i],&lblock_ptr[i]);
    }
}

void UBlock_vector_z_init(vector<UBlock<complex<double>>>& src,UBlock_z_t** dest, int* len){
    int size = src.size();
    *len = size;
    UBlock_z_t* ublock_ptr = (UBlock_z_t*)malloc(size * sizeof(UBlock_z_t));
    *dest = ublock_ptr;
    for(int i = 0; i < size; ++i){
        UBlock_copy_z(src[i],&ublock_ptr[i]);
    }
}

void LBlock_vector_z_destory(LBlock_z_t* dest){
    free(dest);
}

void UBlock_vector_z_destory(UBlock_z_t* dest){
    free(dest);
}

void LBlock_matrix_z_init(vector<vector<LBlock<complex<double>>>>& src,LBlocks_z_t** dest, int* len){
    int size = src.size();
    *len = size;
    LBlocks_z_t* lblocks_ptr = (LBlocks_z_t*)malloc(size * sizeof(LBlocks_z_t));
    *dest = lblocks_ptr;
    for(int i = 0; i < size; ++i){
        LBlock_z_t* lblocks;
        int len;
        LBlock_vector_z_init(src[i], &lblocks, &len);
        lblocks_ptr[i].lblocks = lblocks;
        lblocks_ptr[i].len = len;
    }
}

void LBlock_matrix_z_destroy(LBlocks_z_t* dest, int len){
    for(int i = 0; i< len;++i){
        LBlock_vector_z_destory(dest[i].lblocks);
    }
}

void UBlock_matrix_z_init(vector<vector<UBlock<complex<double>>>>& src,UBlocks_z_t** dest, int* len){
    int size = src.size();
    *len = size;
    UBlocks_z_t* ublocks_ptr = (UBlocks_z_t*)malloc(size * sizeof(UBlocks_z_t));
    *dest = ublocks_ptr;
    for(int i = 0; i < size; ++i){
        UBlock_z_t* ublocks;
        int len;
        UBlock_vector_z_init(src[i], &ublocks, &len);
        ublocks_ptr[i].ublocks = ublocks;
        ublocks_ptr[i].len = len;
    }
}

void UBlock_matrix_z_destroy(UBlocks_z_t* dest, int len){
    for(int i = 0; i< len;++i){
        UBlock_vector_z_destory(dest[i].ublocks);
    }
}

void Matrix_z_init(NumMat<complex<double>>& src,Matrix_z_t* dest){
    dest->val = src.Data();
    dest->m = src.m();
    dest->n = src.n();
    dest->ld = src.m();
}

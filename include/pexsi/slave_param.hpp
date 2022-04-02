#pragma once

#include "pexsi/sw_def.hpp"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int* ids;
} hello_world_param_t;

typedef struct {
    LBlock_z_t* LcolRecvTmp;
    int LcolRecvTmpSize;
    UBlock_z_t* UrowRecvTmp;
    int UrowRecvTmpSize;
    LBlocks_z_t* LTmp;
    int LTmpSize;
    UBlocks_z_t* UTmp;
    int UTmpSize;
    Matrix_z_t* AinvBufTmp;
    int* rowPtr;
    int rowPtrSize;
    int* colPtr;
    int colPtrSize;
    int numProcCol;
    int numProcRow;
    int* superPtr;
    int superPtrSize;
} JB_Loop_param_z_t;

typedef struct {
    LBlock_d_t* LcolRecvTmp;
    int LcolRecvTmpSize;
    UBlock_d_t* UrowRecvTmp;
    int UrowRecvTmpSize;
    LBlocks_d_t* LTmp;
    int LTmpSize;
    UBlocks_d_t* UTmp;
    int UTmpSize;
    Matrix_d_t* AinvBufTmp;
    int* rowPtr;
    int rowPtrSize;
    int* colPtr;
    int colPtrSize;
    int numProcCol;
    int numProcRow;
    int* superPtr;
    int superPtrSize;
} JB_Loop_param_d_t;

#ifdef __cplusplus
}
#endif
#pragma once

#include <cstdlib>
#include <complex>

using std::complex;

typedef struct{
    int blockIdx;
    int numRow;
    int numCol;
    int* rows;
    int rows_size;
    double* nzval;
    int nzval_m;
    int nzval_n;
    int nzval_bufsize;
} LBlock_d_t;

typedef struct{
    int blockIdx;
    int numRow;
    int numCol;
    int* cols;
    int cols_size;
    double* nzval;
    int nzval_m;
    int nzval_n;
    int nzval_bufsize;
} UBlock_d_t;

typedef struct{
    LBlock_d_t* lblocks;
    int len;
}LBlocks_d_t;

typedef struct{
    UBlock_d_t* ublocks;
    int len;
}UBlocks_d_t;

typedef struct{
    double* val;
    int m;
    int n;
    int ld;
} Matrix_d_t;

// z ----------------------------------------------------------------------------------------------

typedef struct{
    int blockIdx;
    int numRow;
    int numCol;
    int* rows;
    int rows_size;
    complex<double>* nzval;
    int nzval_m;
    int nzval_n;
    int nzval_bufsize;
} LBlock_z_t;

typedef struct{
    int blockIdx;
    int numRow;
    int numCol;
    int* cols;
    int cols_size;
    complex<double>* nzval;
    int nzval_m;
    int nzval_n;
    int nzval_bufsize;
} UBlock_z_t;

typedef struct{
    LBlock_z_t* lblocks;
    int len;
}LBlocks_z_t;

typedef struct{
    UBlock_z_t* ublocks;
    int len;
}UBlocks_z_t;

typedef struct{
    complex<double>* val;
    int m;
    int n;
    int ld;
} Matrix_z_t;
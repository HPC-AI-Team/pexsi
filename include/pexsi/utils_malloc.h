#pragma once

#ifdef __cplusplus
extern "C" { 
#endif 

#include <stdlib.h>

#ifndef MAX_BLOCK_SIZE 
#define MAX_BLOCK_SIZE 128
#endif

int get_superlu_env_nsup(){
    char* tmp = getenv("NSUP");
    if(tmp)
        return atoi(tmp);
    return 128;
}

int get_superlu_env_nrel(){
    char* tmp = getenv("NREL");
    if(tmp)
        return atoi(tmp);
    return 20;
}


typedef struct{
    int* segment_ptr;
    int* segment_offset;
    int segment_count;
} indirect_index_segment_compress_t;

static void indirect_index_segment_compress_init(indirect_index_segment_compress_t* segment_compress,const int *indirect_index, int len){
    if(len < 1){
        segment_compress->segment_count = 0;
        return;
    }
    int nsup = get_superlu_env_nsup();
    // count segment
    // int segment_count = 1;
    // int segment_prev = indirect_index[0];
    // for(int i = 1; i < len; i++){
    //     int segment_cur = indirect_index[i] - i;
    //     if(segment_prev != segment_cur){
    //         segment_count += 1;
    //         segment_prev = segment_cur;
    //     }
    // }
    // allocate space
    int* segment_ptr = (int*)malloc((nsup+1) * sizeof(int));
    int* segment_offset = (int*)malloc(nsup * sizeof(int));

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
    segment_compress->segment_ptr = segment_ptr;
    segment_compress->segment_offset = segment_offset;

}

static void indirect_index_segment_compress_destroy(indirect_index_segment_compress_t* segment_compress){
    free(segment_compress->segment_ptr);
    free(segment_compress->segment_offset);
}


void indirect_index_segment_compress_use_example(){
    indirect_index_segment_compress_t segment_compress;
    double* arr;
    for(int ptr = 0; ptr < segment_compress.segment_count; ++ptr){
        int i_start = segment_compress.segment_ptr[ptr];
        int i_end = segment_compress.segment_ptr[ptr+1];
        int offset = segment_compress.segment_offset[ptr];
        double *ARR = arr + offset;
        for(int i = i_start; i < i_end; i++){
            // do some thing for ARR[i]
        } 
    }
}

















#ifdef __cplusplus
} 
#endif 
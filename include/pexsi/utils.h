#pragma once

#ifdef __cplusplus
extern "C" { 
#endif 

#include <stdlib.h>

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
    // count segment
    int segment_count = 1;
    int segment_prev = indirect_index[0];
    for(int i = 1; i < len; i++){
        int segment_cur = indirect_index[i] - i;
        if(segment_prev != segment_cur){
            segment_count += 1;
            segment_prev = segment_cur;
        }
    }
    // allocate space
    int* segment_ptr = (int*)malloc((segment_count+1) * sizeof(int));
    int* segment_offset = (int*)malloc(segment_count * sizeof(int));

    // segment compress
    segment_ptr[0] = 0;
    segment_offset[0] = indirect_index[0];
    segment_prev = indirect_index[0]; 
    int segment_index = 1;
    for(int i = 1; i < len; i++){
        int segment_cur = indirect_index[i] - i;
        if(segment_prev != segment_cur){
            segment_ptr[segment_index] = i;
            segment_offset[segment_index] = segment_cur;
            segment_index += 1;
            segment_prev = segment_cur;
        }
    }
    segment_ptr[segment_index] = len;

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
#pragma once

#include "pexsi/slave_param.hpp"

#ifdef __cplusplus
extern "C" {
#endif

extern void slave_hello_world();
extern void slave_JB_Loop_z();
extern void slave_JB_Loop_index_compress_z();
extern void slave_JB_Loop_index_compress_dma_z();

extern void slave_JB_Loop_d();
extern void slave_JB_Loop_index_compress_d();
extern void slave_JB_Loop_index_compress_dma_d();

#ifdef __cplusplus
}
#endif
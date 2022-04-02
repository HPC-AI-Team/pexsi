#include <cstdio>
#include "pexsi/slave_kernel.hpp"

extern "C" {
#include <slave.h>
}


extern "C" void  hello_world(hello_world_param_t* param){
    for(int id = 0;id < 64;id++){
        if(_MYID == id){
            printf("MYID : %d\n",_MYID);
            printf("PEN : %d\n",_PEN);
            printf("ROW : %d\n",_ROW);
            printf("COL : %d\n",_COL);
            printf("CGN : %d\n",_CGN);
            printf("cache_size : %d\n",_cache_size);
            printf("ldm_share_mode : %d\n",_ldm_share_mode);
            printf("ldm_share_size : %d\n",_ldm_share_size);
            printf("allocatable_size : %ld\n",get_allocatable_size());
        }
        athread_ssync_array();
    }
}

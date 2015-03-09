#! /bin/bash

nproc=$2

cd $1
ls comm_stat* | grep -v tgz | ../parse_avg 0 $nproc 1  > bcastU_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | ../parse_avg 0 $nproc 3  > bcastL_total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | ../parse_avg 0 $nproc 4  > reduceL_total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | ../parse_avg 0 $nproc  > total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | ../parse_avg 1 $nproc  > avg.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | ../parse_avg 0 $nproc 9  > send_L_CD_total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | ../parse_avg_send 0 $nproc  > sender_total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | ../parse_avg_send 0 $nproc 1  > sender_bcastU_total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | ../parse_avg_send 0 $nproc 4  > sender_reduceL_total.dat 2>/dev/null
cd ..

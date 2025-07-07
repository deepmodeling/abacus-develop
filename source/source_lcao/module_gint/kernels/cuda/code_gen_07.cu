#include "vbatch_matrix_mul.cuh"

template void gemm_time_measure<double, 8, 32, 32, 64, 8, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 32, 32, 64, 16, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 32, 32, 64, 24, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 32, 32, 64, 32, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 32, 64, 64, 8, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 32, 64, 64, 16, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 32, 64, 64, 24, 8, 32, 8, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 24, 24, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 24, 32, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 24, 40, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 24, 48, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 24, 56, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 48, 16, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 48, 24, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 8, 48, 32, 12, 12, 8, 12, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 16, 48, 32, 12, 12, 16, 12, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 16, 48, 32, 24, 12, 16, 12, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 16, 48, 48, 12, 12, 16, 12, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 16, 48, 64, 12, 12, 16, 12, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 24, 48, 48, 12, 12, 24, 12, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 12, 24, 48, 48, 24, 12, 24, 12, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 16, 4, 32, 12, 16, 16, 4, 16, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 16, 4, 32, 16, 16, 16, 4, 16, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);


#include "vbatch_matrix_mul.cuh"

template void gemm_time_measure<double, 8, 8, 48, 40, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 56, 16, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 56, 24, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 56, 32, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 64, 16, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 64, 24, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 24, 24, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 24, 24, 16, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 24, 36, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 24, 36, 16, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 24, 48, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 24, 60, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 48, 24, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 48, 36, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 48, 48, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 12, 48, 60, 8, 8, 12, 8, 12>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 16, 48, 8, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 16, 48, 16, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 16, 48, 24, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 16, 64, 8, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 16, 64, 16, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 32, 32, 8, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 16, 32, 32, 16, 8, 16, 8, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);


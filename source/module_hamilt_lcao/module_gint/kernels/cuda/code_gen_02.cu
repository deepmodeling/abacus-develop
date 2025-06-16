#include "vbatch_matrix_mul.cuh"

template void gemm_time_measure<double, 4, 16, 16, 48, 12, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 16, 64, 4, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 16, 64, 8, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 32, 32, 4, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 32, 32, 8, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 32, 32, 12, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 32, 48, 4, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 32, 48, 8, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 48, 32, 4, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 16, 48, 32, 8, 4, 16, 4, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 24, 24, 48, 4, 4, 24, 4, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 24, 24, 48, 8, 4, 24, 4, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 24, 24, 48, 12, 4, 24, 4, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 24, 48, 48, 4, 4, 24, 4, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 24, 48, 48, 8, 4, 24, 4, 24>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 32, 32, 64, 4, 4, 32, 4, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 32, 32, 64, 8, 4, 32, 4, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 32, 32, 64, 12, 4, 32, 4, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 32, 32, 64, 16, 4, 32, 4, 32>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 6, 16, 48, 32, 6, 6, 16, 6, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 6, 16, 48, 32, 12, 6, 16, 6, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 6, 16, 48, 48, 6, 6, 16, 6, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 16, 12, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);


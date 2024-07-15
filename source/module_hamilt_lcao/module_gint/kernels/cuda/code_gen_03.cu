#include "vbatch_matrix_mul.cuh"

template void gemm_time_measure<double, 8, 4, 16, 16, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 16, 20, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 16, 24, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 16, 28, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 16, 32, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 24, 8, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 24, 12, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 24, 16, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 24, 20, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 24, 24, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 32, 8, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 32, 12, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 32, 16, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 4, 40, 8, 8, 8, 4, 8, 4>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 24, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 24, 16, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 32, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 32, 16, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 40, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 48, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 56, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 16, 64, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 8, 8, 24, 16, 8, 8, 8, 8, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);


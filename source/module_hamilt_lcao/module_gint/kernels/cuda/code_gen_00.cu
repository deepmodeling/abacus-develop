#include "vbatch_matrix_mul.cuh"

template void gemm_time_measure<double, 2, 16, 16, 32, 2, 2, 16, 2, 16>(int, int, int*, int*, int*, double**, int*, double**, int*, double**, int*, int, cudaStream_t, float&, matrix_multiple_func_type&, double*, double*, double*);

template void gemm_time_measure<double, 2, 16, 16, 32, 4, 2, 16, 2, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 2, 16, 16, 32, 6, 2, 16, 2, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 2, 16, 16, 32, 8, 2, 16, 2, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 2, 16, 16, 48, 2, 2, 16, 2, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 2, 16, 16, 48, 4, 2, 16, 2, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 2, 16, 16, 48, 6, 2, 16, 2, 16>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 24, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 24, 8, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 24, 12, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 32, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 32, 8, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 40, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 40, 8, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 48, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 56, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 8, 64, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 16, 16, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 16, 16, 8, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 16, 16, 12, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 16, 24, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 16, 24, 8, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);

template void gemm_time_measure<double, 4, 8, 16, 32, 4, 4, 8, 4, 8>(int,int,int*,int*,int*,double**,int*,double**,int*,double**,int*,int,cudaStream_t,float&,matrix_multiple_func_type&,double*,double*,double*);


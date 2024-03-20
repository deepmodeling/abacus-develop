        gemm_time_measure<double, 2,16,16,32,2,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 2,16,16,32,4,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 2,16,16,32,6,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 2,16,16,32,8,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 2,16,16,48,2,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 2,16,16,48,4,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 2,16,16,48,6,2,16,2,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,24,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,24,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,24,12,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,32,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,32,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,40,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,40,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,48,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,56,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,8,64,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,16,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,16,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,16,12,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,24,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,24,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,32,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,32,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,40,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,48,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,16,56,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,24,16,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,24,16,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,24,24,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,24,24,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,24,32,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,24,40,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,32,16,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,32,16,8,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,32,24,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,40,16,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,40,24,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,48,16,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,8,56,16,4,4,8,4,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,32,4,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,32,8,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,32,12,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,32,16,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,48,4,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,48,8,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,48,12,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,64,4,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,16,64,8,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,32,32,4,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,32,32,8,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,32,32,12,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,32,48,4,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,32,48,8,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,48,32,4,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,16,48,32,8,4,16,4,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,24,24,48,4,4,24,4,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,24,24,48,8,4,24,4,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,24,24,48,12,4,24,4,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,24,48,48,4,4,24,4,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,24,48,48,8,4,24,4,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,32,32,64,4,4,32,4,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,32,32,64,8,4,32,4,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,32,32,64,12,4,32,4,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 4,32,32,64,16,4,32,4,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 6,16,48,32,6,6,16,6,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 6,16,48,32,12,6,16,6,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 6,16,48,48,6,6,16,6,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,16,12,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,16,16,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,16,20,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,16,24,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,16,28,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,16,32,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,24,8,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,24,12,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,24,16,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,24,20,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,24,24,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,32,8,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,32,12,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,32,16,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,4,40,8,8,8,4,8,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,24,16,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,32,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,32,16,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,40,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,48,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,56,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,16,64,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,16,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,16,16,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,24,16,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,32,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,40,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,48,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,56,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,24,64,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,16,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,16,16,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,32,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,40,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,48,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,32,56,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,40,16,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,40,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,40,32,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,40,40,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,40,48,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,48,16,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,48,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,48,32,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,48,40,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,56,16,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,56,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,56,32,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,64,16,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,8,64,24,8,8,8,8,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,24,24,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,24,24,16,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,24,36,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,24,36,16,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,24,48,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,24,60,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,48,24,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,48,36,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,48,48,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,12,48,60,8,8,12,8,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,16,48,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,16,48,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,16,48,24,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,16,64,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,16,64,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,32,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,32,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,32,24,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,48,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,48,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,64,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,32,64,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,48,32,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,48,32,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,48,48,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,48,48,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,48,64,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,64,32,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,64,32,16,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,16,64,48,8,8,16,8,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,20,40,40,8,8,20,8,20>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,20,40,40,16,8,20,8,20>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,20,40,60,8,8,20,8,20>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,24,24,48,8,8,24,8,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,24,24,48,16,8,24,8,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,24,24,48,24,8,24,8,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,24,48,48,8,8,24,8,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,24,48,48,16,8,24,8,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,28,56,56,8,8,28,8,28>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,28,56,56,16,8,28,8,28>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,32,64,8,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,32,64,16,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,32,64,24,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,32,64,32,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,64,64,8,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,64,64,16,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 8,32,64,64,24,8,32,8,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,24,24,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,24,32,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,24,40,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,24,48,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,24,56,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,48,16,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,48,24,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,8,48,32,12,12,8,12,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,16,48,32,12,12,16,12,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,16,48,32,24,12,16,12,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,16,48,48,12,12,16,12,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,16,48,64,12,12,16,12,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,24,48,48,12,12,24,12,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 12,24,48,48,24,12,24,12,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,4,32,12,16,16,4,16,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,4,32,16,16,16,4,16,4>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,6,48,12,16,16,6,16,6>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,32,24,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,32,32,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,32,40,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,32,48,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,32,56,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,32,64,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,48,16,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,48,24,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,48,32,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,48,40,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,48,48,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,64,16,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,64,24,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,8,64,32,16,16,8,16,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,12,48,24,16,16,12,16,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,12,48,36,16,16,12,16,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,12,48,48,16,16,12,16,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,12,48,60,16,16,12,16,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,32,48,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,32,48,32,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,32,64,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,32,64,32,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,48,32,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,48,32,32,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,48,48,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,48,48,32,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,48,64,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,64,32,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,64,32,32,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,64,48,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,16,64,64,16,16,16,16,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,24,48,48,16,16,24,16,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,24,48,48,32,16,24,16,24>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,32,64,64,16,16,32,16,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,32,64,64,32,16,32,16,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 16,32,64,64,32,16,32,16,32>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 20,8,40,24,20,20,8,20,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 20,8,40,32,20,20,8,20,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,8,48,24,24,24,8,24,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,8,48,32,24,24,8,24,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,12,48,36,24,24,12,24,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,12,48,48,24,24,12,24,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,12,48,60,24,24,12,24,12>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,16,48,48,24,24,16,24,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 24,16,48,64,24,24,16,24,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 32,8,64,24,32,32,8,32,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 32,8,64,32,32,32,8,32,8>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 32,16,64,48,32,32,16,32,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);

        gemm_time_measure<double, 32,16,64,64,32,32,16,32,16>(max_m, max_n, d_m, d_n, d_k,
                                                              d_global_A_array, d_global_lda,
                                                              d_global_B_array, d_global_ldb,
                                                              d_global_C_array, d_global_ldc,
                                                              batchCount, temp_stream, fastest_time, fastest_algo, cpu_result, h_global_C, d_global_C);


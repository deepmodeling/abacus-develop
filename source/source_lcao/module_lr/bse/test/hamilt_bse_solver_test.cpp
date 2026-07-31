#include <gtest/gtest.h>
#include "../hamilt_bse_solver.h"
#include <mpi.h>
#include "source_base/module_container/base/third_party/blas.h"

#define rand01 (static_cast<double>(rand()) / static_cast<double>(RAND_MAX) - 0.5 ) // [-0.5, 0.5]

std::vector<std::complex<double>> generate_conjugate_matrix(int n) {
    std::vector<std::complex<double>> matrix(n*n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            if(i==j) matrix[i * n + j] = std::complex<double>(20+rand01, 0); // gaurantee {{A,B},{A*,B*}} is positive definite
            else{
                matrix[i * n + j] = std::complex<double>(rand01, rand01);
                matrix[j * n + i] = std::conj(matrix[i * n + j]);
            }
        }
    }
    return matrix;
}
std::vector<std::complex<double>> generate_symmetry_matrix(int n) {
    std::vector<std::complex<double>> matrix(n*n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            if(i==j) matrix[i * n + j] = std::complex<double>(rand01, rand01);
            else{
                matrix[i * n + j] = std::complex<double>(rand01, rand01);
                matrix[j * n + i] = matrix[i * n + j];
            }
        }
    }
    return matrix;
}
void check_eq(std::complex<double>* data1, std::complex<double>* data2, int size, double eps)
    {
        for (int i = 0;i < size;++i)
        {
            EXPECT_NEAR(data1[i].real(), data2[i].real(), eps);
            EXPECT_NEAR(data1[i].imag(), data2[i].imag(), eps);
        }
    };

TEST(BSETest, skewSolver) {
    int my_rank, num_procs;
    Cblacs_pinfo(&my_rank, &num_procs); 

    int nA = 2;
    std::vector<double> ev(2*nA);
    Parallel_2D pM, pA, pA_glb;
    pM.init(2*nA, 2*nA, 1, MPI_COMM_WORLD, false);
    pA.set(nA, nA, 1, pM.blacs_ctxt);
    pA_glb.set(nA, nA, nA, pM.blacs_ctxt);
    std::vector<std::complex<double>> v(pM.get_local_size(), 0.0);
    std::vector<std::complex<double>> A_part(pA.get_local_size(), 0.0);
    std::vector<std::complex<double>> B_part(pA.get_local_size(), 0.0);
    std::vector<std::complex<double>> A_glb = {
    {3.0, 0.0}, {0.5, -1.0}, {0.5, 1.0}, {6.0, 0.0}
    };
    std::vector<std::complex<double>> B_glb = {
    {1.2, 0.6}, {0.4, 0.5}, {0.4, 0.5}, {1.4, 0.3}
    };
    Cpxgemr2d(nA, nA, A_glb.data(), 1, 1, pA_glb.desc,
        A_part.data(), 1, 1, pA.desc,
        pA_glb.blacs_ctxt);
    Cpxgemr2d(nA, nA, B_glb.data(), 1, 1, pA_glb.desc,
        B_part.data(), 1, 1, pA.desc,
        pA_glb.blacs_ctxt);
    BSE::solve_full(my_rank, A_part, B_part, pA, pM, ev, v);
    EXPECT_NEAR(ev[0], -6.127295611, 1e-8); 
    EXPECT_NEAR(ev[1], -2.299184312, 1e-8);
}

TEST(BSETest, TDASolver) {
    int my_rank, num_procs;
    Cblacs_pinfo(&my_rank, &num_procs); 
    int nA = 5;
    Parallel_2D pA_glb, pA;
    pA.init(nA, nA, 1, MPI_COMM_WORLD, false);
    pA_glb.set(nA, nA, nA, pA.blacs_ctxt);

    std::vector<double> ev(nA);
    std::vector<std::complex<double>> A_glb;
    if (my_rank == 0){
        A_glb = generate_conjugate_matrix(nA);
    }
    std::vector<std::complex<double>> A_part(pA.get_local_size(), 0.0);
    std::vector<std::complex<double>> v(pA.get_local_size(), 0.0);

    Cpxgemr2d(nA, nA, A_glb.data(), 1, 1, pA_glb.desc,
                A_part.data(), 1, 1, pA.desc,
                pA_glb.blacs_ctxt);

    std::vector<std::complex<double>> Hv(pA.get_local_size(), 0.0);
    std::vector<std::complex<double>> Ωv(pA.get_local_size(), 0.0);
    
    auto start = std::chrono::high_resolution_clock::now();
    BSE::solve_tda(my_rank, A_part, pA, ev, v);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0){
        std::cout << "BSE::solve_tda execution time: " << elapsed.count() << " seconds" << std::endl;    
    }

    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < pA.get_col_size(); ++i) {
        int col_glb = pA.local2global_col(i);
        for (int j = 0; j < pA.get_row_size(); ++j) {
            Ωv[i * pA.get_row_size() + j] = ev[col_glb] * v[i * pA.get_row_size() + j];
        }
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0)
        std::cout << "Ωv construction time: " << elapsed.count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    ScalapackConnector::gemm('N', 'N', nA, nA, nA, 1.0,
        A_part.data(), 1, 1, pA.desc,
        v.data(), 1, 1, pA.desc,
        0.0,
        Hv.data(), 1, 1, pA.desc);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0)
        std::cout << "Hv execution time: " << elapsed.count() << " seconds" << std::endl;
    check_eq(Hv.data(), Ωv.data(), pA.get_local_size(), 1e-6);
}

TEST(BSETest, skewSolver2) {
    int my_rank, num_procs;
    Cblacs_pinfo(&my_rank, &num_procs); 

    int nA = 5;
    int nM = 2 * nA; // Full matrix dimension
    Parallel_2D pM, pA, pA_glb;
    pM.init(nM, nM, 1, MPI_COMM_WORLD, false);
    pA.set(nA, nA, 1, pM.blacs_ctxt);
    pA_glb.set(nA, nA, nA, pM.blacs_ctxt);

    std::vector<std::complex<double>> A_glb, B_glb;
    if (my_rank==0){
        A_glb = generate_conjugate_matrix(nA);
        B_glb = generate_symmetry_matrix(nA);
    }
    std::vector<std::complex<double>> A_part(pA.get_local_size(), 0.0);
    std::vector<std::complex<double>> B_part(pA.get_local_size(), 0.0);
    std::vector<double> ev(nM, 0.0);
    std::vector<std::complex<double>> v(pM.get_local_size(), 0.0);
    Cpxgemr2d(nA, nA, A_glb.data(), 1, 1, pA_glb.desc,
                A_part.data(), 1, 1, pA.desc,
                pA_glb.blacs_ctxt);
    Cpxgemr2d(nA, nA, B_glb.data(), 1, 1, pA_glb.desc,
                B_part.data(), 1, 1, pA.desc,
                pA_glb.blacs_ctxt);
    // solve_full destroys A_part and B_part; pass copies to preserve originals
    auto A_copy = A_part;
    auto B_copy = B_part;
    auto start = std::chrono::high_resolution_clock::now();
    BSE::solve_full(my_rank, A_copy, B_copy, pA, pM, ev, v);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0){
        std::cout << "BSE::solve_full execution time: " << elapsed.count() << " seconds" << std::endl;    
    }

    std::vector<std::complex<double>> full_H (pM.get_local_size(), 0.0);
    std::vector<std::complex<double>> Hv(pM.get_local_size(), 0.0);
    std::vector<std::complex<double>> Ωv(pM.get_local_size(), 0.0);

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < pM.get_col_size(); ++i) {
        int col_glb = pM.local2global_col(i);
        for (int j = 0; j < pM.get_row_size(); ++j) {
            Ωv[i * pM.get_row_size() + j] = ev[col_glb] * v[i * pM.get_row_size() + j];
        }
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0)
        std::cout << "Ωv construction time: " << elapsed.count() << " seconds" << std::endl;

    BSE::arrayFlatten1(A_part, B_part, full_H, pA, pM);

    start = std::chrono::high_resolution_clock::now();
    ScalapackConnector::gemm('N', 'N', nM, nM, nM, 1.0,
        full_H.data(), 1, 1, pM.desc,
        v.data(), 1, 1, pM.desc,
        0.0,
        Hv.data(), 1, 1, pM.desc);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0)
        std::cout << "Hv execution time: " << elapsed.count() << " seconds" << std::endl;

    check_eq(Hv.data(), Ωv.data(), pM.get_local_size(), 1e-6);
    
    std::vector<std::complex<double>> identity(pM.get_local_size(), 0.0);
    for (int i = 0; i < pM.get_col_size(); ++i) {
        int col_glb = pM.local2global_col(i);
        int row_loc = pM.global2local_row(col_glb);
        if (row_loc != -1)
            identity[i * pM.get_row_size() + row_loc] = 1.0;
    }
    // v = ｢Y* X      left_v = ｢-Y*  X
    //      X* Y｣                X* -Y｣
    std::vector<std::complex<double>> left_v(pM.get_local_size(), 0.0);
    for (int i = 0; i < pM.get_col_size(); ++i) {
        int col_glb = pM.local2global_col(i);
        for (int j = 0; j < pM.get_row_size(); ++j) {
            int row_glb = pM.local2global_row(j);
            if (col_glb < nA && row_glb < nA)
                left_v[i * pM.get_row_size() + j] = -v[i * pM.get_row_size() + j];
            else if (col_glb >= nA && row_glb < nA)
                left_v[i * pM.get_row_size() + j] = v[i * pM.get_row_size() + j];
            else if (col_glb < nA && row_glb >= nA)
                left_v[i * pM.get_row_size() + j] = v[i * pM.get_row_size() + j];
            else // if (col_glb >= nA && row_glb >= nA)
                left_v[i * pM.get_row_size() + j] = -v[i * pM.get_row_size() + j];
        }
    }

    start = std::chrono::high_resolution_clock::now();
    ScalapackConnector::gemm('C', 'N', nM, nM, nM, 1.0,
        left_v.data(), 1, 1, pM.desc,
        v.data(), 1, 1, pM.desc,
        0.0,
        Ωv.data(), 1, 1, pM.desc);// overwrite Ωv by left_v.v
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    if (my_rank == 0)
        std::cout << "left_v.v execution time: " << elapsed.count() << " seconds" << std::endl;
    check_eq(Ωv.data(), identity.data(), pM.get_local_size(), 1e-6);
}


int main(int argc, char **argv) {
    srand(time(nullptr));
    int thread_level = -1;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_level);
    if (thread_level != MPI_THREAD_MULTIPLE)
    {
        std::cerr << "MPI_Init_thread request " << MPI_THREAD_MULTIPLE << " but provide " << thread_level << std::endl;
    }
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
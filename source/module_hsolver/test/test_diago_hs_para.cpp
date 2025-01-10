#ifndef TEST_DIAGO_PXXXGVX_H
#define TEST_DIAGO_PXXXGVX_H

#include <gtest/gtest.h>
#include <vector>
#include <complex>
#include <random>
#include <mpi.h>
#include <chrono>
#include <fstream>

#include "../diag_hs_para.h"
#include "module_hsolver/kernels/dngvd_op.h"

template <typename T>
typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
generate_random_hs_impl(int d, std::mt19937& gen, std::uniform_real_distribution<typename GetTypeReal<T>::type>& dis, std::vector<T>& h_mat, std::vector<T>& s_mat) {
    // For S matrix, firstly we generate a random symmetric matrix s_tmp, then we set S = s_tmp * s_tmp^T + n * I
    std::vector<T> s_tmp(d*d,0);
    for (int i = 0; i < d; ++i) {
        for (int j = i; j < d; ++j) {
            typename GetTypeReal<T>::type value1 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            h_mat[i * d + j] = value1;
            h_mat[j * d + i] = value1;

            // construct a random overlap matrix
            typename GetTypeReal<T>::type value2 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            s_tmp[i * d + j] = value2;
        }
    }

    // set S = s_tmp * s_tmp^T + n * I
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            s_mat[i * d + j] = 0;
            for (int k = 0; k < d; ++k) {
                s_mat[i * d + j] += s_tmp[i * d + k] * s_tmp[j * d + k];
            }
            if (i == j) {
                s_mat[i * d + j] += 2.0;
            }
        }
    }
}

template <typename T>
typename std::enable_if<std::is_same<T, std::complex<double>>::value || std::is_same<T, std::complex<float>>::value>::type
generate_random_hs_impl(int d, std::mt19937& gen, std::uniform_real_distribution<typename GetTypeReal<T>::type>& dis, std::vector<T>& h_mat, std::vector<T>& s_mat) {
    std::vector<T> s_tmp(d*d,0);
    for (int i = 0; i < d; ++i) {
        for (int j = i; j < d; ++j) {
            typename GetTypeReal<T>::type value1 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            typename GetTypeReal<T>::type value2 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            h_mat[i * d + j] = T(value1, value2);
            if (i != j)
            {
                h_mat[j * d + i] = T(value1, -value2);
            }
            else{
                h_mat[j * d + i] = T(value1, 0);
            }

            // construct a random overlap matrix
            value1 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            value2 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            s_tmp[i * d + j] = T(value1, value2);
        }
    }

    // set S = s_tmp * s_tmp^T + n * I
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            s_mat[i * d + j] = T(0, 0);
            for (int k = 0; k < d; ++k) {
                s_mat[i * d + j] += s_tmp[i * d + k] * std::conj(s_tmp[j * d + k]);
            }
            if (i == j) {
                s_mat[i * d + j] += T(2.0, 0);
            }
        }
    }
}

template <typename T>
void generate_random_hs(int d, int random_seed ,std::vector<T>& h_mat, std::vector<T>& s_mat) {
    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<typename GetTypeReal<T>::type> dis(-1.0,1.0);

    h_mat.resize(d * d);
    s_mat.resize(d * d);
    generate_random_hs_impl(d, gen, dis, h_mat, s_mat);
}


template <typename T>
typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
verify_results(const std::vector<T>& h_psi, const std::vector<T>& s_psi, const std::vector<typename GetTypeReal<T>::type>& ekb, int lda, int nbands, double threshold) {
    for (int i = 0; i < lda; ++i) {
        for (int j = 0; j < nbands; ++j) {
            ASSERT_NEAR(h_psi[j * lda + i], ekb[j] * s_psi[j * lda + i], threshold);
        }
    }
}

template <typename T>
typename std::enable_if<std::is_same<T, std::complex<double>>::value || std::is_same<T, std::complex<float>>::value>::type
verify_results(const std::vector<T>& h_psi, const std::vector<T>& s_psi, const std::vector<typename GetTypeReal<T>::type>& ekb, int lda, int nbands, double threshold) {
    for (int i = 0; i < lda; ++i) {
        for (int j = 0; j < nbands; ++j) {
            ASSERT_NEAR(h_psi[j * lda + i].real(), ekb[j] * s_psi[j * lda + i].real(), threshold);
            ASSERT_NEAR(h_psi[j * lda + i].imag(), ekb[j] * s_psi[j * lda + i].imag(), threshold);
        }
    }
}

template <typename T>
void test_diago_hs(int lda, int nb, int random_seed, int nbands, int diag_type, MPI_Comm comm) {
    // diag_type should be 1 (for elpa) or 2 (for scalapack)
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    std::vector<T> h_mat, s_mat, wfc, h_psi, s_psi;
    std::vector<typename GetTypeReal<T>::type> ekb(lda);
    if (my_rank==0)
    {
        h_mat.resize(lda * lda);
        s_mat.resize(lda * lda);
        wfc.resize(lda * lda);
        generate_random_hs(lda, random_seed, h_mat, s_mat);
    }
    hsolver::diago_hs_para<T>(h_mat.data(), s_mat.data(), lda, nbands,ekb.data(), wfc.data(), comm, diag_type, nb);

    // Verify results
    if (my_rank == 0){
        double threshold = 1e-6;
        if (std::is_same<T, std::complex<double>>::value || std::is_same<T, double>::value) {
            threshold = 1e-12;
        }

        h_psi.resize(lda * nbands, 0);
        s_psi.resize(lda * nbands, 0);

        for (int i = 0; i < lda; ++i) {
            for (int j = 0; j < nbands; ++j) {
                for (int k = 0; k < lda; ++k) {
                    h_psi[j * lda + i] += h_mat[k * lda + i] * wfc[j * lda + k];
                    s_psi[j * lda + i] += s_mat[k * lda + i] * wfc[j * lda + k];
                }
            }
        }
        verify_results<T>(h_psi, s_psi, ekb, lda, nbands, threshold);
    }
}

template <typename T>
void test_performance(int lda, int nb, int nbands, MPI_Comm comm,int case_numb, int loop_numb) {
    // generate 10 random H/S, and do the diagonalization 100 times by using elpa/scalapack and lapack.
    int my_rank, nproc;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nproc);

    std::vector<T> h_mat, s_mat, wfc, h_psi, s_psi;
    std::vector<typename GetTypeReal<T>::type> ekb_elpa(lda);
    std::vector<typename GetTypeReal<T>::type> ekb_scalap(lda);
    std::vector<typename GetTypeReal<T>::type> ekb_lapack(lda);

    if (my_rank==0)
    {
        std::cout << "\nMatrix size: " << lda << " x " << lda << std::endl;
        std::cout << "Number of bands: " << nbands << std::endl;
        std::cout << "Number of processors: " << nproc << std::endl;
        std::cout << "Block size of 2D distribution: " << nb << std::endl;
        h_mat.resize(lda * lda);
        s_mat.resize(lda * lda);
        wfc.resize(lda * lda);
    }

    // store all the times in a vector
    std::vector<double> time_elpa(case_numb, 0);
    std::vector<double> time_scalap(case_numb, 0);
    std::vector<double> time_lapack(case_numb, 0);

    if (my_rank == 0) { std::cout << "Random matrix ";
}
    for (int randomi = 0; randomi < case_numb; ++randomi) 
    {
        
        if (my_rank == 0) {
            std::cout << randomi << " ";
            generate_random_hs(lda, randomi, h_mat, s_mat);
        }

        // ELPA
        MPI_Barrier(comm);
        auto start = std::chrono::high_resolution_clock::now();
        for (int j=0;j<loop_numb;j++)
        {
            hsolver::diago_hs_para<T>(h_mat.data(), s_mat.data(), lda, nbands,ekb_elpa.data(), wfc.data(), comm, 1, nb);
            MPI_Barrier(comm);
        }
        MPI_Barrier(comm);
        auto end = std::chrono::high_resolution_clock::now();
        time_elpa[randomi] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


        // scalapack
        start = std::chrono::high_resolution_clock::now();
        for (int j=0;j<loop_numb;j++)
        {
            hsolver::diago_hs_para<T>(h_mat.data(), s_mat.data(), lda, nbands,ekb_scalap.data(), wfc.data(), comm, 2, nb);
            MPI_Barrier(comm);
        }
        MPI_Barrier(comm);
        end = std::chrono::high_resolution_clock::now();
        time_scalap[randomi] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //LApack
        if (my_rank == 0) 
        { 
            std::vector<T> h_tmp, s_tmp;
            start = std::chrono::high_resolution_clock::now();
            base_device::DEVICE_CPU* ctx = {};

            for (int j=0;j<loop_numb;j++)
            {
                h_tmp = h_mat;
                s_tmp = s_mat;
                hsolver::dngvx_op<T,base_device::DEVICE_CPU>()(ctx,
                                      lda,
                                      lda,
                                      h_tmp.data(),
                                      s_tmp.data(),
                                      nbands,
                                      ekb_lapack.data(),
                                      wfc.data());
            }
            end = std::chrono::high_resolution_clock::now();
            time_lapack[randomi] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            //COMPARE EKB
            for (int i = 0; i < nbands; ++i) {
                typename GetTypeReal<T>::type diff_elpa_lapack = std::abs(ekb_elpa[i] - ekb_lapack[i]);
                typename GetTypeReal<T>::type diff_scalap_lapack = std::abs(ekb_scalap[i] - ekb_lapack[i]);
                if (diff_elpa_lapack > 1e-6 || diff_scalap_lapack > 1e-6)
                {
                    std::cout << "eigenvalue " << i << " by ELPA: " << ekb_elpa[i] << std::endl;
                    std::cout << "eigenvalue " << i << " by Scalapack: " << ekb_scalap[i] << std::endl;
                    std::cout << "eigenvalue " << i << " by Lapack: " << ekb_lapack[i] << std::endl;
                }
            }
        }
        MPI_Barrier(comm);   

    }

    if (my_rank == 0)
    {
        std::cout << "\nELPA Time     : ";
        for (int i=0; i < case_numb;i++)
        {std::cout << time_elpa[i] << " ";}
        std::cout << std::endl;

        std::cout << "scalapack Time: ";
        for (int i=0; i < case_numb;i++)
        {std::cout << time_scalap[i] << " ";}
        std::cout << std::endl;

        std::cout << "lapack Time   : ";
        for (int i=0; i < case_numb;i++)
        {std::cout << time_lapack[i] << " ";}
        std::cout << std::endl;

        // print out the average time and speedup
        double avg_time_elpa = 0;
        double avg_time_scalap = 0;
        double avg_time_lapack = 0;
        for (int i=0; i < case_numb;i++)
        {
            avg_time_elpa += time_elpa[i];
            avg_time_scalap += time_scalap[i];
            avg_time_lapack += time_lapack[i];
        }

        avg_time_elpa /= case_numb;
        avg_time_scalap /= case_numb;
        avg_time_lapack /= case_numb;
        std::cout << "Average Lapack Time   : " << avg_time_lapack << " ms" << std::endl;
        std::cout << "Average ELPA Time     : " << avg_time_elpa << " ms, Speedup: " << avg_time_lapack / avg_time_elpa << std::endl;
        std::cout << "Average Scalapack Time: " << avg_time_scalap << " ms, Speedup: " << avg_time_lapack / avg_time_scalap << std::endl;
    }
}

//test_diago_hs(int lda, int nb, int random_seed, int nbands, int diag_type, MPI_Comm comm)
TEST(DiagoPxxxgvxElpaTest, Double) {
    test_diago_hs<double>(16, 4, 0, 10, 1, MPI_COMM_WORLD);
    test_diago_hs<double>(20, 6, 0, 18, 1, MPI_COMM_WORLD);
}

TEST(DiagoPxxxgvxElpaTest, ComplexDouble) {
    test_diago_hs<std::complex<double>>(16, 6, 0, 10, 1, MPI_COMM_WORLD);
    test_diago_hs<std::complex<double>>(20, 8, 0, 18, 1, MPI_COMM_WORLD);
}

TEST(DiagoPxxxgvxScalapackTest, Double) {
    test_diago_hs<double>(16, 4, 0, 10, 2,MPI_COMM_WORLD);
    test_diago_hs<double>(20, 6, 0, 18, 2,MPI_COMM_WORLD);
}

TEST(DiagoPxxxgvxScalapackTest, ComplexDouble) {
    test_diago_hs<std::complex<double>>(16, 4, 0, 10, 2, MPI_COMM_WORLD);
}
TEST(DiagoPxxxgvxScalapackTest, Float) {
    test_diago_hs<float>(16, 4, 0, 10,2,MPI_COMM_WORLD);
}

TEST(DiagoPxxxgvxScalapackTest, ComplexFloat) {
    test_diago_hs<std::complex<float>>(16, 4, 0, 10,2,MPI_COMM_WORLD);
}

//TEST(DiagoPxxxgvxPerformanceTest, Double) {
//    int ndim = 200;
//    int nband = 180;
//    int case_numb = 10;
//    int loop_numb = 10;
//    test_performance<std::complex<double>>(ndim, 32,  nband, MPI_COMM_WORLD, case_numb, loop_numb);
//}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    if (myrank != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        MPI_Finalize();
        return result;
	}

    MPI_Finalize();

	return 0;
}

#endif // TEST_DIAGO_PXXXGVX_H
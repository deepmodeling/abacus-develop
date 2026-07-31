#include "hamilt_bse_solver.h"

namespace BSE
{
template <typename T>
void printM(const std::vector<T>& A, int m, int n, std::string file, std::string name)
{
    std::ofstream ofs(file);
    ofs << name << std::endl;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ofs << std::setw(25) << A[j * m + i];
        }
        ofs << std::endl;
    }
    ofs.close();
}

void arrayFlatten1(const std::vector<std::complex<double>>& A,
                   const std::vector<std::complex<double>>& B,
                   std::vector<std::complex<double>>& result, // result={{A,B},{-B*,-A*}}
                   const Parallel_2D& pA,
                   const Parallel_2D& pM)
{
    int nA = pA.get_global_row_size();
#ifdef __MPI
    Cpxgemr2d(nA, nA, const_cast<std::complex<double>*>(A.data()), 1, 1, const_cast<int*>(pA.desc),
                result.data(), 1, 1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
    Cpxgemr2d(nA, nA, const_cast<std::complex<double>*>(B.data()), 1, 1, const_cast<int*>(pA.desc),
                result.data(), 1, nA+1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
    std::vector<std::complex<double>> temp_A(pA.get_local_size(), 0.0);
    std::vector<std::complex<double>> temp_B(pA.get_local_size(), 0.0);
    #pragma omp parallel for
    for (int i = 0; i < pA.get_local_size(); i++ ){
        temp_A[i] = -std::conj(A[i]);
        temp_B[i] = -std::conj(B[i]);
    }
    Cpxgemr2d(nA, nA, temp_A.data(), 1, 1, const_cast<int*>(pA.desc),
                result.data(), nA+1, nA+1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
    Cpxgemr2d(nA, nA, temp_B.data(), 1, 1, const_cast<int*>(pA.desc),
                result.data(), nA+1, 1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
#else
    for (int j = 0; j < nA; j++ ){
        for (int i = 0; i < nA; i++ ){
            result[ i + j*2*nA ] = A[i+j*nA];
            result[ nA+i + j*2*nA ] = -std::conj(B[i+j*nA]);
            result[ i + (nA+j)*2*nA ] = B[i+j*nA];
            result[ nA+i + (nA+j)*2*nA ] = -std::conj(A[i+j*nA]);
        }
    }
#endif
}

void arrayFlatten2(const std::vector<std::complex<double>>& A,
                   const std::vector<std::complex<double>>& B,
                   std::vector<double>& result, // result={{Re(A+B), Im(A-B)},{-Im(A+B), Re(A-B)}}
                   const Parallel_2D& pA,
                   const Parallel_2D& pM)
{
    int nA = pA.get_global_row_size();
#ifdef __MPI
    std::vector<double> temp(pA.get_local_size(), 0.0);
    for (int i = 0; i < pA.get_local_size(); i++ ){
        temp[i] = A[i].real() + B[i].real();
    }
    Cpxgemr2d(nA, nA, temp.data(), 1, 1, const_cast<int*>(pA.desc),
                result.data(), 1, 1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
    for (int i = 0; i < pA.get_local_size(); i++ ){
        temp[i] = A[i].imag() - B[i].imag();
    }
    Cpxgemr2d(nA, nA, temp.data(), 1, 1, const_cast<int*>(pA.desc),
                result.data(), 1, nA+1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
    for (int i = 0; i < pA.get_local_size(); i++ ){
        temp[i] = - A[i].imag() - B[i].imag();
    }
    Cpxgemr2d(nA, nA, temp.data(), 1, 1, const_cast<int*>(pA.desc),
                result.data(), nA+1, 1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
    for (int i = 0; i < pA.get_local_size(); i++ ){
        temp[i] = A[i].real() - B[i].real();
    }
    Cpxgemr2d(nA, nA, temp.data(), 1, 1, const_cast<int*>(pA.desc),
                result.data(), nA+1, nA+1, const_cast<int*>(pM.desc),
                pA.blacs_ctxt);
#else
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0; j < nA; j++ ){
        for (int i = 0; i < nA; i++ ){
            result[ i + j*2*nA ] = A[i+j*nA].real() + B[i+j*nA].real();
            result[ nA+i + j*2*nA ] = - A[i+j*nA].imag() - B[i+j*nA].imag();
            result[ i + (nA+j)*2*nA ] = A[i+j*nA].imag() - B[i+j*nA].imag();
            result[ nA+i + (nA+j)*2*nA ] = A[i+j*nA].real() - B[i+j*nA].real();
        }
    }
#endif
}

void solve_full(const int my_rank,
                std::vector<std::complex<double>>& A_part,
                std::vector<std::complex<double>>& B_part,
                const Parallel_2D& pA,
                const Parallel_2D& pM,
                std::vector<double>& ev,
                std::vector<std::complex<double>>& v)
{
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full");
    ModuleBase::timer::start("HamiltBSE", "elpa_solve_full");

    assert(pM.get_global_row_size() == pM.get_global_col_size());
    int n = pM.get_global_row_size(); // full matrix size
    int nA = pA.get_global_row_size();
    assert(n == 2*nA);

    elpa_t elpaInstance;
    int status;

    if (elpa_init(20210430) != ELPA_OK)
    {
        fprintf(stderr, "Error: ELPA API version not supported");
        exit(1);
    }

    elpaInstance = elpa_allocate(&status);
    if (status != ELPA_OK)
    {
        std::cout << "Could not allocate elpa instance" << std::endl;
        exit(1);
    }
    elpa_set(elpaInstance, "na", n, &status);
    elpa_set(elpaInstance, "nev", n, &status); // have to solve all ev, see step6: v = SQLzΩ^(-1/2)
    elpa_set(elpaInstance, "local_nrows", pM.get_row_size(), &status);
    elpa_set(elpaInstance, "local_ncols", pM.get_col_size(), &status);
    elpa_set(elpaInstance, "nblk", pM.get_block_size(), &status);
    elpa_set(elpaInstance, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &status);
    elpa_set(elpaInstance, "process_row", pM.get_coord_row(), &status);
    elpa_set(elpaInstance, "process_col", pM.get_coord_col(), &status);
    elpa_set(elpaInstance, "solver", ELPA_SOLVER_2STAGE, &status);
#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif
    elpa_set(elpaInstance, "omp_threads", num_threads, &status);
    status = elpa_setup(elpaInstance);
    if (status != ELPA_OK)
    {
        fprintf(stderr, "Could not set up the ELPA object");
    }

    // step1: construct M
    // M = {{Re(A_part+B_part), Im(A_part-B_part)}, {-Im(A_part+B_part), Re(A_part-B_part)}}
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full1");
    std::vector<double> M(pM.get_local_size());
    arrayFlatten2(A_part, B_part, M, pA, pM);
    std::vector<std::complex<double>>().swap(A_part);
    std::vector<std::complex<double>>().swap(B_part);

    // step2: construct J = {{0, I}, {-I, 0}}
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full2");
    std::vector<double> J(pM.get_local_size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0;j < pM.get_col_size();++j)
    {
        int col_glb = pM.local2global_col(j);
        for (int i = 0;i < pM.get_row_size();++i)
        {
            int row_glb = pM.local2global_row(i);
            if (col_glb - row_glb == nA)
                J[j * pM.get_row_size() + i] = 1.0;
            else if (row_glb - col_glb == nA)
                J[j * pM.get_row_size() + i] = -1.0;
        }
    }

    // step3: Cholesky factorization M = U^T U
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full3");
    elpa_cholesky(elpaInstance, M.data(), &status); // output M is upper triangular matrix U, M=U^T U

    // std::vector<double> global_U(n*n, 0.0);
    // Cpxgemr2d(n, n, M.data(), 1, 1, const_cast<int*>(pM.desc),
    //             global_U.data(), 1, 1, pM_glb.desc,
    //             pM.blacs_ctxt);
    // if (my_rank == 0) {std::cout<<"U:"<<std::endl; printM( global_U, n, n, "U", "U" );}

    // step4: compute anti-symmetric matrix UJL = U J U^T
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full4");
    std::vector<double> UJ (pM.get_local_size(), 0.0);

    ScalapackConnector::gemm('N', 'N', n, n, n, 1.0,
        M.data(), 1, 1, pM.desc, 
        J.data(), 1, 1, pM.desc,
        0.0,
        UJ.data(), 1, 1, pM.desc);

    // std::vector<double> global_UJ (n*n);
    // Cpxgemr2d(n, n, UJ.data(), 1, 1, const_cast<int*>(pM.desc),
    //             global_UJ.data(), 1, 1, pM_glb.desc,
    //             pM.blacs_ctxt);    
    //if (my_rank == 0) {std::cout<<"UJ:"<<std::endl; printM( global_UJ, n, n, "M", "M" );}

    ScalapackConnector::gemm('N', 'T', n, n, n, 1.0,
        UJ.data(), 1, 1, pM.desc,
        M.data(), 1, 1, pM.desc,
        0.0,
        J.data(), 1, 1, pM.desc); // J is UJL here
    std::vector<double>().swap(UJ);

    // std::vector<double> global_UJL(n*n);
    // Cpxgemr2d(n, n, J.data(), 1, 1, const_cast<int*>(pM.desc),
    //             global_UJL.data(), 1, 1, pM_glb.desc,
    //             pM.blacs_ctxt);
    //if (my_rank == 0) {std::cout<<"UJL:"<<std::endl; printM( global_UJL, n, n, "UJL", "UJL" );}

    // step5: compute eigenvalues ev and eigenvectors z of UJL
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full5");
    std::vector<double> z(2 * pM.get_local_size()); // 2 for elpa_skew stores complex as 2 double

    elpa_skew_eigenvectors(elpaInstance, J.data()/*UJL*/, ev.data(), z.data(), &status);
    assert(status == ELPA_OK);
    std::vector<double>().swap(J);
    // std::vector<std::complex<double>> global_z(n * n, 0.0);
    // Cpxgemr2d(n, n, z.data(), 1, 1, const_cast<int*>(pM.desc),
    //             global_z.data(), 1, 1, pM_glb.desc,
    //             pM.blacs_ctxt);

    // if (my_rank == 0) {
        // std::cout << "Eigenvalues computed successfully." << std::endl;
        // std::cout << "Eigenvalues:" << std::endl;
        // for (int i = 0; i < n; ++i) {
        //     std::cout << ev[i] << std::endl;
        // }
        //std::cout << "Eigenvectors (right):" << std::endl;
        //printM(global_z, n, n);
    //}
    // step6: compute normalized eigenvectors v = SQLzΩ^(-1/2), where Ω = diag(ev)
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full6.1");
    std::vector<double> Lz_real(pM.get_local_size(), 0.0);
    std::vector<double> Lz_imag(pM.get_local_size(), 0.0);

// 6.1: zΩ^(-1/2). NOTE: both positive and negative eigenvalues are handled here
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0; j < pM.get_col_size(); ++j)
    {
        double ev_sqrt = std::sqrt(std::abs(ev[pM.local2global_col(j)]));
        for (int i = 0; i < pM.get_row_size(); ++i)
        {
            z[j * pM.get_row_size() + i] /= ev_sqrt;                       // real part
            z[j * pM.get_row_size() + i + pM.get_local_size()] /= ev_sqrt; // imaginary part
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // 6.2: LzΩ^(-1/2), combine real and imaginary part
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full6.2");
    ScalapackConnector::gemm('T', 'N', n, n, n, 1.0,
        M.data(), 1, 1, pM.desc,
        z.data(), 1, 1, pM.desc,
        0.0,
        Lz_real.data(), 1, 1, pM.desc);
    ScalapackConnector::gemm('T', 'N', n, n, n, 1.0,
        M.data(), 1, 1, pM.desc,
        z.data() + pM.get_local_size(), 1, 1, pM.desc,
        0.0,
        Lz_imag.data(), 1, 1, pM.desc);

    std::vector<std::complex<double>> Lz(pM.get_local_size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0; j < pM.get_col_size(); ++j)
    {
        for (int i = 0; i < pM.get_row_size(); ++i)
        {
            Lz[j * pM.get_row_size() + i]
                = std::complex<double>(Lz_real[j * pM.get_row_size() + i], Lz_imag[j * pM.get_row_size() + i]);
        }
    }
    std::vector<double>().swap(Lz_real);
    std::vector<double>().swap(Lz_imag);

    std::vector<std::complex<double>> SQ(pM.get_local_size(), 0.0);// SQ = {{I, -i I}, {-I, -i I}} / sqrt(2)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int j = 0; j < pM.get_col_size(); ++j) {
        int col_glb = pM.local2global_col(j);
        for (int i = 0; i < pM.get_row_size(); ++i) {
            int row_glb = pM.local2global_row(i);
            if ( (col_glb < nA) && (row_glb < nA) && (col_glb == row_glb) ) {
                SQ[j * pM.get_row_size() + i] = 1/std::sqrt(2);
            }
            else if ( (col_glb >= nA) && (row_glb >= nA) && (col_glb == row_glb) ) {
                SQ[j * pM.get_row_size() + i] = {0, -1/std::sqrt(2)};
            }
            else if (col_glb == row_glb - nA) {
                SQ[j * pM.get_row_size() + i] = -1/std::sqrt(2);
            }
            else if (col_glb - nA == row_glb) {
                SQ[j * pM.get_row_size() + i] = {0, -1/std::sqrt(2)};
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // 6.3: v=SQLzΩ^{-1/2}
    ModuleBase::TITLE("HamiltBSE", "elpa_solve_full6.3");
    ScalapackConnector::gemm('N', 'N', n, n, n, 1.0,
        SQ.data(), 1, 1, pM.desc,
        Lz.data(), 1, 1, pM.desc,
        0.0,
        v.data(), 1, 1, pM.desc);

    elpa_deallocate(elpaInstance, &status);
    elpa_uninit(&status);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "elpa_solve_full");
    ModuleBase::timer::end("HamiltBSE", "elpa_solve_full");
}
} // namespace BSE
#ifdef __PEXSI
#include <mpi.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "module_hsolver/diago_pexsi.h"
#include "module_hsolver/module_pexsi/pexsi_solver.h"
#include "module_hsolver/test/diago_elpa_utils.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/parallel_global.h"
#include "module_base/global_variable.h"

#define PASSTHRESHOLD 5e-4
#define DETAILINFO    false
#define PRINT_HS      false
#define REPEATRUN     1

template <typename T> class HamiltTEST : public hamilt::Hamilt<T>
{
    public:
    int desc[9];
    int nrow, ncol;
    std::vector<T> h_local;
    std::vector<T> s_local;

    void matrix(hamilt::MatrixBlock<T> &hk_in, hamilt::MatrixBlock<T> &sk_in) 
    {
        hk_in = hamilt::MatrixBlock<T>{this->h_local.data(),
                                    (size_t)this->nrow,
                                    (size_t)this->ncol,
                                    this->desc};
        sk_in = hamilt::MatrixBlock<T>{this->s_local.data(),
                                    (size_t)this->nrow,
                                    (size_t)this->ncol,
                                    this->desc};
    }

    void constructHamilt(const int iter, const hamilt::MatrixBlock<double> rho) {}
    void updateHk(const int ik) {}
};



template<class T> class PexsiPrepare 
{
  public:
    PexsiPrepare(int nlocal, int nbands, int nb2d, int sparsity, std::string hfname, std::string sfname, std::string dmname)
        : nlocal(nlocal), nbands(nbands), nb2d(nb2d), sparsity(sparsity), hfname(hfname), sfname(sfname), dmname(dmname)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    }

    int dsize, myrank;
    int nlocal, nbands, nb2d, sparsity;
    int nprows, npcols, myprow, mypcol;
    int nrow, ncol;
    double hsolver_time = 0.0;
    std::string sfname, hfname, dmname;
    HamiltTEST<T> hmtest;
    std::vector<T> h;
    std::vector<T> s;
    std::vector<T> h_local;
    std::vector<T> s_local;
    psi::Psi<T> psi;
    hsolver::DiagoPexsi<T>* dh = nullptr;
    Parallel_Orbitals po;
    std::vector<double> abc;
    int icontxt;

    double mu;

    // density matrix
    std::vector<T*> dm_local;
    std::vector<T*> edm_local;

    std::vector<T> dm;

    bool read_HS()
    {
        bool readhfile = false;
        bool readsfile = false;
        if (this->myrank == 0)
        {
            int hdim, sdim;
            readhfile = LCAO_DIAGO_TEST::read_hs<std::vector<T>>(hfname, this->h);
            readsfile = LCAO_DIAGO_TEST::read_hs<std::vector<T>>(sfname, this->s);
            hdim = sqrt(this->h.size());
            sdim = sqrt(this->s.size());
            if (hdim != sdim)
            {
                printf("Error: dimensions of H and S are not equal, %d, %d\n", hdim, sdim);
                readhfile = readsfile = false;
            }
            nlocal = hdim;
        }
        MPI_Bcast(&nlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readhfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readsfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        nbands = nlocal/2;
        if (readhfile && readsfile)
            return true;
        return false;
    }

    bool produce_HS()
    {
        bool ok = this->read_HS();
        return ok;
    }

    void pb2d()
    {
        LCAO_DIAGO_TEST::process_2d(nprows, npcols, myprow, mypcol, icontxt);

        hmtest.nrow = LCAO_DIAGO_TEST::na_rc(nlocal, nb2d, nprows, myprow); // the number of row of the new_matrix in each process
        hmtest.ncol = LCAO_DIAGO_TEST::na_rc(nlocal, nb2d, npcols, mypcol); // the number of column of the new_matrix in each process


        int ISRC = 0, info;
        descinit_(hmtest.desc, &nlocal, &nlocal, &nb2d, &nb2d, &ISRC, &ISRC, &icontxt, &(hmtest.nrow), &info);
        if (info != 0)
        {
            printf("Invalid blacs-distribution. Abort!\n");
            exit(1);
        }

        
        // set po variables
        po.ncol = hmtest.ncol;
        po.nrow = hmtest.nrow;
        po.nb = nb2d;
        po.blacs_ctxt = icontxt;
        po.comm_2D = MPI_COMM_WORLD;
        po.dim0 = nprows;
        po.dim1 = npcols;
        po.testpb = true;

        if (DETAILINFO && myrank == 0)
        {
            std::cout << "nrow: " << hmtest.nrow << ", ncol: " << hmtest.ncol << ", nb: " << nb2d << std::endl;
        }

        dh = new hsolver::DiagoPexsi<T>(&po);
    }

    void distribute_data()
    {

        int local_size = hmtest.nrow * hmtest.ncol;
        this->h_local.resize(local_size);
        this->s_local.resize(local_size);

        LCAO_DIAGO_TEST::distribute_data<T>(this->h.data(),this->h_local.data(),nlocal,nb2d,hmtest.nrow,hmtest.ncol,icontxt);
        LCAO_DIAGO_TEST::distribute_data<T>(this->s.data(),this->s_local.data(),nlocal,nb2d,hmtest.nrow,hmtest.ncol,icontxt);
    }

    void set_env()
    {
        GlobalV::NLOCAL = nlocal;
        GlobalV::NBANDS = nbands;
        GlobalV::DSIZE = dsize;
        GlobalV::NSPIN = 1;
        DIAG_WORLD = MPI_COMM_WORLD;
        GlobalV::NPROC = dsize;

        psi.fix_k(0);
    }

    void set_pexsi_vars()
    {
        pexsi::PEXSI_Solver::set_pexsi_vars();
        pexsi::PEXSI_Solver::pexsi_mu = mu;
    }

    void diago()
    {
        if (DETAILINFO && myrank == 0)
        {
            std::cout << "Start to solve the KS equation using PEXSI" << std::endl;
        }
        this->pb2d();
        if (DETAILINFO && myrank == 0)
        {
            std::cout << "Finish the 2D parallelization" << std::endl;
        }
        this->distribute_data();
        if (DETAILINFO && myrank == 0)
        {
            std::cout << "Finish the data distribution" << std::endl;
        }
        this->set_env();
        if (DETAILINFO && myrank == 0)
        {
            std::cout << "Finish the environment setting" << std::endl;
        }
        double starttime = 0.0, endtime = 0.0;
        this->set_pexsi_vars();
        if (DETAILINFO && myrank == 0)
        {
            std::cout << "Finish the PEXSI setting" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        starttime = MPI_Wtime();
        for(int i=0;i<REPEATRUN;i++)
        {
            hmtest.h_local = this->h_local;
            hmtest.s_local = this->s_local;
            dh->diag(&hmtest, psi, nullptr);

            // copy the density matrix to dm_local
            dm_local = dh->DM;
            edm_local = dh->EDM;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (DETAILINFO && myrank == 0)
        {
            std::cout << "Finish the KS equation solving" << std::endl;
        }
        endtime = MPI_Wtime();
        hsolver_time = (endtime - starttime)/REPEATRUN;
    }

    bool read_ref()
    {
        auto f_dm = std::ifstream(dmname);
        if (!f_dm.is_open())
        {
            std::cout << "Error: cannot open the reference file " << dmname << std::endl;
            return false;
        }
        int nread = 0;
        f_dm >> nread;
        if (nread != nlocal)
        {
            std::cout << "Error: the number of global orbitals in the reference file is not equal to the current calculation" << std::endl;
            return false;
        }
        f_dm >> nread;
        if (nread != nlocal)
        {
            std::cout << "Error: the number of global orbitals in the reference file is not equal to the current calculation" << std::endl;
            return false;
        }

        f_dm >> GlobalV::nelec >> mu;
        
        dm.resize(nread*nread);
        // T* edm = new T[nglobal*nglobal];
        for (int i = 0; i < nread; i++)
        {
            for (int j = 0; j < nread; j++)
            {
                f_dm >> dm[i*nread+j];
            }
        }
        return true;
    }


    bool compare_ref(std::stringstream &out_info)
    {
        double maxerror = 0.0;
        int iindex = 0;
        bool pass = true;

        auto ofs = std::ofstream("dm_local" + std::to_string(myprow) + std::to_string(mypcol) + ".dat");

        int SENDPROW = 0, SENDPCOL = 0, tag = 0;

        // do iteration for matrix, distribute old_matrix to each process, pass a block each time
        for (int row = 0; row < nlocal; row++)
        {
            int recv_prow = (row / nb2d) % nprows; // the row number of recive process
            int nm_row = ((row / nb2d) / nprows) * nb2d + row % nb2d; // row number of block in new_matrix
            for (int col = 0; col < nlocal; col += nb2d)
            {
                int recv_pcol = (col / nb2d) % npcols; // the column number of recive process
                int nm_col = ((col / nb2d) / npcols) * nb2d + col % nb2d;
                int pass_length = std::min(nlocal - col, nb2d); // nlocal may not be devided by nb2d;
                // at global: nlocal * row + col + i
                // at local: (nm_col + i) * hmtest.nrow + nm_row

                if (myprow == recv_prow && mypcol == recv_pcol)
                {
                    double diff = 0;
                    for (int i = 0; i < pass_length; i++)
                    {
                        diff = std::abs(dm_local[0][(nm_col + i) * hmtest.nrow + nm_row] - dm[nlocal * row + col + i]);
                        if (diff > maxerror)
                        {
                            maxerror = diff;
                        }
                        if (diff > PASSTHRESHOLD)
                        {
                            pass = false;
                        }
                    }
                }
            }
        }

        bool pass_all = true;
        double maxerror_all = 0.0;
        MPI_Allreduce(&pass, &pass_all, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        MPI_Allreduce(&maxerror, &maxerror_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (myrank == 0)
        {
            std::cout << "H/S matrix are read from " << hfname << ", " << sfname << std::endl;
            std::cout << "Density matrix are read from " << dmname << std::endl;
            std::cout << std::endl;
            out_info << "Maximum difference between ks_hsolver and ref is " << maxerror_all << ", the pass threshold is " << PASSTHRESHOLD << std::endl;

            if (DETAILINFO)
            {
                std::cout << out_info.str();
                out_info.str("");
                out_info.clear();
            }
        }
        delete dh;
        return pass_all;
    }
};

class PexsiGammaOnlyTest : public ::testing::TestWithParam<PexsiPrepare<double>> {};

TEST_P(PexsiGammaOnlyTest, LCAO)
{
    std::stringstream out_info;
    PexsiPrepare<double> dp = GetParam();
    if (DETAILINFO && dp.myrank == 0)
    {
        std::cout << "nlocal: " << dp.nlocal << ", nbands: " << dp.nbands << ", nb2d: " << dp.nb2d << ", sparsity: " << dp.sparsity << std::endl;
    
    }
    ASSERT_TRUE(dp.produce_HS());
    if (DETAILINFO && dp.myrank == 0)
    {
        std::cout << "H/S matrix are read from " << dp.hfname << ", " << dp.sfname << std::endl;
    }
    ASSERT_TRUE(dp.read_ref());
    if (DETAILINFO && dp.myrank == 0)
    {
        std::cout << "Density matrix are read from " << dp.dmname << std::endl;
    }
    dp.diago();
    if (DETAILINFO && dp.myrank == 0)
    {
        std::cout << "Time for hsolver: " << dp.hsolver_time << "s" << std::endl;
    }

    bool pass = dp.compare_ref(out_info);
    EXPECT_TRUE(pass) << out_info.str();
    
    MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(
    DiagoTest,
    PexsiGammaOnlyTest,
    ::testing::Values( //int nlocal, int nbands, int nb2d, int sparsity, std::string ks_solver_in, std::string hfname, std::string sfname
        PexsiPrepare<double>(0, 0, 2, 0, "PEXSI-H-GammaOnly-Si2.dat", "PEXSI-S-GammaOnly-Si2.dat", "PEXSI-DM-GammaOnly-Si2.dat"),
        PexsiPrepare<double>(0, 0, 1, 0, "PEXSI-H-GammaOnly-Si2.dat", "PEXSI-S-GammaOnly-Si2.dat", "PEXSI-DM-GammaOnly-Si2.dat")

    ));


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int mypnum, dsize;
    MPI_Comm_size(MPI_COMM_WORLD, &dsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);

    testing::InitGoogleTest(&argc, argv);
    //Parallel_Global::split_diag_world(dsize);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (mypnum != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }
    int result = RUN_ALL_TESTS();

    if (mypnum == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }
    else
    {
        MPI_Finalize();
        return 0;
    }
}


#endif // __PEXSI
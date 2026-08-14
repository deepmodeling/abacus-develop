#include <chrono>
#include <array>
#include <cmath>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_cell/klist.h"

namespace
{
// Fill a Hermitian matrix H (row-major, nw*nw) with deterministic pseudo-random
// values. If real_only is true, H is real symmetric (used for k-points that are
// their own time-reversal partner).
void fill_hermitian(std::vector<std::complex<double>>& H, const int nw, const unsigned seed, const bool real_only = false)
{
    H.assign(nw * nw, std::complex<double>(0.0, 0.0));
    for (int i = 0; i < nw; ++i)
    {
        for (int j = i; j < nw; ++j)
        {
            const double re = 0.5 * std::sin(double(i * 7 + j * 13 + seed * 3));
            const double im = (i == j || real_only) ? 0.0 : 0.5 * std::cos(double(i * 11 + j * 5 + seed * 7));
            H[i * nw + j] = std::complex<double>(re, im);
            H[j * nw + i] = std::complex<double>(re, -im);
        }
    }
}

// Build an HContainer with the given (iat1, iat2, R) entries.
// If paired is true, both (iat1,iat2,R) and (iat2,iat1,-R) are inserted.
template <typename T>
hamilt::HContainer<T> build_r_container(
    const Parallel_Orbitals* paraV,
    const std::vector<std::pair<int, int>>& pairs,
    const std::vector<ModuleBase::Vector3<int>>& r_list,
    const bool paired)
{
    hamilt::HContainer<T> hR(paraV);
    for (const auto& pr: pairs)
    {
        for (const auto& R: r_list)
        {
            hR.insert_pair(hamilt::AtomPair<T>(pr.first, pr.second, R, paraV));
            if (paired)
            {
                hR.insert_pair(hamilt::AtomPair<T>(pr.second, pr.first, -R, paraV));
            }
        }
    }
    hR.allocate(nullptr, true);
    return hR;
}
} // namespace

/************************************************
 *  unit test of DensityMatrix constructor
 ***********************************************/

/**
 * This unit test construct a DensityMatrix object
 */

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;

class DMTest : public testing::Test
{
  protected:
    Parallel_Orbitals* paraV;
    int dsize;
    int my_rank = 0;
    UnitCell ucell;
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        // set up a unitcell, with one element and test_size atoms, each atom has test_nw orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw);
        ucell.atoms[0].iw2m.resize(test_nw);
        ucell.atoms[0].iw2n.resize(test_nw);
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 0;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();

        // set paraV
        init_parav();
    }

    void TearDown() override
    {
        delete paraV;
        delete[] ucell.atoms;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 2;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif
};

TEST_F(DMTest, cal_DMR_full)
{
    // get my rank of this process
    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    // output dim and nrow, ncol
    if (my_rank == 0)
    {
        std::cout << "my rank: " << my_rank << " dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    else
    {
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    // initalize a kvectors, Gamma-only
    K_Vectors* kv = nullptr;
    int nspin = 4;
    int nks = 2; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    // construct DM
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks());
    // set this->_DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks(); ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, std::complex<double>(0.77, 0.77));
                }
            }
        }
    }
    // initialize dmR_full
    hamilt::HContainer<std::complex<double>> dmR_full(ucell, paraV);
    // calculate this->_DMR
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_DMR_full(&dmR_full);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result
    for (int i = 0; i < dmR_full.size_atom_pairs(); i++)
    {
        const std::complex<double>* ptr1 = dmR_full.get_atom_pair(i).get_HR_values(0, 0, 0).get_pointer();
        //
        for (int j = 0; j < dmR_full.get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j].real(), 1.54, 1e-10);
            EXPECT_NEAR(ptr1[j].imag(), 1.54, 1e-10);
        }
    }
    delete kv;
}

TEST_F(DMTest, cal_DMR_blas_double)
{
    // get my rank of this process
    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    // output dim and nrow, ncol
    if (my_rank == 0)
    {
        std::cout << "my rank: " << my_rank << " dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    else
    {
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    // initalize a kvectors, Gamma-only
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 2; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    // construct DM
    elecstate::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    // set this->_DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, 0.77);
                }
            }
        }
    }
    // initialize this->_DMR
    Grid_Driver gd(0, 0);
    DM.init_DMR(&gd, &ucell);
    // set Gamma-only
    for (int is = 1; is <= nspin; is++)
    {
        DM.get_DMR_pointer(is)->fix_gamma();
    }
    // calculate this->_DMR
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_DMR();
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(0, 0, 0).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], 0.77, 1e-10);
        }
    }
    delete kv;
}

TEST_F(DMTest, cal_DMR_blas_complex)
{
    // get my rank of this process
    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    // output dim and nrow, ncol
    if (my_rank == 0)
    {
        std::cout << "my rank: " << my_rank << " dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    else
    {
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct DM
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    // set this->_DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, is * 0.77 * (ik + 1));
                }
            }
        }
    }
    // initialize this->_DMR
    Grid_Driver gd(0, 0);
    DM.init_DMR(&gd, &ucell);
    // calculate this->_DMR
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_DMR();
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result for spin-up
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77, 1e-10);
        }
    }
    // compare the result for spin-down
    for (int i = 0; i < DM.get_DMR_pointer(2)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(2)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(2)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77 * 2, 1e-10);
        }
    }
    // calculate DMR_total
    DM.switch_dmr(1);
    // compare the result for spin-up after sum
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77 * 3, 1e-10);
        }
    }
    // restore to normal DMR 
    DM.switch_dmr(0);
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77, 1e-10);
        }
    }
    // calculate DMR_differenct
    DM.switch_dmr(2);
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], 0.77, 1e-10);
        }
    }
    delete kv;
}

// T1: Fourier round-trip consistency between cal_DMR (e^{-ikR}) and
// folding_HR (e^{+ikR}). This is the test that locks the Fourier sign of the
// inverse transform: with the k-phase flipped (density_matrix.cpp: -sinp ->
// +sinp), D~(k) = DMK(-k) = DMK(k)† and the assertions below turn red.
// See docs/dm_dmk_dmr_action_plan.md (T1).
TEST_F(DMTest, T1_fourier_round_trip)
{
    // k-grid {0, 1/4, 1/2, 3/4} along x: contains non-Gamma k-points and k and
    // -k as distinct points (a {0, 1/2} grid is blind to the sign because
    // -k == k there).
    const int nk = 4;
    const int nR = 4;
    std::vector<ModuleBase::Vector3<double>> kvec_d(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        kvec_d[ik] = ModuleBase::Vector3<double>(0.25 * ik, 0.0, 0.0);
    }

    // Time-reversal constrained Hermitian DMK per k-point: k=0 and k=1/2 are
    // their own TR partner (real symmetric), k=1/4 and k=3/4 are conjugate
    // pairs. This makes D(R) real, so the real DMR path (which stores
    // Re[D(R)]) is lossless and the full complex round trip is exact.
    const int nw = test_nw;
    std::vector<std::vector<std::complex<double>>> M(nk);
    fill_hermitian(M[0], nw, 1, true);
    fill_hermitian(M[2], nw, 2, true);
    fill_hermitian(M[1], nw, 3, false);
    M[3].assign(nw * nw, std::complex<double>(0.0, 0.0));
    for (int i = 0; i < nw * nw; ++i)
    {
        M[3][i] = std::conj(M[1][i]);
    }

    // ---- real DMR path (nspin<4): DensityMatrix<complex, double> ----
    {
        elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, 1, kvec_d, nk);
        // set_DMK stores the transposed element (see density_matrix.h): to set
        // the element read by cal_DMR as (mu,nu), pass (nu,mu) to the setter.
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int i = 0; i < nw; ++i)
            {
                for (int j = 0; j < nw; ++j)
                {
                    DM.set_DMK(1, ik, j, i, M[ik][i * nw + j]);
                }
            }
        }
        // synthetic HContainer: pair (0,0) with the complete set of R
        // representatives {0,1,2,3} of the k-grid period
        hamilt::HContainer<double> hR(paraV);
        for (int rx = 0; rx < nR; ++rx)
        {
            hR.insert_pair(hamilt::AtomPair<double>(0, 0, rx, 0, 0, paraV));
        }
        hR.allocate(nullptr, true);
        DM.init_DMR(hR);
        DM.cal_DMR();

        for (int ik = 0; ik < nk; ++ik)
        {
            std::vector<std::complex<double>> hk(paraV->nrow * paraV->ncol, std::complex<double>(0.0, 0.0));
            hamilt::folding_HR(*DM.get_DMR_pointer(1), hk.data(), kvec_d[ik], paraV->ncol, 0);
            for (int i = 0; i < nw; ++i)
            {
                for (int j = 0; j < nw; ++j)
                {
                    // folding_HR has no weight normalization: on a complete
                    // R-representative set D~(k) = nR * DMK(k) (in production
                    // the k weights are embedded in DMK).
                    const std::complex<double> expect = double(nR) * M[ik][i * nw + j];
                    EXPECT_NEAR(hk[i * paraV->ncol + j].real(), expect.real(), 1e-10)
                        << "real DMR path, ik=" << ik << " mu=" << i << " nu=" << j;
                    EXPECT_NEAR(hk[i * paraV->ncol + j].imag(), expect.imag(), 1e-10)
                        << "real DMR path, ik=" << ik << " mu=" << i << " nu=" << j;
                }
            }
        }
    }

    // ---- complex DMR path (nspin=4 via cal_DMR_full) ----
    {
        elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, 4, kvec_d, nk);
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int i = 0; i < nw; ++i)
            {
                for (int j = 0; j < nw; ++j)
                {
                    DM.set_DMK(1, ik, j, i, M[ik][i * nw + j]);
                }
            }
        }
        hamilt::HContainer<std::complex<double>> hR(paraV);
        for (int rx = 0; rx < nR; ++rx)
        {
            hR.insert_pair(hamilt::AtomPair<std::complex<double>>(0, 0, rx, 0, 0, paraV));
        }
        hR.allocate(nullptr, true);
        DM.cal_DMR_full(&hR);

        for (int ik = 0; ik < nk; ++ik)
        {
            std::vector<std::complex<double>> hk(paraV->nrow * paraV->ncol, std::complex<double>(0.0, 0.0));
            hamilt::folding_HR(hR, hk.data(), kvec_d[ik], paraV->ncol, 0);
            for (int i = 0; i < nw; ++i)
            {
                for (int j = 0; j < nw; ++j)
                {
                    const std::complex<double> expect = double(nR) * M[ik][i * nw + j];
                    EXPECT_NEAR(hk[i * paraV->ncol + j].real(), expect.real(), 1e-10)
                        << "complex DMR path, ik=" << ik << " mu=" << i << " nu=" << j;
                    EXPECT_NEAR(hk[i * paraV->ncol + j].imag(), expect.imag(), 1e-10)
                        << "complex DMR path, ik=" << ik << " mu=" << i << " nu=" << j;
                }
            }
        }
    }
}

// T2: DMR Hermiticity. For a Hermitian DMK and a paired HContainer
// ((iat1,iat2,R) and (iat2,iat1,-R)), cal_DMR must satisfy
//   D(iat2,iat1,-R) = D(iat1,iat2,R)^T   (real DMR, nspin<4)
//   D(iat2,iat1,-R) = D(iat1,iat2,R)^\dagger (complex DMR, nspin=4)
// See docs/dm_dmk_dmr_action_plan.md (T2).
TEST_F(DMTest, T2_dmr_hermiticity)
{
    // two k-points, both non-Gamma, k and -k distinct
    std::vector<ModuleBase::Vector3<double>> kvec_d = {
        ModuleBase::Vector3<double>(0.25, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.75, 0.0, 0.0)};
    const int nk = 2;
    const int nw = test_nw;

    // Hermitian DMK (arbitrary, no TRS constraint needed for Hermiticity)
    std::vector<std::vector<std::complex<double>>> M(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        fill_hermitian(M[ik], nw, 10 + ik);
    }

    // pairs: cross (0,1)/(1,0) and self (0,0); R: 3 non-zero + origin
    const std::vector<std::pair<int, int>> pairs = {{0, 1}, {0, 0}};
    const std::vector<ModuleBase::Vector3<int>> r_list = {
        ModuleBase::Vector3<int>(0, 0, 0),
        ModuleBase::Vector3<int>(1, 0, 0),
        ModuleBase::Vector3<int>(0, 1, 0),
        ModuleBase::Vector3<int>(0, 0, 1)};

    // ---- real DMR path ----
    {
        elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, 1, kvec_d, nk);
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int i = 0; i < nw; ++i)
            {
                for (int j = 0; j < nw; ++j)
                {
                    DM.set_DMK(1, ik, j, i, M[ik][i * nw + j]);
                }
            }
        }
        hamilt::HContainer<double> hR = build_r_container<double>(paraV, pairs, r_list, true);
        DM.init_DMR(hR);
        DM.cal_DMR();
        const hamilt::HContainer<double>* dmr = DM.get_DMR_pointer(1);
        for (int iap = 0; iap < dmr->size_atom_pairs(); ++iap)
        {
            const hamilt::AtomPair<double>& ap = dmr->get_atom_pair(iap);
            const int iat1 = ap.get_atom_i();
            const int iat2 = ap.get_atom_j();
            for (int ir = 0; ir < ap.get_R_size(); ++ir)
            {
                const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                const hamilt::BaseMatrix<double>* mat = ap.find_matrix(R);
                const hamilt::BaseMatrix<double>* mat_rev = dmr->find_matrix(iat2, iat1, -R.x, -R.y, -R.z);
                ASSERT_NE(mat_rev, nullptr) << "missing reverse pair (" << iat2 << "," << iat1 << "," << -R.x << "," << -R.y << "," << -R.z << ")";
                const int row_size = ap.get_row_size();
                const int col_size = ap.get_col_size();
                for (int i = 0; i < row_size; ++i)
                {
                    for (int j = 0; j < col_size; ++j)
                    {
                        EXPECT_NEAR(mat_rev->get_pointer()[j * col_size + i],
                                    mat->get_pointer()[i * col_size + j],
                                    1e-10)
                            << "real DMR hermiticity: iat1=" << iat1 << " iat2=" << iat2
                            << " R=(" << R.x << "," << R.y << "," << R.z << ")";
                    }
                }
            }
        }
    }

    // ---- complex DMR path (cal_DMR_full) ----
    {
        elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, 4, kvec_d, nk);
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int i = 0; i < nw; ++i)
            {
                for (int j = 0; j < nw; ++j)
                {
                    DM.set_DMK(1, ik, j, i, M[ik][i * nw + j]);
                }
            }
        }
        hamilt::HContainer<std::complex<double>> hR = build_r_container<std::complex<double>>(paraV, pairs, r_list, true);
        DM.cal_DMR_full(&hR);
        for (int iap = 0; iap < hR.size_atom_pairs(); ++iap)
        {
            const hamilt::AtomPair<std::complex<double>>& ap = hR.get_atom_pair(iap);
            const int iat1 = ap.get_atom_i();
            const int iat2 = ap.get_atom_j();
            for (int ir = 0; ir < ap.get_R_size(); ++ir)
            {
                const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                const hamilt::BaseMatrix<std::complex<double>>* mat = ap.find_matrix(R);
                const hamilt::BaseMatrix<std::complex<double>>* mat_rev = hR.find_matrix(iat2, iat1, -R.x, -R.y, -R.z);
                ASSERT_NE(mat_rev, nullptr) << "missing reverse pair (" << iat2 << "," << iat1 << "," << -R.x << "," << -R.y << "," << -R.z << ")";
                const int row_size = ap.get_row_size();
                const int col_size = ap.get_col_size();
                for (int i = 0; i < row_size; ++i)
                {
                    for (int j = 0; j < col_size; ++j)
                    {
                        EXPECT_NEAR(mat_rev->get_pointer()[j * col_size + i].real(),
                                    mat->get_pointer()[i * col_size + j].real(),
                                    1e-10)
                            << "complex DMR hermiticity real: iat1=" << iat1 << " iat2=" << iat2
                            << " R=(" << R.x << "," << R.y << "," << R.z << ")";
                        EXPECT_NEAR(mat_rev->get_pointer()[j * col_size + i].imag(),
                                    -mat->get_pointer()[i * col_size + j].imag(),
                                    1e-10)
                            << "complex DMR hermiticity imag: iat1=" << iat1 << " iat2=" << iat2
                            << " R=(" << R.x << "," << R.y << "," << R.z << ")";
                    }
                }
            }
        }
    }
}

// T8: full-direction pairing storage guard. Every (iat1,iat2,R) block of the
// DMR HContainer built by init_DMR must have a (iat2,iat1,-R) counterpart.
// gint_rho and the force/stress paths rely on this pairing (closed-trace
// protection, see docs/dm_dmk_dmr_action_plan.md T8). A future half-set
// storage optimization will turn this test red.
TEST_F(DMTest, T8_full_direction_pairing_guard)
{
    // multi-k DM (TK=complex) keeps the R structure of init_DMR(Record_adj);
    // the TK=double path collapses to gamma-only by design.
    const std::vector<ModuleBase::Vector3<double>> kvec = {ModuleBase::Vector3<double>(0.1, 0.0, 0.0)};
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, 1, kvec, 1);

    // synthetic full symmetric neighbor list (as produced by Record_adj::cal_adj)
    const int nat = test_size;
    Record_adj ra;
    ra.na_each = new int[nat];
    ra.info = new int**[nat];
    const std::vector<ModuleBase::Vector3<int>> r_list = {
        ModuleBase::Vector3<int>(0, 0, 0),
        ModuleBase::Vector3<int>(1, 0, 0),
        ModuleBase::Vector3<int>(0, 1, 0),
        ModuleBase::Vector3<int>(0, 0, 1)};
    std::vector<std::vector<std::array<int, 5>>> lists(nat);
    for (int iat1 = 0; iat1 < nat; ++iat1)
    {
        for (int iat2 = 0; iat2 < nat; ++iat2)
        {
            for (const auto& R: r_list)
            {
                // both directions: (iat1,iat2,R) and (iat2,iat1,-R)
                lists[iat1].push_back({R.x, R.y, R.z, 0, iat2});
                lists[iat2].push_back({-R.x, -R.y, -R.z, 0, iat1});
            }
        }
    }
    for (int iat = 0; iat < nat; ++iat)
    {
        ra.na_each[iat] = lists[iat].size();
        ra.info[iat] = new int*[ra.na_each[iat]];
        for (int ad = 0; ad < ra.na_each[iat]; ++ad)
        {
            ra.info[iat][ad] = new int[5];
            for (int k = 0; k < 5; ++k)
            {
                ra.info[iat][ad][k] = lists[iat][ad][k];
            }
        }
    }

    DM.init_DMR(ra, &ucell);

    const hamilt::HContainer<double>* dmr = DM.get_DMR_pointer(1);
    EXPECT_GT(dmr->size_atom_pairs(), 0);
    for (int iap = 0; iap < dmr->size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<double>& ap = dmr->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
            EXPECT_NE(dmr->find_matrix(iat2, iat1, -R.x, -R.y, -R.z), nullptr)
                << "pairing guard: missing (iat1,iat2,R)=(" << iat1 << "," << iat2 << ","
                << R.x << "," << R.y << "," << R.z << ") reverse";
        }
    }

    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ad = 0; ad < ra.na_each[iat]; ++ad)
        {
            delete[] ra.info[iat][ad];
        }
        delete[] ra.info[iat];
    }
    delete[] ra.info;
    delete[] ra.na_each;
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}

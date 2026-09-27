#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_cell/klist.h"

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

TEST_F(DMTest, cal_dmr_full)
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
    module_dm::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks());
    // set this->dmk
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks(); ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_dmk(is, ik, i, j, std::complex<double>(0.77, 0.77));
                }
            }
        }
    }
    // initialize dmR_full
    hamilt::HContainer<std::complex<double>> dmR_full(ucell, paraV);
    // calculate this->dmr
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_dmr_full(&dmR_full, -1);
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
    module_dm::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    // set this->dmk
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_dmk(is, ik, i, j, 0.77);
                }
            }
        }
    }
    // initialize this->dmr
    Grid_Driver gd(0, 0);
    DM.init_dmr(&gd, &ucell);
    // set Gamma-only
    for (int is = 1; is <= nspin; is++)
    {
        DM.get_dmr_ptr(is)->fix_gamma();
    }
    // calculate this->dmr
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_dmr(-1);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result
    for (int i = 0; i < DM.get_dmr_ptr(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_dmr_ptr(1)->get_atom_pair(i).get_HR_values(0, 0, 0).get_pointer();
        //
        for (int j = 0; j < DM.get_dmr_ptr(1)->get_atom_pair(i).get_size(); j++)
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
    module_dm::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    // set this->dmk
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_dmk(is, ik, i, j, is * 0.77 * (ik + 1));
                }
            }
        }
    }
    // initialize this->dmr
    Grid_Driver gd(0, 0);
    DM.init_dmr(&gd, &ucell);
    // calculate this->dmr
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_dmr(-1);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result for spin-up
    for (int i = 0; i < DM.get_dmr_ptr(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_dmr_ptr(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_dmr_ptr(1)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77, 1e-10);
        }
    }
    // compare the result for spin-down
    for (int i = 0; i < DM.get_dmr_ptr(2)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_dmr_ptr(2)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_dmr_ptr(2)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77 * 2, 1e-10);
        }
    }
    // calculate DMR_total
    DM.switch_dmr(1);
    // compare the result for spin-up after sum
    for (int i = 0; i < DM.get_dmr_ptr(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_dmr_ptr(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_dmr_ptr(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77 * 3, 1e-10);
        }
    }
    // restore to normal DMR 
    DM.switch_dmr(0);
    for (int i = 0; i < DM.get_dmr_ptr(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_dmr_ptr(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_dmr_ptr(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77, 1e-10);
        }
    }
    // calculate DMR_differenct
    DM.switch_dmr(2);
    for (int i = 0; i < DM.get_dmr_ptr(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_dmr_ptr(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_dmr_ptr(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], 0.77, 1e-10);
        }
    }
    delete kv;
}

// Regression test for the SOC/noncollinear (global nspin==4) cal_dmr path.
//
// Background: in a real SOC run allocate_dm.cpp constructs the DensityMatrix with
//   spin_mult = 1   (the 2x2 spin block is stored as ONE doubled matrix),
// while the GLOBAL physical nspin is 4. cal_dmr must still take the
// spin-resolved (Pauli) branch, which folds each 2x2 complex spin block into
// (rho_0, rho_x, rho_y, rho_z) via xyz_to_updown(). That branch used to be
// selected by the global PARAM.inp.nspin==4; a refactor (commit dcad8913d)
// switched the condition to dm.spin_mult==4, which is never true in SOC
// (spin_mult==1), silently dropping the rho_x/y/z spin channels and producing a
// wrong charge density (tests/03_NAO_multik/*spin4* failed by ~41 eV).
//
// This test reproduces the real SOC construction (spin_mult=1, nspin=4)
// and fills the spin-diagonal DMK entries (uu, dd) with (a + i b), leaving the
// spin off-diagonal entries (ud, du) zero. It then checks that cal_dmr selects
// the Pauli branch:
//   * correct (Pauli) branch : rho_0 = (uu+dd).real() = 2a, rho_z = (uu-dd).real() = 0
//   * wrong   (real-project) : every element = a  (imaginary part b dropped)
// With the pre-fix condition (spin_mult==4 never taken) this test FAILS because
// rho_0 would come out as a instead of 2a.
TEST_F(DMTest, cal_DMR_soc_pauli_branch)
{
    // SOC doubles the orbital dimension (npol = 2): each atom carries nw*npol rows/cols.
    // The fixture's ucell has nw = test_nw, so the doubled global dimension is used here.
    const int npol = 2;
    const int global_dim_soc = test_size * test_nw * npol;
#ifdef __MPI
    Parallel_Orbitals* pv_soc = new Parallel_Orbitals();
    pv_soc->init(global_dim_soc, global_dim_soc, 2, MPI_COMM_WORLD);
    // build the iat2iwt map for the doubled (spinor) orbital count
    std::vector<int> iat2iwt_soc(test_size);
    for (int iat = 0; iat < test_size; ++iat)
    {
        iat2iwt_soc[iat] = iat * test_nw * npol;
    }
    pv_soc->set_atomic_trace(iat2iwt_soc.data(), test_size, global_dim_soc);
#else
    Parallel_Orbitals* pv_soc = paraV; // fallback; MPI path is the supported configuration
#endif

    // a single Gamma k-point; construct exactly like the real SOC setup_dm does:
    // nspin_dm = 1, but the global physical nspin = 4.
    std::vector<ModuleBase::Vector3<double>> kvec_d(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    const int nspin_dm = 1;
    const int nspin_global = 4;
    module_dm::DensityMatrix<std::complex<double>, double> DM(pv_soc, nspin_dm, kvec_d, 1, nspin_global);

    // fill the single DMK: spin-diagonal entries (uu, dd) are (a + i b),
    // spin off-diagonal entries (ud, du) stay zero. With a constant fill the
    // 2x2 spin block would also have ud = du = (a + i b), and rho_x =
    // Re(ud + du) would correctly be 2a instead of the asserted 0.
    const double a = 0.5;
    const double b = 0.25;
    for (int i = 0; i < pv_soc->nrow; i++)
    {
        for (int j = 0; j < pv_soc->ncol; j++)
        {
            // global spinor indices determine the spin parity; local indices
            // need not preserve parity under a 2D block-cyclic distribution
            const bool same_spin = (pv_soc->local2global_row(i) % npol)
                                   == (pv_soc->local2global_col(j) % npol);
            const std::complex<double> dmk_value = same_spin
                                                       ? std::complex<double>(a, b)
                                                       : std::complex<double>(0.0, 0.0);
            DM.set_dmk(1, 0, i, j, dmk_value);
        }
    }

    // build the real-space DMR
    Grid_Driver gd(0, 0);
    DM.init_dmr(&gd, &ucell);
    // Gamma-only: reduce R vectors to (0, 0, 0), as cal_DMR_blas_double does
    DM.get_dmr_ptr(1)->fix_gamma();
    DM.cal_dmr(-1);

    // check the Gamma (R = 0) block: rho_0 must be 2a (Pauli), NOT a (real projection);
    // rho_x = rho_y = rho_z = 0 for uu == dd and zero spin off-diagonals.
    hamilt::HContainer<double>* dmr = DM.get_dmr_ptr(1);
    for (int i = 0; i < dmr->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& ap = dmr->get_atom_pair(i);
        double* rho = ap.get_HR_values(0, 0, 0).get_pointer();
        const int col_size = ap.get_col_size();
        const int row_size = ap.get_row_size();
        // walk the 2x2 spin blocks (step_trace = {0, 1, col_size, col_size+1})
        for (int irow = 0; irow < row_size; irow += 2)
        {
            for (int icol = 0; icol < col_size; icol += 2)
            {
                const double* blk = rho + irow * col_size + icol;
                const double rho_0 = blk[0];            // step_trace[0]
                const double rho_x = blk[1];            // step_trace[1]
                const double rho_y = blk[col_size];     // step_trace[2]
                const double rho_z = blk[col_size + 1]; // step_trace[3]
                EXPECT_NEAR(rho_0, 2.0 * a, 1e-10)
                    << "rho_0 wrong: cal_dmr did NOT take the SOC Pauli branch (nspin_global==4)";
                EXPECT_NEAR(rho_x, 0.0, 1e-10);
                EXPECT_NEAR(rho_y, 0.0, 1e-10);
                EXPECT_NEAR(rho_z, 0.0, 1e-10);
            }
        }
    }
#ifdef __MPI
    delete pv_soc;
#endif
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

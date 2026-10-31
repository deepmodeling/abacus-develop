#include "gmock/gmock.h"
#include "gtest/gtest.h"

#undef __LCAO

#include "source_base/module_out/filename.h" // mohan add 2025-05-17
#include "source_base/parallel_grid.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_parallel.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_wf/read_wf2rho_pw.h"
#include "source_io/module_wf/write_wfc_pw.h"
#include "source_psi/psi.h"

#ifdef __MPI
#include "mpi.h"
#include "source_base/parallel_global.h"
#include "source_basis/module_pw/test/test_tool.h"
#endif

/// @brief Friend helper to mutate PARAM private members in unit tests.
/// @details Parameter grants friend access to TestParameters so the test can
/// modify input/sys fields without `#define private public`. The class must
/// stay at global scope to match the friend declaration in parameter.h;
/// an anonymous-namespace class would not be the friend.
class TestParameters
{
  public:
    static Input_para& input() { return PARAM.input; }
    static System_para& sys() { return PARAM.sys; }
};

UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
int XC_Functional::func_type = 0;
bool XC_Functional::ked_flag = false;

namespace module_charge
{
void cal_rhog_symm(const int& spin_now,
                   const Charge& CHR,
                   const ModulePW::PW_Basis* rho_basis,
                   ModuleSymmetry::Symmetry& symm)
{
    return;
}
} // namespace module_charge

void cal_ik2iktot(std::vector<int>& ik2iktot, const int& nks, const int& nkstot)
{
    if (TestParameters::input().kpar == 1)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            ik2iktot[ik] = ik;
        }
    }
    else if (TestParameters::input().kpar == 2)
    {
        if (GlobalV::MY_POOL == 0)
        {
            for (int ik = 0; ik < nks; ++ik)
            {
                ik2iktot[ik] = ik;
            }
        }
        else if (GlobalV::MY_POOL == 1)
        {
            for (int ik = 0; ik < nks; ++ik)
            {
                ik2iktot[ik] = ik + 2; // only works for this test
            }
        }
    }
}

/**
 * - Tested Functions:
 *  - write_wfc_pw()
 *  - read_wf2rho_pw()
 */

class ReadWfcRhoTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis_K* wfcpw = nullptr;
    ModulePW::PW_Basis* rhopw = nullptr;
    K_Vectors* kv = nullptr;
    psi::Psi<std::complex<double>>* psi = nullptr;
    Charge chg;
    ModuleSymmetry::Symmetry symm;
    Parallel_Grid pgrid;
    int world_size = 1;
    virtual void SetUp()
    {
        wfcpw = new ModulePW::PW_Basis_K;
        rhopw = new ModulePW::PW_Basis;
        kv = new K_Vectors;
        // output .dat file
        TestParameters::input().out_wfc_pw = 2;
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    }
    virtual void TearDown()
    {
        delete wfcpw;
        delete rhopw;
        delete kv;
    }
};

TEST_F(ReadWfcRhoTest, ReadWfcRho)
{
    // kpar=1, if nproc=1
    // kpar=2, if nproc>1
    const int kpar = GlobalV::KPAR;
    const int nspin = 1;
    const int nbands = 4;
    const int my_pool = GlobalV::MY_POOL;
    const int my_rank = GlobalV::MY_RANK;
    const int nks = 2;
    const int nkstot = GlobalV::KPAR * nks;

    //-------------------------
    // Initialize the k-points
    //-------------------------
    kv->set_nkstot(nkstot);
    kv->set_nks(nks);
    kv->isk = {0, 0};
    const double shift = my_pool * 0.1;
    kv->kvec_d = {ModuleBase::Vector3<double>(shift, shift, shift), ModuleBase::Vector3<double>(0.5 + shift, 0.5 + shift, 0.5 + shift)};
    kv->ik2iktot.resize(nks);
    cal_ik2iktot(kv->ik2iktot, nks, nkstot);

    //-------------------------
    // Initialize the pw basis
    //-------------------------
#ifdef __MPI
    wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
    rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    rhopw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 80);
    rhopw->initparameters(false, 80);
    rhopw->setuptransform();
    rhopw->collect_local_pw();
    pgrid.init(rhopw->nx, rhopw->ny, rhopw->nz, rhopw->nplane, rhopw->nrxx, rhopw->nz, 1, world_size);

    wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 80);
    wfcpw->initparameters(false, 20, nks, kv->kvec_d.data());
    wfcpw->setuptransform();
    wfcpw->collect_local_pw();

    //-------------------------
    // Initialize k points
    //-------------------------
    kv->kvec_c.clear();
    for (int ik = 0; ik < nks; ++ik)
    {
        kv->kvec_c.push_back(wfcpw->kvec_c[ik]);
    }
    kv->ngk = {wfcpw->npwk[0], wfcpw->npwk[1]};
    kv->wk = {1.0, 1.0};

    //----------------------------------------
    // Initialize weights for wave functions
    //----------------------------------------
    ModuleBase::matrix wg(nkstot, nbands);
    wg.fill_out(1.0);
    if (GlobalV::MY_RANK == 0)
    {
        // mohan update 2025-06-02
        std::ofstream ofs("eig.txt");
        ofs << " Electronic state energy (eV) and occupations" << std::endl;
        ofs << " Spin number " << nspin << std::endl;
        ofs << std::setprecision(8);
        ofs << std::setiosflags(std::ios::showpoint);

        const int is = 0; // nspin is 1
        for (int ik = 0; ik < nkstot; ++ik)
        {
            ofs << " spin=" << is + 1 << " k-point=" << ik + 1 << "/" << nkstot << " Cartesian=" << kv->kvec_c[ik].x << " "
                << kv->kvec_c[ik].y << " " << kv->kvec_c[ik].z << " (" << kv->ngk[ik] << " plane wave)" << std::endl;

            ofs << std::setprecision(16);
            ofs << std::setiosflags(std::ios::showpoint);

            double ekb = -1.23456; // energy
            for (int ib = 0; ib < nbands; ib++)
            {
                ofs << " " << ib + 1 << " " << ekb << " " << wg(ik, ib) << std::endl;
            }

            ofs << std::endl;
        }
        ofs.close();
    }

    //----------------------------------------
    // Initialize wave functions Psi
    //----------------------------------------
    psi = new psi::Psi<std::complex<double>>(nks, nbands, wfcpw->npwk_max, kv->ngk, true);
    std::complex<double>* ptr = psi->get_pointer();
    for (int i = 0; i < nks * nbands * wfcpw->npwk_max; i++)
    {
        ptr[i] = std::complex<double>((i + my_pool * 100) / 100.0, (i + my_pool) / 100.0);
    }

    //----------------------------------------
    // Initialize charge density
    //----------------------------------------
    chg.set_rhopw(rhopw);
    chg.allocate(nspin, false, false, 0);
    chg.pgrid = &pgrid;

    //----------------------------------------
    // set charge_ref
    //----------------------------------------
    Charge chg_ref;
    chg_ref.set_rhopw(rhopw);
    chg_ref.allocate(nspin, false, false, 0);
    chg_ref.pgrid = &pgrid;
    std::vector<std::complex<double>> rho_tmp(rhopw->nrxx);

    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const std::complex<double>* wfc_ib = ptr + ik * nbands * wfcpw->npwk_max + ib * wfcpw->npwk_max;
            wfcpw->recip2real(wfc_ib, rho_tmp.data(), ik);

            const double w1 = wg(ik, ib) / wfcpw->omega;

            for (int ir = 0; ir < rhopw->nrxx; ir++)
            {
                chg_ref.rho[0][ir] += w1 * std::norm(rho_tmp[ir]);
            }
        }
    }

#ifdef __MPI
    module_charge::reduce_diff_pools(chg_ref.rho[0], chg_ref, GlobalV::KPAR,
                                     PARAM.globalv.all_ks_run, PARAM.inp.bndpar);
#endif

    // for spin=1 or 2, npol=1
    const int npol = 1;

    // Write the wave functions to file
    const std::string out_dir = "./";

    // Read the wave functions to charge density
    std::stringstream ss;
    ss << "running_log" << GlobalV::MY_RANK << ".txt";
    std::ofstream running_log(ss.str().c_str());

    running_log << " rank=" << GlobalV::MY_RANK << std::endl;

    const double ecutwfc = 20; // this is a fake number

    const int istep = -1; // -1 means ionic iteration number will not appear in file name
    const int iter = -1;  // -1 means electronic iteration number will not appear in file name
    const int append_step = -1; // single-point write: truncate, do not append

    ModuleIO::write_wfc_pw(istep,
                           iter,
                           append_step,
                           kpar,
                           my_pool,
                           my_rank,
                           nbands,
                           nspin,
                           npol,
                           GlobalV::RANK_IN_POOL,
                           GlobalV::NPROC_IN_POOL,
                           TestParameters::input().out_wfc_pw,
                           0,
                           ecutwfc,
                           out_dir,
                           *psi,
                           *kv,
                           wfcpw,
                           running_log);

    ModuleIO::read_wf2rho_pw(wfcpw,
                             symm,
                             chg,
                             out_dir,
                             kpar,
                             my_pool,
                             my_rank,
                             GlobalV::NPROC_IN_POOL,
                             GlobalV::RANK_IN_POOL,
                             nbands,
                             nspin,
                             npol,
                             nkstot,
                             kv->ik2iktot,
                             kv->isk,
                             running_log);

    // compare the charge density
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        EXPECT_NEAR(chg.rho[0][ir], chg_ref.rho[0][ir], 1e-8);
    }

    if (GlobalV::NPROC == 1)
    {
        EXPECT_NEAR(chg.rho[0][0], 8617.076357957576, 1e-8);
    }
    else if (GlobalV::NPROC == 4)
    {
        const std::vector<double> ref = {8207.849135313403, 35.34776105132742, 8207.849135313403, 35.34776105132742};
        EXPECT_NEAR(chg.rho[0][0], ref[GlobalV::MY_RANK], 1e-8);
        // for (int ip = 0; ip < GlobalV::NPROC; ++ip)
        // {
        //     if (GlobalV::MY_RANK == ip)
        //     {
        //         std::cout.precision(16);
        //         std::cout << GlobalV::MY_RANK << " " << chg.rho[0][0] << std::endl;
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        // }
    }

    delete psi;

    if (GlobalV::MY_RANK == 0)
    {
        remove("running_log0.txt");
        remove("eig.txt");
        remove("wfs1k1_pw.dat");
        remove("wfs1k2_pw.dat");
        if (GlobalV::KPAR == 2)
        {
            remove("wfs1k3_pw.dat");
            remove("wfs1k4_pw.dat");
            remove("running_log1.txt");
            remove("running_log2.txt");
            remove("running_log3.txt");
        }
    }
}

// Size of a file in bytes; returns -1 when the file cannot be opened.
static long file_size_bytes(const std::string& path)
{
    std::ifstream ifs(path.c_str(), std::ios::binary | std::ios::ate);
    if (!ifs.is_open())
    {
        return -1;
    }
    return static_cast<long>(ifs.tellg());
}

// In append mode (out_app_flag == 1) with a fixed filename (istep == -1), the
// first output step of a run must truncate the file while a later step must
// append a new frame. This is driven by append_step, independent of the
// filename sentinel istep. Before the fix, istep == -1 always truncated, so a
// second write never grew the file.
TEST_F(ReadWfcRhoTest, WriteWfcAppendMode)
{
    const int kpar = GlobalV::KPAR;
    const int nspin = 1;
    const int nbands = 4;
    const int my_pool = GlobalV::MY_POOL;
    const int my_rank = GlobalV::MY_RANK;
    const int nproc_in_pool = GlobalV::NPROC_IN_POOL;
    const int rank_in_pool = GlobalV::RANK_IN_POOL;
    const int nks = 2;
    const int nkstot = kpar * nks;

    // K points and PW basis
    kv->set_nkstot(nkstot);
    kv->set_nks(nks);
    kv->isk = {0, 0};
    const double shift = my_pool * 0.1;
    kv->kvec_d = {ModuleBase::Vector3<double>(shift, shift, shift), ModuleBase::Vector3<double>(0.5 + shift, 0.5 + shift, 0.5 + shift)};
    kv->ik2iktot.resize(nks);
    cal_ik2iktot(kv->ik2iktot, nks, nkstot);

#ifdef __MPI
    wfcpw->initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 80);
    wfcpw->initparameters(false, 20, nks, kv->kvec_d.data());
    wfcpw->setuptransform();
    wfcpw->collect_local_pw();

    kv->kvec_c.clear();
    for (int ik = 0; ik < nks; ++ik)
    {
        kv->kvec_c.push_back(wfcpw->kvec_c[ik]);
    }
    kv->ngk = {wfcpw->npwk[0], wfcpw->npwk[1]};
    kv->wk = {1.0, 1.0};

    psi = new psi::Psi<std::complex<double>>(nks, nbands, wfcpw->npwk_max, kv->ngk, true);
    std::complex<double>* ptr = psi->get_pointer();
    for (int i = 0; i < nks * nbands * wfcpw->npwk_max; i++)
    {
        ptr[i] = std::complex<double>((i + my_pool * 100) / 100.0, (i + my_pool) / 100.0);
    }

    const std::string out_dir = "./";
    std::ofstream running_log("running_log_app0.txt");
    const double ecutwfc = 20;
    const int npol = 1;
    const int istep = -1;   // filename sentinel: no step suffix
    const int iter = -1;
    const int out_app_flag = 1; // append mode

    auto write_once = [&](const int append_step) {
        ModuleIO::write_wfc_pw(istep,
                               iter,
                               append_step,
                               kpar,
                               my_pool,
                               my_rank,
                               nbands,
                               nspin,
                               npol,
                               rank_in_pool,
                               nproc_in_pool,
                               TestParameters::input().out_wfc_pw,
                               out_app_flag,
                               ecutwfc,
                               out_dir,
                               *psi,
                               *kv,
                               wfcpw,
                               running_log);
    };

    write_once(0); // fresh first step: truncate
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const long size_first = file_size_bytes("wfk1_pw.dat");

    write_once(1); // later step: append a frame
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const long size_second = file_size_bytes("wfk1_pw.dat");

    write_once(0); // fresh run again: truncate back to a single frame
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const long size_third = file_size_bytes("wfk1_pw.dat");

    if (my_rank == 0)
    {
        EXPECT_GT(size_first, 0);
        EXPECT_GT(size_second, size_first);
        EXPECT_EQ(size_third, size_first);
    }

    running_log.close();
    delete psi;

    if (my_rank == 0)
    {
        remove("running_log_app0.txt");
        remove("wfk1_pw.dat");
        remove("wfk2_pw.dat");
        if (kpar == 2)
        {
            remove("wfk3_pw.dat");
            remove("wfk4_pw.dat");
        }
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    setupmpi(argc, argv, GlobalV::NPROC, GlobalV::MY_RANK);

    // when kpar == 2, nspin == 2
    TestParameters::input().kpar = (GlobalV::NPROC > 1) ? 2 : 1;
    GlobalV::KPAR = TestParameters::input().kpar;
    TestParameters::input().bndpar = 1;
    Parallel_Global::divide_pools(GlobalV::NPROC,
                                  GlobalV::MY_RANK,
                                  TestParameters::input().bndpar,
                                  GlobalV::KPAR,
                                  GlobalV::NPROC_IN_BNDGROUP,
                                  GlobalV::RANK_IN_BPGROUP,
                                  GlobalV::MY_BNDGROUP,
                                  GlobalV::NPROC_IN_POOL,
                                  GlobalV::RANK_IN_POOL,
                                  GlobalV::MY_POOL);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    finishmpi();
#endif
    return result;
}

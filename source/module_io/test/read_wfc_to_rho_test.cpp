#include "gmock/gmock.h"
#include "gtest/gtest.h"

#undef __LCAO

#define private public
#include "module_elecstate/module_charge/charge.h"
#include "module_parameter/parameter.h"
#undef private

#ifdef __MPI
#include "module_base/parallel_global.h"
#include "module_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate_getters.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/read_wfc_to_rho.h"
#include "module_io/write_wfc_pw.h"
#include "module_psi/psi.h"

Parallel_Grid::Parallel_Grid()
{
}
Parallel_Grid::~Parallel_Grid()
{
}
Charge::Charge()
{
}
Charge::~Charge()
{
}
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
int elecstate::get_xc_func_type()
{
    return 0;
}
int XC_Functional::get_func_type()
{
    return 0;
}
Symmetry_rho::Symmetry_rho()
{
}
Symmetry_rho::~Symmetry_rho()
{
}
void Symmetry_rho::begin(const int& spin_now,
                         const Charge& CHR,
                         const ModulePW::PW_Basis* rho_basis,
                         ModuleSymmetry::Symmetry& symm) const
{
    return;
}

void cal_ik2iktot(std::vector<int>& ik2iktot, const int& nks, const int& nkstot)
{
    for (int ik = 0; ik < nks; ++ik)
    {
        int nkp = nkstot / PARAM.inp.kpar;
        int rem = nkstot % PARAM.inp.kpar;
        if (GlobalV::MY_POOL < rem)
        {
            ik2iktot[ik] = GlobalV::MY_POOL * nkp + GlobalV::MY_POOL + ik;
        }
        else
        {
            ik2iktot[ik] = GlobalV::MY_POOL * nkp + rem + ik;
        }
    }
}

namespace GlobalC
{
Parallel_Grid Pgrid;
} // namespace GlobalC

/**
 * - Tested Functions:
 *  - write_wfc_pw()
 *  - read_wfc_to_rho()
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
    virtual void SetUp()
    {
        wfcpw = new ModulePW::PW_Basis_K;
        rhopw = new ModulePW::PW_Basis;
        kv = new K_Vectors;
        PARAM.input.nbands = 4;
        PARAM.input.nspin = 1;
        PARAM.input.out_wfc_pw = 2;
    }
    virtual void TearDown()
    {
        delete wfcpw;
        delete rhopw;
        delete kv;
    }
};

// Test the read_wfc_pw function
TEST_F(ReadWfcRhoTest, ReadWfcRho)
{
    // Init K_Vectors
    const int my_pool = GlobalV::MY_POOL;
    const int nbands = PARAM.input.nbands;
    const int nks = 2;
    const int nkstot = GlobalV::KPAR * nks;
    kv->set_nkstot(nkstot);
    kv->set_nks(nks);
    kv->isk = {0, 0};
    const double shift = my_pool * 0.1;
    kv->kvec_d = {ModuleBase::Vector3<double>(shift, shift, shift),
                  ModuleBase::Vector3<double>(0.5 + shift, 0.5 + shift, 0.5 + shift)};
    kv->ik2iktot.resize(nks);
    cal_ik2iktot(kv->ik2iktot, nks, nkstot);

    // Init the pw basis
#ifdef __MPI
    wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
    rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    rhopw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 80);
    rhopw->initparameters(false, 80);
    rhopw->setuptransform();
    rhopw->collect_local_pw();

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

    // Init wg
    ModuleBase::matrix wg(nkstot, nbands);
    wg.fill_out(1.0);
    if (GlobalV::MY_RANK == 0)
    {
        std::ofstream ofs("istate.info");
        for (int ik = 0; ik < nkstot; ++ik)
        {
            ofs << "BAND               Energy(ev)               Occupation                Kpoint" << std::endl;
            for (int ib = 0; ib < nbands; ++ib)
            {
                ofs << "  " << ib + 1 << "                  0.0000000                " << 1.0 << std::endl;
            }
            ofs << std::endl;
        }
        ofs.close();
    }

    // Init Psi
    psi = new psi::Psi<std::complex<double>>(nks, nbands, wfcpw->npwk_max, kv->ngk, true);
    std::complex<double>* ptr = psi->get_pointer();
    for (int i = 0; i < nks * nbands * wfcpw->npwk_max; i++)
    {
        ptr[i] = std::complex<double>((i + my_pool * 100) / 100.0, (i + my_pool) / 100.0);
    }

    // Init charge
    chg.rho = new double*[1];
    chg._space_rho = new double[rhopw->nrxx];
    chg.rho[0] = chg._space_rho;
    ModuleBase::GlobalFunc::ZEROS(chg.rho[0], rhopw->nrxx);
    chg.rhopw = rhopw;
    chg.nrxx = rhopw->nrxx;
    // set charge_ref
    Charge chg_ref;
    chg_ref.rho = new double*[1];
    chg_ref._space_rho = new double[rhopw->nrxx];
    chg_ref.rho[0] = chg_ref._space_rho;
    ModuleBase::GlobalFunc::ZEROS(chg_ref.rho[0], rhopw->nrxx);
    std::vector<std::complex<double>> rho_tmp(rhopw->nrxx);
    chg_ref.nrxx = rhopw->nrxx;

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
    chg_ref.init_chgmpi();
    chg_ref.reduce_diff_pools(chg_ref.rho[0]);
#endif

    // Write the wave functions to file
    ModuleIO::write_wfc_pw("WAVEFUNC", *psi, *kv, wfcpw);

    // Read the wave functions to charge density
    ModuleIO::read_wfc_to_rho(wfcpw, symm, kv->ik2iktot.data(), nkstot, kv->isk, chg);

    // compare the charge density
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        EXPECT_NEAR(chg.rho[0][ir], chg_ref.rho[0][ir], 1e-8);
    }
    // std::cout.precision(16);
    // std::cout<<chg.rho[0][0]<<std::endl;
    if (GlobalV::NPROC == 1) {
        EXPECT_NEAR(chg.rho[0][0], 8617.076357957576, 1e-8);
    } else if (GlobalV::NPROC == 4)
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

    delete[] chg.rho;
    delete[] chg._space_rho;
    delete[] chg_ref.rho;
    delete[] chg_ref._space_rho;
    delete psi;
    if (GlobalV::MY_RANK == 0)
    {
        remove("istate.info");
        remove("WAVEFUNC1.dat");
        remove("WAVEFUNC2.dat");
        if (GlobalV::KPAR > 1)
        {
            remove("WAVEFUNC3.dat");
            remove("WAVEFUNC4.dat");
        }
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    setupmpi(argc, argv, GlobalV::NPROC, GlobalV::MY_RANK);
    PARAM.input.kpar = (GlobalV::NPROC > 1) ? 2 : 1;
    GlobalV::KPAR = PARAM.input.kpar;
    PARAM.input.bndpar = 1;
    Parallel_Global::divide_pools(GlobalV::NPROC,
                                  GlobalV::MY_RANK,
                                  PARAM.inp.bndpar,
                                  GlobalV::KPAR,
                                  GlobalV::NPROC_IN_STOGROUP,
                                  GlobalV::RANK_IN_STOGROUP,
                                  GlobalV::MY_STOGROUP,
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
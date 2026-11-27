#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/magnetism.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_symm.h"

#include <complex>
#include <vector>

// charge.cpp references Magnetism; provide a lightweight stub.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_symm.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - symmetrize_rho: dispatch per nspin to cal_rhog_symm / cal_rhog_symm_soc
 *   - cal_rhog_symm: no-op when symm_flag != 1; otherwise FFT + psymmg + FFT back
 *   - cal_rhog_symm (raw array overload): same no-op behavior
 *   - cal_rhog_symm_soc: no-op when symm_flag != 1
 *
 * The no-op paths (symm_flag == 0) are fully covered without a real symmetry
 * group: the functions return immediately, leaving rho/rhog unchanged.
 */

class ChgSymmTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis pw_basis;
    Charge charge;

    void SetUp() override
    {
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
    }

    void setup_charge(int nspin)
    {
        charge.set_rhopw(&pw_basis);
        const bool kin_den = false;
        const bool meta_gga = false;
        charge.allocate(nspin, kin_den, meta_gga, 0);
    }
};

// ---------------------------------------------------------------------------
// no-op path: symm_flag == 0 leaves density untouched
// ---------------------------------------------------------------------------

TEST_F(ChgSymmTest, SymmetrizeRhoSymmFlagOffIsNoopNspin1)
{
    const int nspin = 1;
    setup_charge(nspin);
    ModuleSymmetry::Symmetry symm;
    ModuleSymmetry::Symmetry::symm_flag = 0;

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        charge.rho[0][ir] = static_cast<double>(ir + 1);
    }
    std::vector<double> rho_before(charge.rho[0], charge.rho[0] + pw_basis.nrxx);

    module_charge::symmetrize_rho(nspin, charge, &pw_basis, symm);

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_EQ(charge.rho[0][ir], rho_before[ir]);
    }
}

TEST_F(ChgSymmTest, SymmetrizeRhoSymmFlagOffIsNoopNspin4)
{
    const int nspin = 4;
    setup_charge(nspin);
    ModuleSymmetry::Symmetry symm;
    ModuleSymmetry::Symmetry::symm_flag = 0;

    std::vector<std::vector<double>> rho_before(nspin);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw_basis.nrxx; ++ir)
        {
            charge.rho[is][ir] = static_cast<double>(is * 100 + ir);
        }
        rho_before[is].assign(charge.rho[is], charge.rho[is] + pw_basis.nrxx);
    }

    module_charge::symmetrize_rho(nspin, charge, &pw_basis, symm);

    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw_basis.nrxx; ++ir)
        {
            EXPECT_EQ(charge.rho[is][ir], rho_before[is][ir]);
        }
    }
}

TEST_F(ChgSymmTest, CalRhogSymmRawArrayNoop)
{
    const int nspin = 1;
    setup_charge(nspin);
    ModuleSymmetry::Symmetry symm;
    ModuleSymmetry::Symmetry::symm_flag = 0;

    std::vector<double> rho_buf(pw_basis.nrxx, 3.0);
    std::vector<std::complex<double>> rhog_buf(pw_basis.npw, std::complex<double>(0.0, 0.0));
    double* rho_ptrs[1] = {rho_buf.data()};
    std::complex<double>* rhog_ptrs[1] = {rhog_buf.data()};

    module_charge::cal_rhog_symm(0, rho_ptrs, rhog_ptrs, pw_basis.npw, nullptr, &pw_basis, symm);

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_EQ(rho_buf[ir], 3.0);
    }
}

TEST_F(ChgSymmTest, CalRhogSymmSocNoop)
{
    const int nspin = 4;
    setup_charge(nspin);
    ModuleSymmetry::Symmetry symm;
    ModuleSymmetry::Symmetry::symm_flag = 0;

    std::vector<std::vector<double>> rho_buf(nspin, std::vector<double>(pw_basis.nrxx, 1.0));
    std::vector<std::vector<std::complex<double>>> rhog_buf(
        nspin, std::vector<std::complex<double>>(pw_basis.npw, std::complex<double>(0.0, 0.0)));
    double* rho_ptrs[4] = {rho_buf[0].data(), rho_buf[1].data(), rho_buf[2].data(), rho_buf[3].data()};
    std::complex<double>* rhog_ptrs[4] = {rhog_buf[0].data(), rhog_buf[1].data(),
                                          rhog_buf[2].data(), rhog_buf[3].data()};

    module_charge::cal_rhog_symm_soc(rho_ptrs, rhog_ptrs, &pw_basis, symm);

    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw_basis.nrxx; ++ir)
        {
            EXPECT_EQ(rho_buf[is][ir], 1.0);
        }
    }
}

#include "source_pw/module_pwdft/op_pw_meta.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

int XC_Functional::func_type = 1;
bool XC_Functional::ked_flag = false;

TEST(MetaPW, ConstantTauPotentialActsOnMetaGgaAndHybridMetaGga)
{
    ModulePW::PW_Basis_K basis("cpu", "double");
    const ModuleBase::Matrix3 lattice(1, 0, 0, 0, 1, 0, 0, 0, 1);
    const ModuleBase::Vector3<double> kpoint(0.2, 0.1, 0.0);
    basis.initgrids(6.0, lattice, 6, 6, 6);
    basis.initparameters(false, 3.0, 1, &kpoint, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
    const int npw = basis.npwk[0];
    const int spin = 0;
    const double vtau = 0.7;
    std::vector<double> potential(basis.nrxx, vtau);
    const std::complex<double> coefficient(0.3, -0.4);
    const std::complex<double> initial(1.0, 0.5);
    std::vector<std::complex<double>> wavefunction(npw, coefficient);
    std::vector<std::complex<double>> hpsi(npw, initial);
    hamilt::Meta<hamilt::OperatorPW<std::complex<double>>> meta(
        basis.tpiba, &spin, potential.data(), 1, basis.nrxx, &basis);
    meta.init(0);

    // A constant v_tau gives -div(v_tau grad) psi = v_tau |G+k|^2 psi.
    // Types 1/2/4 must leave the accumulated Hamiltonian unchanged.
    for (int type = 1; type <= 5; ++type)
    {
        const bool needs_tau = type == 3 || type == 5;
        XC_Functional::set_func_type(type);
        XC_Functional::set_ked_flag(needs_tau);
        hpsi.assign(npw, initial);
        meta.act(1, npw, 1, wavefunction.data(), hpsi.data(), npw, false);
        for (int ig = 0; ig < npw; ++ig)
        {
            const double multiplier = vtau * basis.tpiba2 * basis.gk2[ig];
            const std::complex<double> contribution = needs_tau ? multiplier * coefficient : 0.0;
            const std::complex<double> expected = initial + contribution;
            EXPECT_NEAR(hpsi[ig].real(), expected.real(), 1.0e-12) << type << ":" << ig;
            EXPECT_NEAR(hpsi[ig].imag(), expected.imag(), 1.0e-12) << type << ":" << ig;
        }
    }
    XC_Functional::set_func_type(1);
    XC_Functional::set_ked_flag(false);
}

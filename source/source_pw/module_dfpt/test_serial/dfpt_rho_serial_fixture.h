#ifndef DFPT_RHO_SERIAL_FIXTURE_H
#define DFPT_RHO_SERIAL_FIXTURE_H

// Shared serial-side fixture for the DFPT_Rho test files split along the
// mixing/response boundary (test_dfpt_rho.cpp, test_dfpt_rho_mixing.cpp).
// Everything runs without __MPI on the DFPTSerialBase shared-grid wiring.

#include "dfpt_serial_fixture.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"
#include "source_pw/module_dfpt/dfpt_rho.h"

#include <cmath>
#include <complex>
#include <vector>

namespace dfpt_rho_test_rng
{

// C++11-safe shared seed: a function-local static inside an inline
// function has exactly one instance across the translation units of the
// test target (an inline VARIABLE would require C++17; the unit-test
// configuration only raises the language level to C++14)
inline unsigned& g_seed()
{
    static unsigned seed = 20260815u;
    return seed;
}

inline double test_rand()
{
    unsigned& seed = g_seed();
    seed = seed * 1664525u + 1013904223u;
    return ((seed >> 8) & 0xffffff) / 16777216.0 * 2.0 - 1.0;
}

inline std::complex<double> crand()
{
    return std::complex<double>(test_rand(), test_rand());
}

} // namespace dfpt_rho_test_rng

class DFPTRhoSerialTest : public DFPTSerialBase
{
  protected:
    ModuleDFPT::DFPT_Rho rho_;

    static const int nbands_ = 2;

    void SetUp() override
    {
        DFPTSerialBase::SetUp();
        rho_.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_plain"), 0.4, 0.0});
    }

    void FillRandomStates(psi::Psi<std::complex<double>>& psi,
                          std::vector<std::vector<std::vector<std::complex<double>>>>& dpsi)
    {
        using namespace dfpt_rho_test_rng;
        psi::Psi<std::complex<double>> p(1, nbands_, pw_wfc_.npwk_max, pw_wfc_.npwk[0], true);
        ModuleDFPT::DFPT_KQ_Basis kq;
        kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
        dpsi.assign(1, std::vector<std::vector<std::complex<double>>>(nbands_));
        for (int ib = 0; ib < nbands_; ++ib)
        {
            for (int igl = 0; igl < pw_wfc_.npwk[0]; ++igl)
            {
                p(0, ib, igl) = crand();
            }
            dpsi[0][ib].assign(kq.get_npwk(), std::complex<double>(0.0, 0.0));
            for (int igl = 0; igl < kq.get_npwk(); ++igl)
            {
                dpsi[0][ib][igl] = crand();
            }
            data_.set_dpsi(0, 0, ib, dpsi[0][ib]);
        }
        psi = p;
    }
};

#endif // DFPT_RHO_SERIAL_FIXTURE_H

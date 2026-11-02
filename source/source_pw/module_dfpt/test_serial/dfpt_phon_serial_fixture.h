#ifndef DFPT_PHON_SERIAL_FIXTURE_H
#define DFPT_PHON_SERIAL_FIXTURE_H

// Shared serial-side fixture for the DFPT_Phon test files split along the
// source boundaries (test_dfpt_phon{,_ewald,_elec,_asr}.cpp). Everything
// runs without __MPI on the DFPTSerialBase cell/basis/data wiring.

#include "dfpt_serial_fixture.h"
#include "source_base/complexmatrix.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_phon.h"
#include "source_cell/stru_fac.h"

class DFPTPhonSerialTest : public DFPTSerialBase
{
  protected:
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Phon phon_;

    void SetUp() override
    {
        DFPTSerialBase::SetUp();
        SetupPhon(k_d_, q_d_);
    }

    // (re)initialize the bases and the pert/phon wiring for a given (k, q)
    // pair; SetUp uses the default fixture values
    void SetupPhon(const ModuleBase::Vector3<double>& k_d, const ModuleBase::Vector3<double>& q_d)
    {
        SetupBases(k_d, q_d, 2);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        phon_.init(ucell_, &pw_rho_, &pert_);
    }

    // independent Ry/bohr^2/amu -> cm^-1 conversion used by diagonalize
    double RyBohr2AmuToCm1() const
    {
        const double amu_kg = 1.66053906660e-27; // CODATA amu in kg
        return std::sqrt(ModuleBase::RYDBERG_SI / amu_kg) / (0.529177210903e-10 * 2.0 * ModuleBase::PI * 2.99792458e10);
    }

    // common setup of the isotropic loto closed-form tests: zero 6x6
    // dynamical matrix, eps_inf = 3I, Z*_1 = 1, Z*_2 = 2, and the
    // two-atom 12/4 mass table for the mass lookup
    void SetupIsotropicLoto()
    {
        data_.set_dynmat(0, ModuleBase::ComplexMatrix(6, 6, true));
        ModuleBase::matrix eps(3, 3, true);
        for (int d = 0; d < 3; ++d)
        {
            eps(d, d) = 3.0;
        }
        data_.set_dielectric(eps);
        ModuleBase::matrix z1(3, 3, true);
        ModuleBase::matrix z2(3, 3, true);
        z1(0, 0) = z1(1, 1) = z1(2, 2) = 1.0;
        z2(0, 0) = z2(1, 1) = z2(2, 2) = 2.0;
        data_.set_born(0, z1);
        data_.set_born(1, z2);
        MakeTwoAtomCell();
    }

    // band 0 occupied with wg = 2, band 1 unoccupied
    static ModuleBase::matrix MakeOccWeights()
    {
        ModuleBase::matrix wg(1, 2, true);
        wg(0, 0) = 2.0;
        wg(0, 1) = 0.0;
        return wg;
    }

    // psi: band 0 = single plane wave at G' = 0 (c = 1, occupied),
    // band 1 unoccupied. The buffer is allocated uninitialized; zero it so
    // only the component set below is nonzero regardless of heap history
    // from earlier tests. getgpluskcar returns the cartesian k+G, so the
    // G = 0 entry is the one with k+G = k_cart.
    psi::Psi<std::complex<double>> MakeSinglePlaneWavePsi(const ModuleBase::Vector3<double>& k_d)
    {
        const int npwk = pw_wfc_.npwk[0];
        psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
        psi.zero_out();
        const ModuleBase::Vector3<double> k_cart = k_d * ucell_.G;
        for (int ig = 0; ig < npwk; ++ig)
        {
            const ModuleBase::Vector3<double> gk = pw_wfc_.getgpluskcar(0, ig);
            if (std::abs(gk.x - k_cart.x) < 1e-10 && std::abs(gk.y - k_cart.y) < 1e-10
                && std::abs(gk.z - k_cart.z) < 1e-10)
            {
                psi(0, 0, ig) = std::complex<double>(1.0, 0.0);
                break;
            }
        }
        return psi;
    }

    // analytic cross term sum_G'' conj(dpsi_G'') sum_i c_i
    // RHS^a(G''-G'_i): the first-order local Coulomb potential on the k+q
    // basis, RHS^a(w) = -i tpiba w_a Vloc(|w|^2) e^{-i 2pi w.tau} with
    // w = G'' - G' + q (GS structure-factor phase convention); the
    // Delta + q = 0 component is dropped by dVloc
    std::complex<double> AnalyticCrossTerm(const ModuleDFPT::DFPT_KQ_Basis& kq,
                                           const std::vector<std::complex<double>>& psi_coef,
                                           const std::vector<ModuleBase::Vector3<double>>& psi_gcart,
                                           const std::vector<std::complex<double>>& dpsi_inj,
                                           int adir) const
    {
        std::complex<double> cross(0.0, 0.0);
        for (int igl = 0; igl < kq.get_npwk(); ++igl)
        {
            const ModuleBase::Vector3<double> gpp = kq.get_gpluskq(igl);
            for (size_t ic = 0; ic < psi_coef.size(); ++ic)
            {
                // AnalyticDVloc returns 0 at w = 0 (dVloc drop)
                cross += psi_coef[ic] * std::conj(dpsi_inj[igl])
                         * AnalyticDVloc(adir, gpp - psi_gcart[ic] + q_cart_);
            }
        }
        return cross;
    }
};

#endif // DFPT_PHON_SERIAL_FIXTURE_H

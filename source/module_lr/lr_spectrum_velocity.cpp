#include "lr_spectrum.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_lr/utils/lr_util_hcontainer.h"
#include "math.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
namespace LR
{
    /// get the velocity matrix v(R)
    inline TD_current get_velocity_matrix_R(const UnitCell& ucell,
        const Grid_Driver& gd,
        const Parallel_Orbitals& pmat,
        const TwoCenterBundle& two_center_bundle)
    {
        // convert the orbital object to the old class for TD_current
        LCAO_Orbitals orb;
        const auto& inp = PARAM.inp;
        two_center_bundle.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);
        // actually this class calculates the velocity matrix v(R) at A=0
        TD_current vR(&ucell, &gd, &pmat, orb, two_center_bundle.overlap_orb.get());
        vR.calculate_vcomm_r(); // $<\mu, 0|[Vnl, r]|\nu, R>$
        vR.calculate_grad_term();   // $<\mu, 0|\nabla|\nu, R>$
        return vR;
    }

    inline double lorentz_delta(const double freq_diff, const double eta)
    {
        return eta / (freq_diff * freq_diff + eta * eta) / M_PI;
    }

    template<typename T> inline ModuleBase::Vector3<T> convert_vector_to_vector3(const std::vector<std::complex<double>>& vec);
    template<> inline ModuleBase::Vector3<double> convert_vector_to_vector3(const std::vector<std::complex<double>>& vec)
    {
        assert(vec.size() == 3);
        return ModuleBase::Vector3<double>(vec[0].real(), vec[1].real(), vec[2].real());
    }
    template<> inline ModuleBase::Vector3<std::complex<double>> convert_vector_to_vector3(const std::vector<std::complex<double>>& vec)
    {
        assert(vec.size() == 3);
        return ModuleBase::Vector3<std::complex<double>>(vec[0], vec[1], vec[2]);
    }

    /// this algorithm has bug in multi-k cases, just for test
    template<typename T>
    ModuleBase::Vector3<T> LR::LR_Spectrum<T>::cal_transition_dipole_istate_velocity_R(const int istate, const TD_current& vR)
    {
        // transition density matrix D(R)
        const elecstate::DensityMatrix<T, T>& DM_trans = this->cal_transition_density_matrix(istate);

        std::vector<std::complex<double>> trans_dipole(3, 0.0);    // $=\sum_{uvR} v(R) D(R) = \sum_{iak}X_{iak}<ck|v|vk>$
        const std::complex<double> fac = ModuleBase::IMAG_UNIT / (eig[istate] / ModuleBase::e2);    // eV to Hartree
        for (int i = 0; i < 3; i++)
        {
            for (int is = 0;is < this->nspin_x; ++is)
            {
                trans_dipole[i] += LR_Util::dot_R_matrix(*vR.get_current_term_pointer(i), *DM_trans.get_DMR_pointer(is + 1), ucell.nat) * fac;
            }   // end for spin_x, only matter in open-shell system
            trans_dipole[i] *= static_cast<double>(this->nk);  // nk is divided inside DM_trans, now recover it
            if (this->nspin_x == 1) { trans_dipole[i] *= sqrt(2.0); } // *2 for 2 spins, /sqrt(2) for the halfed dimension of X in the normalizaiton
            Parallel_Reduce::reduce_all(trans_dipole[i]);
        }   // end for direction
        return convert_vector_to_vector3<T>(trans_dipole);
    }

    // this algorithm is actually in use
    template<typename T>
    ModuleBase::Vector3<T> LR::LR_Spectrum<T>::cal_transition_dipole_istate_velocity_k(const int istate, const TD_current& vR)
    {
        // transition density matrix D(R)
        const elecstate::DensityMatrix<T, T>& DM_trans = this->cal_transition_density_matrix(istate, this->X, false);

        std::vector<std::complex<double>> trans_dipole(3, 0.0);    // $=\sum_{uvR} v(R) D(R) = \sum_{iak}X_{iak}<ck|v|vk>$
        const std::complex<double> fac = ModuleBase::IMAG_UNIT / (eig[istate] / ModuleBase::e2);    // eV to Hartree
        for (int i = 0; i < 3; i++)
        {
            for (int is = 0;is < this->nspin_x;++is)
            {
                for (int ik = 0;ik < nk;++ik)
                {
                    std::vector<std::complex<double>> vk(pmat.get_local_size(), 0.0);
                    hamilt::folding_HR(*vR.get_current_term_pointer(i), vk.data(), kv.kvec_d[ik], pmat.get_row_size(), 1);
                    trans_dipole[i] += std::inner_product(vk.begin(), vk.end(), DM_trans.get_DMK_pointer(is * nk + ik), std::complex<double>(0., 0.)) * fac;
                }
            }   // end for spin_x, only matter in open-shell system
            trans_dipole[i] *= static_cast<double>(this->nk);  // nk is divided inside DM_trans, now recover it
            if (this->nspin_x == 1) { trans_dipole[i] *= sqrt(2.0); } // *2 for 2 spins, /sqrt(2) for the halfed dimension of X in the normalizaiton
            Parallel_Reduce::reduce_all(trans_dipole[i]);
        }   // end for direction
        return convert_vector_to_vector3<T>(trans_dipole);
    }

    template<typename T>
    void LR::LR_Spectrum<T>::cal_transition_dipoles_velocity()
    {
        const TD_current& vR = get_velocity_matrix_R(ucell, gd_, pmat, two_center_bundle_);     // velocity matrix v(R)
        transition_dipole_.resize(nstate);
        this->mean_squared_transition_dipole_.resize(nstate);
        for (int istate = 0;istate < nstate;++istate)
        {
            transition_dipole_[istate] = cal_transition_dipole_istate_velocity_k(istate, vR);
            mean_squared_transition_dipole_[istate] = cal_mean_squared_dipole(transition_dipole_[istate]);
        }
    }

    template<typename T>
    void LR::LR_Spectrum<T>::optical_absorption_method2(const std::vector<double>& freq, const double eta)
    {
        ModuleBase::TITLE("LR::LR_Spectrum", "optical_absorption_velocity");
        // 4*pi^2/V * mean_squared_dipole *delta(w-Omega_S)
        std::ofstream ofs(PARAM.globalv.global_out_dir + "absorption.dat");
        if (GlobalV::MY_RANK == 0) { ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; }
        const double fac = 4 * M_PI * M_PI / ucell.omega * ModuleBase::e2 / this->nk;  // e2: Ry to Hartree in the denominator
        for (int f = 0;f < freq.size();++f)
        {
            double abs_value = 0.0;
            for (int i = 0;i < nstate;++i)
            {
                abs_value += this->mean_squared_transition_dipole_[i] * lorentz_delta((freq[f] - eig[i]), eta);
            }
            abs_value *= fac;
            if (GlobalV::MY_RANK == 0) { ofs << freq[f] * ModuleBase::Ry_to_eV << "\t" << 91.126664 / freq[f] << "\t" << abs_value << std::endl; }
        }
    }

    inline void cal_eig_ks_diff(double* const eig_ks_diff, const double* const eig_ks, const Parallel_2D& px, const int nk, const int nocc, const int nvirt)
    {
        for (int ik = 0;ik < nk;++ik)
        {
            const int& start_k = ik * (nocc + nvirt);
            for (int io = 0;io < px.get_col_size();++io)    //nocc_local
            {
                for (int iv = 0;iv < px.get_row_size();++iv)    //nvirt_local
                {
                    int io_g = px.local2global_col(io);
                    int iv_g = px.local2global_row(iv);
                    eig_ks_diff[ik * px.get_local_size() + io * px.get_row_size() + iv] = (eig_ks[start_k + nocc + iv_g] - eig_ks[start_k + io_g]) / ModuleBase::e2;  // eV to Hartree
                }
            }
        }
    }

    template<typename T>
    void LR::LR_Spectrum<T>::test_transition_dipoles_velocity_ks(const double* const ks_eig)
    {
        // velocity matrix v(R)
        const TD_current& vR = get_velocity_matrix_R(ucell, gd_, pmat, two_center_bundle_);
        //  (e_c-e_v) of KS eigenvalues
        std::vector<double> eig_ks_diff(this->ldim);
        for (int is = 0;is < this->nspin_x;++is)
        {
            cal_eig_ks_diff(eig_ks_diff.data() + is * nk * pX[0].get_local_size(), ks_eig, pX[is], nk, nocc[is], nvirt[is]);
        }
        //  X/(ec-ev)
        std::vector<T> X_div_ks_eig(nstate * this->ldim);
        for (int istate = 0;istate < nstate;++istate)
        {
            const int st = istate * this->ldim;
            std::transform(X + st, X + st + ldim, eig_ks_diff.begin(), X_div_ks_eig.data() + st, std::divides<T>());
        }

        this->transition_dipole_.resize(nstate);
        this->mean_squared_transition_dipole_.resize(nstate);
        for (int istate = 0;istate < nstate;++istate)
        {
            // transition density matrix D(R)
            const elecstate::DensityMatrix<T, T>& DM_trans = this->cal_transition_density_matrix(istate, X_div_ks_eig.data());
            std::vector<std::complex<double>> tmp_trans_dipole(3, 0.0);
            for (int i = 0; i < 3; i++)
            {
                for (int is = 0;is < this->nspin_x; ++is)
                {
                    tmp_trans_dipole[i] += LR_Util::dot_R_matrix(*vR.get_current_term_pointer(i), *DM_trans.get_DMR_pointer(is + 1), ucell.nat) * ModuleBase::IMAG_UNIT;
                }   // end for spin_x, only matter in open-shell system
            }   // end for direction
            this->transition_dipole_[istate] = convert_vector_to_vector3<T>(tmp_trans_dipole);
            this->mean_squared_transition_dipole_[istate] = cal_mean_squared_dipole(transition_dipole_[istate]);
        }
    }
}
template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;
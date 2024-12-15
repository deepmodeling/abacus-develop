#include "lr_spectrum.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_hamilt_lcao/module_tddft/td_current.h"
#include "module_lr/utils/lr_util_hcontainer.h"
#include "math.h"
#include "module_parameter/parameter.h"
namespace LR
{
    /// get the velocity matrix v(R)
    inline TD_current get_velocity_matrix_R(const UnitCell& ucell,
        Grid_Driver& gd,
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

    template<typename T>
    ModuleBase::Vector3<T> LR::LR_Spectrum<T>::cal_transition_dipole_istate_velocity(const int istate)
    {
        // velocity matrix v(R)
        const TD_current& vR = get_velocity_matrix_R(ucell, gd_, pmat, two_center_bundle_);
        // transition density matrix D(R)
        const elecstate::DensityMatrix<T, T>& DM_trans = this->cal_transition_density_matrix(istate);

        std::vector<std::complex<double>> trans_dipole(3, 0.0);    // $=\sum_{uvR} v(R) D(R) = \sum_{iak}X_{iak}<ck|v|vk>$
        double mean_transition_velocity_norm2 = 0.0;    // $= |v|^2/3$
        for (int i = 0; i < 3; i++)
        {
            const std::complex<double> fac = ModuleBase::IMAG_UNIT / (eig[istate] / ModuleBase::e2);    // eV to Hartree
            for (int is = 0;is < this->nspin_x; ++is)
            {
                trans_dipole[i] += LR_Util::dot_R_matrix(*vR.get_current_term_pointer(i), *DM_trans.get_DMR_pointer(is + 1), ucell.nat) * fac;
            }   // end for spin_x, only matter in open-shell system
        }   // end for direction
        return convert_vector_to_vector3<T>(trans_dipole);
    }

    template<typename T>
    void LR::LR_Spectrum<T>::cal_transition_dipoles_velocity()
    {
        transition_dipole_.resize(nstate);
        this->mean_squared_transition_dipole_.resize(nstate);
        for (int istate = 0;istate < nstate;++istate)
        {
            transition_dipole_[istate] = cal_transition_dipole_istate_velocity(istate);
            mean_squared_transition_dipole_[istate] = cal_mean_squared_dipole(transition_dipole_[istate]);
        }
    }

    template<typename T>
    void LR::LR_Spectrum<T>::optical_absorption_method2(const std::vector<double>& freq, const double eta, const std::string& spintype)
    {
        ModuleBase::TITLE("LR::LR_Spectrum", "optical_absorption_velocity");
        // 4*pi^2/V * mean_squared_dipole *delta(w-Omega_S)
        std::ofstream ofs(PARAM.globalv.global_out_dir + "absorption_" + spintype + ".dat");
        if (GlobalV::MY_RANK == 0) { ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; }
        const double fac = 4 * M_PI * M_PI / ucell.omega * ModuleBase::e2;  // e2: Ry to Hartree in the denominator
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
}
template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;
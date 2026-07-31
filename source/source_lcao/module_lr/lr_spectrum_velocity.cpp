#include "lr_spectrum.h"
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "math.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
namespace LR
{
    /// get the velocity matrix v(R)
    inline Velocity_op<std::complex<double>> get_velocity_matrix_R(const UnitCell& ucell,
        const Grid_Driver& gd,
        const Parallel_Orbitals& pmat,
        const TwoCenterBundle& two_center_bundle)
    {
        // convert the orbital object to the old class for Velocity_op
        LCAO_Orbitals orb;
        const auto& inp = PARAM.inp;
        two_center_bundle.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax,
                                           inp.out_element_info, inp.cal_force);
        // actually this class calculates the velocity matrix v(R) at A=0
        Velocity_op<std::complex<double>> vR(&ucell, &gd, &pmat, orb, two_center_bundle.overlap_orb.get());
        vR.calculate_vcomm_r(); // $<\mu, 0|i[Vnl, r]|\nu, R>$
        vR.calculate_grad_term();   // $<\mu, 0|-i∇r|\nu, R>$
        return vR;
    }

    inline double lorentz_delta(const double dfreq_au, const double eta_au)
    {
        return eta_au / (dfreq_au * dfreq_au + eta_au * eta_au) / M_PI;
    }
    inline double gauss_delta(const double dfreq_au, const double eta_au)
    {
        const double c = eta_au / std::sqrt(2. * std::log(2.));
        return std::exp(-dfreq_au * dfreq_au / (2 * c * c)) / (std::sqrt(2 * M_PI) * c);
    }
    template<typename T>
    void LR::LR_Spectrum<T>::optical_absorption_method2(const std::vector<double>& freq, const double eta)
    {
        ModuleBase::TITLE("LR::LR_Spectrum", "optical_absorption_method2");
        // 4*pi^2/V * mean_squared_dipole *delta(w-Omega_S)
        std::ofstream ofs(PARAM.globalv.global_out_dir + "absorption.dat");
        if (GlobalV::MY_RANK == 0) { ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; }
        const double fac = 4 * M_PI * M_PI / ucell.omega / this->nk;
        for (int f = 0;f < freq.size();++f)
        {
            double abs_value = 0.0;
            for (int i = 0;i < nstate;++i)
            {
                abs_value += this->mean_squared_transition_dipole_[i] * lorentz_delta((freq[f] - omega[i]) / ModuleBase::e2, eta / ModuleBase::e2); // e2: Ry to Hartree 
            }
            abs_value *= fac;
            if (GlobalV::MY_RANK == 0) { ofs << freq[f] * ModuleBase::Ry_to_eV << "\t" << 91.126664 / freq[f] << "\t" << abs_value << std::endl; }
        }
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
    inline ModuleBase::Vector3<T> convert_ptr_to_vector3(const std::complex<double>* ptr);
    template<>
    inline ModuleBase::Vector3<double> convert_ptr_to_vector3<double>(const std::complex<double>* ptr)
    {
        return ModuleBase::Vector3<double>(ptr[0].real(), ptr[1].real(), ptr[2].real());
    }
    template<>
    inline ModuleBase::Vector3<std::complex<double>> convert_ptr_to_vector3<std::complex<double>>(const std::complex<double>* ptr)
    {
        return ModuleBase::Vector3<std::complex<double>>(ptr[0], ptr[1], ptr[2]);
    }

    /// this algorithm is stored for reference
    template<typename T>
    ModuleBase::Vector3<T> LR::LR_Spectrum<T>::cal_transition_dipole_istate_velocity_R(const int istate, const Velocity_op<std::complex<double>>& vR)
    {
        // transition density matrix D(R)
        const elecstate::DensityMatrix<T, T>& DM_trans = this->cal_transition_density_matrix(istate);

        std::vector<std::complex<double>> trans_dipole(3, 0.0);    // $=\sum_{uvR} v(R) D(R) = \sum_{aik}X_{aik}<ik|v|ak>$
        const std::complex<double> fac = ModuleBase::IMAG_UNIT / (omega[istate] / ModuleBase::e2);    // Ry to Hartree
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

    // this algorithm is stored for reference
    template<typename T>
    ModuleBase::Vector3<T> LR::LR_Spectrum<T>::cal_transition_dipole_istate_velocity_k(const int istate, const Velocity_op<std::complex<double>>& vR)
    {
        // transition density matrix D(R)
        const elecstate::DensityMatrix<T, T>& DM_trans = this->cal_transition_density_matrix(istate, this->X, false);

        std::vector<std::complex<double>> trans_dipole(3, 0.0);    // $=\sum_{uvk} v(k) D(k) = \sum_{aik}X_{aik}<ik|v|ak>$
        const std::complex<double> fac = ModuleBase::IMAG_UNIT / (omega[istate] / ModuleBase::e2);    // Ry to Hartree
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

    // this algorithm is faster since velocity_mo is calculated and each transition state only need to contract with X
    template<typename T>
    void LR::LR_Spectrum<T>::cal_transition_dipole_istate_velocity_mo(DipoleEnergyType method, const std::vector<double>& eig_ks_diff)
    {
        ModuleBase::timer::start("LR_Spectrum", "cal_transition_dipole_istate_velocity_mo");
        if (this->vmo_ptr == nullptr)
        {
            ModuleBase::WARNING_QUIT("LR_Spectrum", "velocity_mo is null. Please pass a valid pointer.");
        }
        const int nbands = this->nocc[0] + this->nvirt[0];
        assert(nbands == this->pc.get_global_col_size());
        const bool use_ks_gap = (method == DipoleEnergyType::KS_GAP);
        
        std::vector<std::complex<double>> trans_dipole_buf(3 * nstate, 0.0); // $= \sum_{aik} i (<ik|v|ak>X_{aik}/(Ea-Ei) + <ak|v|ik>Y_{aik}/(Ei-Ea))$
        // vmo is global [spin, direction, kpoint, nbands, nbands], X is local [spin, kpoint, nocc_local, nvirt_local]

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int istate = 0; istate < nstate; ++istate)
        {
            const std::complex<double> fac = use_ks_gap ? ModuleBase::IMAG_UNIT :
                                            (ModuleBase::IMAG_UNIT / (omega[istate] / 2.0)); // Ry to Hartree;
            const std::size_t loffset_X_b = std::size_t(istate) * std::size_t(this->ldim);
            for (int id = 0; id < 3; ++id)
            {
                std::complex<double> td = 0.0; // short name of transition dipole
                for (int is = 0; is < this->nspin_x; ++is)
                {
                    const std::size_t loffset_X_bs = loffset_X_b + is * nk * pX[0].get_local_size();
                    const int goffset_v_ds = (is * 3 + id) * nk * nbands * nbands;
                    for (int ik = 0; ik < nk; ++ik)
                    {
                        const double wk = this->kv.wk[ik] * nk; // k-point weight, normalized as sum = nk
                        const std::size_t loffset_X = loffset_X_bs + ik * pX[is].get_local_size();
                        const int goffset_v = goffset_v_ds + ik * nbands * nbands;
                        for (int io = 0; io < pX[is].get_col_size(); ++io)    // nocc_local
                        {
                            int io_g = pX[is].local2global_col(io);
                            for (int iv = 0; iv < pX[is].get_row_size(); ++iv)    // nvirt_local
                            {
                                int iv_g = pX[is].local2global_row(iv);
                                const std::size_t X_index = loffset_X + io * pX[is].get_row_size() + iv;
                                const int v_index = goffset_v + (iv_g+nocc[is]) * nbands + io_g;
                                if (use_ks_gap)
                                {
                                    td += this->vmo_ptr[v_index] * X[X_index] / eig_ks_diff[X_index - loffset_X_b] * wk;
                                    if (this->is_full)
                                    {
                                        // THE HERMITIAN CONJUGATE OF VMO HAS BEEN VERIFIED
                                        // const int v_index2 = goffset_v + io_g * nbands + iv_g + nocc[is];
                                        // if (std::abs(vmo_ptr[v_index2] - std::conj(this->vmo_ptr[v_index])) > 1e-5){
                                        //     std::cout<<"io:"<<io_g<<" iv:"<<iv_g<<" v:"<<vmo_ptr[v_index]
                                        //         <<" v^T:" << vmo_ptr[v_index2]<<std::endl;
                                        //     }
                                        // <ik|v|ak>X_{aik}/(Ea-Ei) + <ak|v|ik>Y_{aik}/(Ei-Ea)
                                        td -= std::conj(this->vmo_ptr[v_index]) * Y[X_index] / eig_ks_diff[X_index - loffset_X_b] * wk;
                                    }
                                }
                                else
                                {
                                    td += this->vmo_ptr[v_index] * X[X_index] * wk;
                                    if (this->is_full)
                                    {
                                        td += std::conj(this->vmo_ptr[v_index]) * Y[X_index] * wk;
                                    }
                                }
                            }
                        }
                    }
                }   // end for spin_x, only matter in open-shell system
                td *= fac;
                if (this->nspin_x == 1) { td *= sqrt(2.0); } // *2 for 2 spins, /sqrt(2) for the halfed dimension of X in the normalizaiton
                trans_dipole_buf[3 * istate + id] = td;
            }   // end for direction
        }
        Parallel_Reduce::reduce_all(trans_dipole_buf.data(), 3 * nstate);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int istate = 0; istate < nstate; ++istate)
        {
            std::complex<double>* ptr = &trans_dipole_buf[3 * istate];
            this->transition_dipole_[istate] = convert_ptr_to_vector3<T>(ptr);
            this->mean_squared_transition_dipole_[istate] = cal_mean_squared_dipole(transition_dipole_[istate]);
        }
        ModuleBase::timer::end("LR_Spectrum", "cal_transition_dipole_istate_velocity_mo");
    }

    template<typename T>
    void LR::LR_Spectrum<T>::test_transition_dipoles_velocity_omega()
    {
        ModuleBase::timer::start("LR_Spectrum", "test_transition_dipoles_velocity_omega");
        this->transition_dipole_.resize(nstate);
        this->mean_squared_transition_dipole_.resize(nstate);
        this->cal_transition_dipole_istate_velocity_mo(DipoleEnergyType::LR_EIG, {});
        this->oscillator_strength();
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LR::LR_Spectrum::test_transition_dipoles_velocity_omega");
        ModuleBase::timer::end("LR_Spectrum", "test_transition_dipoles_velocity_omega");
    }

    inline void cal_eig_ks_diff(double* const eig_ks_diff, const double* const eig_ks, const Parallel_2D& px, const int nk, const int nocc, const int nvirt)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ik = 0;ik < nk;++ik)
        {
            const int start_k = ik * (nocc + nvirt);
            for (int io = 0;io < px.get_col_size();++io)    //nocc_local
            {
                for (int iv = 0;iv < px.get_row_size();++iv)    //nvirt_local
                {
                    int io_g = px.local2global_col(io);
                    int iv_g = px.local2global_row(iv);
                    eig_ks_diff[ik * px.get_local_size() + io * px.get_row_size() + iv] = (eig_ks[start_k + nocc + iv_g] - eig_ks[start_k + io_g]) / ModuleBase::e2;  // Ry to Hartree
                }
            }
        }
    }

    template<typename T>
    void LR::LR_Spectrum<T>::cal_transition_dipoles_velocity(const double* const eig_ks)
    {
        ModuleBase::timer::start("LR_Spectrum", "cal_transition_dipoles_velocity");

        this->transition_dipole_.resize(nstate);
        this->mean_squared_transition_dipole_.resize(nstate);

        //  (e_c-e_v) of KS eigenvalues
        std::vector<double> eig_ks_diff(this->ldim);
        for (int is = 0;is < this->nspin_x;++is)
        {
            cal_eig_ks_diff(eig_ks_diff.data() + is * nk * pX[0].get_local_size(), eig_ks, pX[is], nk, nocc[is], nvirt[is]);
        }

        //  X/(ec-ev)
        // std::vector<T> X_div_ks_eig(nstate * this->ldim);
        // for (int istate = 0;istate < nstate;++istate)
        // {
        //     const int st = istate * this->ldim;
        //     std::transform(X + st, X + st + ldim, eig_ks_diff.begin(), X_div_ks_eig.data() + st, std::divides<T>());
        // }
        
        this->cal_transition_dipole_istate_velocity_mo(DipoleEnergyType::KS_GAP, eig_ks_diff);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LR::LR_Spectrum::cal_transition_dipoles_velocity");
        ModuleBase::timer::end("LR_Spectrum", "cal_transition_dipoles_velocity");
    }

} // namespace LR

template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;
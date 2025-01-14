#ifdef __DEEPKS
#include "LCAO_deepks_interface.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/cal_dm.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_parameter/parameter.h"

template <typename TK, typename TR>
LCAO_Deepks_Interface<TK, TR>::LCAO_Deepks_Interface(std::shared_ptr<LCAO_Deepks> ld_in) : ld(ld_in)
{
}

template <typename TK, typename TR>
void LCAO_Deepks_Interface<TK, TR>::out_deepks_labels(const double& etot,
                                                      const int& nks,
                                                      const int& nat,
                                                      const int& nlocal,
                                                      const ModuleBase::matrix& ekb,
                                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                      const UnitCell& ucell,
                                                      const LCAO_Orbitals& orb,
                                                      const Grid_Driver& GridD,
                                                      const Parallel_Orbitals* ParaV,
                                                      const psi::Psi<TK>& psi,
                                                      const elecstate::DensityMatrix<TK, double>* dm,
                                                      hamilt::HamiltLCAO<TK, TR>* p_ham)
{
    ModuleBase::TITLE("LCAO_Deepks_Interface", "out_deepks_labels");
    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");

    // Note: out_deepks_labels does not support equivariant version now!

    // define TH for different types
    using TH = std::conditional_t<std::is_same<TK, double>::value, ModuleBase::matrix, ModuleBase::ComplexMatrix>;

    // These variables are frequently used in the following code
    const int inlmax = orb.Alpha[0].getTotal_nchi() * nat;
    const int lmaxd = orb.get_lmax_d();
    const int nmaxd = ld->nmaxd;

    const int des_per_atom = ld->des_per_atom;
    const int* inl_l = ld->inl_l;
    const ModuleBase::IntArray* inl_index = ld->inl_index;
    const std::vector<hamilt::HContainer<double>*> phialpha = ld->phialpha;

    std::vector<torch::Tensor> pdm = ld->pdm;
    bool init_pdm = ld->init_pdm;
    double E_delta = ld->E_delta;
    double e_delta_band = ld->e_delta_band;

    const int my_rank = GlobalV::MY_RANK;
    const int nspin = PARAM.inp.nspin;

    // Note : update PDM and all other quantities with the current dm
    // DeePKS PDM and descriptor
    if (PARAM.inp.deepks_out_labels || PARAM.inp.deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        DeePKS_domain::cal_pdm<TK>
                (init_pdm, inlmax, lmaxd, inl_l, inl_index, dm, phialpha, ucell, orb, GridD, *ParaV, pdm);

        DeePKS_domain::check_pdm(inlmax, inl_l, pdm); // print out the projected dm for NSCF calculaiton

        std::vector<torch::Tensor> descriptor;
        DeePKS_domain::cal_descriptor(nat, inlmax, inl_l, pdm, descriptor,
                                      des_per_atom); // final descriptor
        DeePKS_domain::check_descriptor(inlmax, des_per_atom, inl_l, ucell, PARAM.globalv.global_out_dir, descriptor);

        if (PARAM.inp.deepks_out_labels)
        {
            LCAO_deepks_io::save_npy_d(nat,
                                       des_per_atom,
                                       inlmax,
                                       inl_l,
                                       PARAM.inp.deepks_equiv,
                                       descriptor,
                                       PARAM.globalv.global_out_dir,
                                       GlobalV::MY_RANK); // libnpy needed
        }

        if (PARAM.inp.deepks_scf)
        {
            // update E_delta and gedm
            // new gedm is also useful in cal_f_delta, so it should be ld->gedm
            DeePKS_domain::cal_edelta_gedm(nat,
                        lmaxd,
                        nmaxd,
                        inlmax,
                        des_per_atom,
                        inl_l,
                        descriptor,
                        pdm,
                        ld->model_deepks,
                        ld->gedm,
                        E_delta);
        }
    }

    // Used for deepks_bandgap == 1 and deepks_v_delta > 0
    std::vector<std::vector<TK>>* h_delta = nullptr;
    if constexpr (std::is_same<TK, double>::value)
    {
        h_delta = &ld->H_V_delta;
    }
    else
    {
        h_delta = &ld->H_V_delta_k;
    }

    // calculating deepks correction and save the results
    if (PARAM.inp.deepks_out_labels)
    {
        // Used for deepks_scf == 1
        std::vector<torch::Tensor> gevdm;
        if (PARAM.inp.deepks_scf)
        {
            DeePKS_domain::cal_gevdm(nat, inlmax, inl_l, pdm, gevdm);
        }

        // Energy Part
        const std::string file_etot = PARAM.globalv.global_out_dir + "deepks_etot.npy";
        const std::string file_ebase = PARAM.globalv.global_out_dir + "deepks_ebase.npy";

        LCAO_deepks_io::save_npy_e(etot, file_etot, my_rank);

        if (PARAM.inp.deepks_scf)
        {
            /// ebase :no deepks E_delta including
            LCAO_deepks_io::save_npy_e(etot - E_delta, file_ebase, my_rank);
        }
        else // deepks_scf = 0; base calculation
        {
            /// no scf, e_tot=e_base
            LCAO_deepks_io::save_npy_e(etot, file_ebase, my_rank);
        }

        // Force Part
        if (PARAM.inp.cal_force)
        {
            if (PARAM.inp.deepks_scf
                && !PARAM.inp.deepks_equiv) // training with force label not supported by equivariant version now
            {
                std::vector<std::vector<TK>> dm_vec = dm->get_DMK_vector();
                torch::Tensor gdmx;
                DeePKS_domain::cal_gdmx<
                    TK>(lmaxd, inlmax, nks, kvec_d, phialpha, inl_index, dm_vec, ucell, orb, *ParaV, GridD, gdmx);

                torch::Tensor gvx;
                DeePKS_domain::cal_gvx(ucell.nat, inlmax, des_per_atom, inl_l, gevdm, gdmx, gvx);
                const std::string file_gradvx = PARAM.globalv.global_out_dir + "deepks_gradvx.npy";
                LCAO_deepks_io::save_tensor2npy<double>(file_gradvx, gvx, my_rank);

                if (PARAM.inp.deepks_out_unittest)
                {
                    DeePKS_domain::check_gdmx(gdmx);
                    DeePKS_domain::check_gvx(gvx);
                }
            }
        }

        // Stress Part
        if (PARAM.inp.cal_stress)
        {
            if (PARAM.inp.deepks_scf
                && !PARAM.inp.deepks_equiv) // training with stress label not supported by equivariant version now
            {
                std::vector<std::vector<TK>> dm_vec = dm->get_DMK_vector();
                torch::Tensor gdmepsl;
                DeePKS_domain::cal_gdmepsl<
                    TK>(lmaxd, inlmax, nks, kvec_d, phialpha, inl_index, dm_vec, ucell, orb, *ParaV, GridD, gdmepsl);

                torch::Tensor gvepsl;
                DeePKS_domain::cal_gvepsl(ucell.nat, inlmax, des_per_atom, inl_l, gevdm, gdmepsl, gvepsl);
                const std::string file_gvepsl = PARAM.globalv.global_out_dir + "deepks_gvepsl.npy";
                LCAO_deepks_io::save_tensor2npy<double>(file_gvepsl, gvepsl, my_rank);

                if (PARAM.inp.deepks_out_unittest)
                {
                    DeePKS_domain::check_gdmepsl(gdmepsl);
                    DeePKS_domain::check_gvepsl(gvepsl);
                }
            }
        }

        // Bandgap Part
        if (PARAM.inp.deepks_bandgap)
        {
            const int nocc = (PARAM.inp.nelec + 1) / 2;
            ModuleBase::matrix o_tot(nks, 1);
            for (int iks = 0; iks < nks; ++iks)
            {
                // record band gap for each k point (including spin)
                o_tot(iks, 0) = ekb(iks, nocc) - ekb(iks, nocc - 1);
            }

            const std::string file_otot = PARAM.globalv.global_out_dir + "deepks_otot.npy";
            LCAO_deepks_io::save_matrix2npy(file_otot, o_tot, my_rank); // Unit: Ry

            if (PARAM.inp.deepks_scf)
            {
                ModuleBase::matrix wg_hl;
                std::vector<TH> dm_bandgap;

                // Calculate O_delta
                if constexpr (std::is_same<TK, double>::value) // for gamma only
                {
                    wg_hl.create(nspin, PARAM.inp.nbands);
                    dm_bandgap.resize(nspin);
                    for (int is = 0; is < nspin; ++is)
                    {
                        wg_hl.zero_out();
                        wg_hl(is, nocc - 1) = -1.0;
                        wg_hl(is, nocc) = 1.0;
                        elecstate::cal_dm(ParaV, wg_hl, psi, dm_bandgap);
                    }
                }
                else // for multi-k
                {
                    wg_hl.create(nks, PARAM.inp.nbands);
                    dm_bandgap.resize(nks);
                    wg_hl.zero_out();
                    for (int ik = 0; ik < nks; ik++)
                    {
                        wg_hl(ik, nocc - 1) = -1.0;
                        wg_hl(ik, nocc) = 1.0;
                    }
                    elecstate::cal_dm(ParaV, wg_hl, psi, dm_bandgap);
                }

                ModuleBase::matrix o_delta(nks, 1);

                // calculate and save orbital_precalc: [nks,NAt,NDscrpt]
                torch::Tensor orbital_precalc;
                DeePKS_domain::cal_orbital_precalc<TK, TH>(dm_bandgap,
                                                           lmaxd,
                                                           inlmax,
                                                           nat,
                                                           nks,
                                                           inl_l,
                                                           kvec_d,
                                                           phialpha,
                                                           gevdm,
                                                           inl_index,
                                                           ucell,
                                                           orb,
                                                           *ParaV,
                                                           GridD,
                                                           orbital_precalc);
                DeePKS_domain::cal_o_delta<TK, TH>(dm_bandgap, *h_delta, o_delta, *ParaV, nks);

                // save obase and orbital_precalc
                const std::string file_orbpre = PARAM.globalv.global_out_dir + "deepks_orbpre.npy";
                LCAO_deepks_io::save_tensor2npy<double>(file_orbpre, orbital_precalc, my_rank);

                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_matrix2npy(file_obase, o_tot - o_delta, my_rank); // Unit: Ry
            }                                                                          // end deepks_scf == 1
            else                                                                       // deepks_scf == 0
            {
                const std::string file_obase = PARAM.globalv.global_out_dir + "deepks_obase.npy";
                LCAO_deepks_io::save_matrix2npy(file_obase, o_tot, my_rank); // no scf, o_tot=o_base
            }                                                                // end deepks_scf == 0
        }                                                                    // end bandgap label

        // H(R) matrix part, not realized now
        if (true) // should be modified later!
        {
            const std::string file_hr = PARAM.globalv.global_out_dir + "deepks_hr.npy";
            const hamilt::HContainer<TR>& hR = *(p_ham->getHR());

            // How to save H(R)?
        }

        // H(k) matrix part
        if (PARAM.inp.deepks_v_delta)
        {
            std::vector<TH> h_tot(nks);
            std::vector<std::vector<TK>> h_mat(nks, std::vector<TK>(ParaV->nloc));
            for (int ik = 0; ik < nks; ik++)
            {
                h_tot[ik].create(nlocal, nlocal);
                p_ham->updateHk(ik);
                const TK* hk_ptr = p_ham->getHk();
                for (int i = 0; i < ParaV->nloc; i++)
                {
                    h_mat[ik][i] = hk_ptr[i];
                }
            }

            DeePKS_domain::collect_h_mat<TK, TH>(*ParaV, h_mat, h_tot, nlocal, nks);

            const std::string file_htot = PARAM.globalv.global_out_dir + "deepks_htot.npy";
            LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_htot, nlocal, nks, my_rank);

            if (PARAM.inp.deepks_scf)
            {
                std::vector<TH> v_delta(nks);
                std::vector<TH> h_base(nks);
                for (int ik = 0; ik < nks; ik++)
                {
                    v_delta[ik].create(nlocal, nlocal);
                    h_base[ik].create(nlocal, nlocal);
                }
                DeePKS_domain::collect_h_mat<TK, TH>(*ParaV, *h_delta, v_delta, nlocal, nks);

                // save v_delta and h_base
                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                for (int ik = 0; ik < nks; ik++)
                {
                    h_base[ik] = h_tot[ik] - v_delta[ik];
                }
                LCAO_deepks_io::save_npy_h<TK, TH>(h_base, file_hbase, nlocal, nks, my_rank);

                const std::string file_vdelta = PARAM.globalv.global_out_dir + "deepks_vdelta.npy";
                LCAO_deepks_io::save_npy_h<TK, TH>(v_delta, file_vdelta, nlocal, nks, my_rank);

                if (PARAM.inp.deepks_v_delta == 1) // v_delta_precalc storage method 1
                {
                    torch::Tensor v_delta_precalc;
                    DeePKS_domain::cal_v_delta_precalc<TK>(nlocal,
                                                           lmaxd,
                                                           inlmax,
                                                           nat,
                                                           nks,
                                                           inl_l,
                                                           kvec_d,
                                                           phialpha,
                                                           gevdm,
                                                           inl_index,
                                                           ucell,
                                                           orb,
                                                           *ParaV,
                                                           GridD,
                                                           v_delta_precalc);

                    const std::string file_vdpre = PARAM.globalv.global_out_dir + "deepks_vdpre.npy";
                    LCAO_deepks_io::save_tensor2npy<TK>(file_vdpre, v_delta_precalc, my_rank);
                }
                else if (PARAM.inp.deepks_v_delta == 2) // v_delta_precalc storage method 2
                {
                    torch::Tensor phialpha_out;
                    DeePKS_domain::prepare_phialpha<
                        TK>(nlocal, lmaxd, inlmax, nat, nks, kvec_d, phialpha, ucell, orb, *ParaV, GridD, phialpha_out);
                    const std::string file_phialpha = PARAM.globalv.global_out_dir + "deepks_phialpha.npy";
                    LCAO_deepks_io::save_tensor2npy<TK>(file_phialpha, phialpha_out, my_rank);

                    torch::Tensor gevdm_out;
                    DeePKS_domain::prepare_gevdm(nat, lmaxd, inlmax, orb, gevdm, gevdm_out);
                    const std::string file_gevdm = PARAM.globalv.global_out_dir + "deepks_gevdm.npy";
                    LCAO_deepks_io::save_tensor2npy<double>(file_gevdm, gevdm_out, my_rank);
                }
            }
            else // deepks_scf == 0
            {
                const std::string file_hbase = PARAM.globalv.global_out_dir + "deepks_hbase.npy";
                LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_hbase, nlocal, nks, my_rank);
            }
        } // end v_delta label

    } // end deepks_out_labels

    /// print out deepks information to the screen
    if (PARAM.inp.deepks_scf)
    {
        DeePKS_domain::cal_e_delta_band(dm->get_DMK_vector(), *h_delta, nks, ParaV, e_delta_band);
        std::cout << "E_delta_band = " << std::setprecision(8) << e_delta_band << " Ry"
                  << " = " << std::setprecision(8) << e_delta_band * ModuleBase::Ry_to_eV << " eV" << std::endl;
        std::cout << "E_delta_NN = " << std::setprecision(8) << E_delta << " Ry"
                  << " = " << std::setprecision(8) << E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
        if (PARAM.inp.deepks_out_unittest)
        {
            LCAO_deepks_io::print_dm(nks, PARAM.globalv.nlocal, ParaV->nrow, dm->get_DMK_vector());

            DeePKS_domain::check_gedm(inlmax, inl_l, ld->gedm);

            std::ofstream ofs("E_delta_bands.dat");
            ofs << std::setprecision(10) << e_delta_band;

            std::ofstream ofs1("E_delta.dat");
            ofs1 << std::setprecision(10) << E_delta;
        }
    }
    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");
}

template class LCAO_Deepks_Interface<double, double>;
template class LCAO_Deepks_Interface<std::complex<double>, double>;
template class LCAO_Deepks_Interface<std::complex<double>, std::complex<double>>;

#endif

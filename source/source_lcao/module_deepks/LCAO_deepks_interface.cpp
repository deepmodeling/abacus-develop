#ifdef __MLALGO
#include "LCAO_deepks_interface.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "source_estate/cal_dm.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_variable.h"
#include "source_base/tool_title.h"

#include <unordered_map>

template <typename TK, typename TR>
LCAO_Deepks_Interface<TK, TR>::LCAO_Deepks_Interface(std::shared_ptr<LCAO_Deepks<TK>> ld_in) : ld(ld_in)
{
}

// Helper function to map file_type to true names
std::string true_file_type(const std::string& file_type) 
{
    static const std::unordered_map<std::string, std::string> file_type_map = {
        {"etot", "energy"},
        {"ftot", "force"},
        {"stot", "stress"},
        {"otot", "orbital"},
        {"htot", "hamiltonian"}
    };

    auto it = file_type_map.find(file_type);
    return it != file_type_map.end() ? it->second : file_type;
}

// global_out_dir/deepks_*.npy for iter=-1 (called in after_scf)
// global_out_dir/DeePKS_Labels_Elec/*_e*.npy for iter>0 (called during electronic steps)
template <typename TK, typename TR>
std::string LCAO_Deepks_Interface<TK, TR>::get_filename(const std::string& file_type,
                                                        const int& label_type,
                                                        const int& iter) 
{
    std::ostringstream file_name;
    file_name << (iter == -1 ? PARAM.globalv.global_out_dir : PARAM.globalv.global_deepks_label_elec_dir);
    if (iter == -1) 
    {
        file_name << "deepks_";
    }
    file_name << (label_type == 2 ? true_file_type(file_type) : file_type);
    if (iter != -1)
    {
        file_name << "_e" << iter;
    }
    file_name << ".npy";
    return file_name.str();
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
                                                      hamilt::HamiltLCAO<TK, TR>* p_ham,
                                                      const int& iter,
                                                      const bool& conv_esolver,
                                                      const int rank,
                                                      std::ostream& ofs_running)
{
    ModuleBase::TITLE("LCAO_Deepks_Interface", "out_deepks_labels");
    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");

    // Note: out_deepks_labels does not support equivariant version now!

    // define TH for different types
    using TH = std::conditional_t<std::is_same<TK, double>::value, ModuleBase::matrix, ModuleBase::ComplexMatrix>;

    // These variables are frequently used in the following code
    const int nlmax = orb.Alpha[0].getTotal_nchi();
    const int inlmax = nlmax * nat;
    const int lmaxd = orb.get_lmax_d();
    const int nmaxd = ld->nmaxd;

    const int des_per_atom = ld->des_per_atom;
    const std::vector<int> inl2l = ld->inl2l;
    const ModuleBase::IntArray* inl_index = ld->inl_index;
    const std::vector<hamilt::HContainer<double>*> phialpha = ld->phialpha;

    std::vector<torch::Tensor> pdm = ld->pdm;
    bool init_pdm = ld->init_pdm;
    double E_delta = ld->E_delta;
    double e_delta_band = ld->e_delta_band;
    hamilt::HContainer<double>* dmr = ld->dm_r;
    // Used for deepks_bandgap == 1 and deepks_v_delta > 0
    std::vector<std::vector<TK>>* h_delta = &ld->V_delta;

    const int nspin = PARAM.inp.nspin;
    const int nk = nks / nspin;

    const bool is_after_scf = (iter == -1); // called in after_scf, not in electronic steps
    const bool output_base = PARAM.inp.deepks_out_labels == 1 && is_after_scf; // not output when deepks_out_labels=2 and in electronic step (output true base elsewhere)
    const bool output_precalc = (PARAM.inp.deepks_out_labels == 1) && (PARAM.inp.deepks_scf || PARAM.inp.deepks_out_freq_elec);

//================================================================================
// 1. Update real-space density matrix (DMR) for deepks, projected density matrix (PDM)
//    and descriptor. Output descriptor if needed.
//================================================================================

    // Update DMR in any case of deepks_out_labels/deepks_scf
    DeePKS_domain::update_dmr(kvec_d, dm->get_DMK_vector(), ucell, orb, *ParaV, GridD, dmr);

    // Note : update PDM and all other quantities with the current dm
    // DeePKS PDM and descriptor
    if (PARAM.inp.deepks_out_labels == 1 || PARAM.inp.deepks_scf)
    {
        // this part is for integrated test of deepks
        // so it is printed no matter even if deepks_out_labels is not used
        DeePKS_domain::cal_pdm<
            TK>(init_pdm, inlmax, lmaxd, inl2l, inl_index, kvec_d, dmr, phialpha, ucell, orb, GridD, *ParaV, pdm);

        DeePKS_domain::check_pdm(inlmax, inl2l, pdm); // print out the projected dm for NSCF calculaiton

        std::vector<torch::Tensor> descriptor;
        DeePKS_domain::cal_descriptor(nat, inlmax, inl2l, pdm, descriptor,
                                      des_per_atom); // final descriptor
        DeePKS_domain::check_descriptor(inlmax,
                                        des_per_atom,
                                        inl2l,
                                        ucell,
                                        PARAM.globalv.global_out_dir,
                                        descriptor,
                                        rank);
        
        const std::string file_d = get_filename("dm_eig", PARAM.inp.deepks_out_labels, iter);
        LCAO_deepks_io::save_npy_d(nat,
                                des_per_atom,
                                inlmax,
                                inl2l,
                                PARAM.inp.deepks_equiv,
                                descriptor,
                                file_d,
                                rank); // libnpy needed


        if (PARAM.inp.deepks_scf)
        {
            // update E_delta and gedm
            // new gedm is also useful in cal_f_delta, so it should be ld->gedm
            if (PARAM.inp.deepks_equiv)
            {
                DeePKS_domain::cal_edelta_gedm_equiv(nat,
                                                     lmaxd,
                                                     nmaxd,
                                                     inlmax,
                                                     des_per_atom,
                                                     inl2l,
                                                     descriptor,
                                                     ld->gedm,
                                                     E_delta,
                                                     rank);
            }
            else
            {
                DeePKS_domain::cal_edelta_gedm(nat,
                                               inlmax,
                                               des_per_atom,
                                               inl2l,
                                               descriptor,
                                               pdm,
                                               ld->model_deepks,
                                               ld->gedm,
                                               E_delta);
            }
        }
    }

    // calculating deepks correction and save the results
    if (PARAM.inp.deepks_out_labels)
    {
        // Used for deepks_scf == 1 or deepks_out_freq_elec!=0, for *precalc items, not for deepks_out_labels=2
        std::vector<torch::Tensor> gevdm;
        if ( output_precalc )
        {
            DeePKS_domain::cal_gevdm(nat, inlmax, inl2l, pdm, gevdm);
        }

//================================================================================
// 2. Energy
//================================================================================

        // etot
        const std::string file_etot = get_filename("etot", PARAM.inp.deepks_out_labels, iter);
        LCAO_deepks_io::save_npy_e(etot, file_etot, rank);

        // ebase
        if ( output_base )
        {
            const std::string file_ebase = get_filename("ebase", PARAM.inp.deepks_out_labels, iter);
            if (PARAM.inp.deepks_scf)
            {
                /// ebase :no deepks E_delta including
                LCAO_deepks_io::save_npy_e(etot - E_delta, file_ebase, rank);
            }
            else // deepks_scf = 0; base calculation
            {
                /// no scf, e_tot=e_base
                LCAO_deepks_io::save_npy_e(etot, file_ebase, rank);
            }
        }

//================================================================================
// 3. Force and Stress
//================================================================================

        if ( is_after_scf )
        {
            // Force Part
            if (PARAM.inp.cal_force)
            {
                // these items are not related to model, so can output without deepks_scf
                if (output_precalc // don't need these when deepks_out_labels == 2
                    && !PARAM.inp.deepks_equiv) // training with force label not supported by equivariant version now
                {
                    torch::Tensor gdmx;
                    DeePKS_domain::cal_gdmx<
                        TK>(lmaxd, inlmax, nks, kvec_d, phialpha, inl_index, dmr, ucell, orb, *ParaV, GridD, gdmx);

                    torch::Tensor gvx;
                    DeePKS_domain::cal_gvx(ucell.nat, inlmax, des_per_atom, inl2l, gevdm, gdmx, gvx, rank);
                    const std::string file_gradvx = get_filename("gradvx", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_tensor2npy<double>(file_gradvx, gvx, rank);

                    if (PARAM.inp.deepks_out_unittest)
                    {
                        DeePKS_domain::check_tensor<double>(gdmx, "gdmx.dat", rank);
                        DeePKS_domain::check_tensor<double>(gvx, "gvx.dat", rank);
                    }
                }
            }

            // Stress Part
            if (PARAM.inp.cal_stress)
            {
                // these items are not related to model, so can output without deepks_scf
                if (output_precalc // don't need these when deepks_out_labels == 2
                    && !PARAM.inp.deepks_equiv) // training with stress label not supported by equivariant version now
                {
                    torch::Tensor gdmepsl;
                    DeePKS_domain::cal_gdmepsl<
                        TK>(lmaxd, inlmax, nks, kvec_d, phialpha, inl_index, dmr, ucell, orb, *ParaV, GridD, gdmepsl);

                    torch::Tensor gvepsl;
                    DeePKS_domain::cal_gvepsl(ucell.nat, inlmax, des_per_atom, inl2l, gevdm, gdmepsl, gvepsl, rank);
                    const std::string file_gvepsl = get_filename("gvepsl", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_tensor2npy<double>(file_gvepsl, gvepsl, rank);

                    if (PARAM.inp.deepks_out_unittest)
                    {
                        DeePKS_domain::check_tensor<double>(gdmepsl, "gdmepsl.dat", rank);
                        DeePKS_domain::check_tensor<double>(gvepsl, "gvepsl.dat", rank);
                    }
                }
            }
        }

//================================================================================
// 4. Bandgap
//================================================================================

        if (PARAM.inp.deepks_bandgap > 0)
        {
            // Get the number of the occupied bands
            // Notice that the index of band starts from 0, so actually (nocc - 1) is the index of HOMO state
            int nocc = (PARAM.inp.nelec + 1) / 2;
            if (PARAM.inp.deepks_bandgap == 3)
            {
                int natom_H = 0;
                for (int it = 0; it < ucell.ntype; it++)
                {
                    if (ucell.atoms[it].label == "H")
                    {
                        natom_H = ucell.atoms[it].na;
                        break;
                    }
                }
                nocc = (PARAM.inp.nelec - natom_H) / 2;
            }

            // Get the number of bandgap to be recorded
            int range = 1;                     // normally use only one gap
            if (PARAM.inp.deepks_bandgap == 2) // for bandgap label 2, use multi bandgap
            {
                range = PARAM.inp.deepks_band_range[1] - PARAM.inp.deepks_band_range[0] + 1;
                // For cases where HOMO(nocc - 1) is included in the range
                if ((PARAM.inp.deepks_band_range[0] <= -1) && (PARAM.inp.deepks_band_range[1] >= -1))
                {
                    range -= 1;
                }
            }

            // Calculate the bandgap for each k point
            ModuleBase::matrix o_tot(nks, range);
            for (int iks = 0; iks < nks; ++iks)
            {
                int ib = 0;
                if (PARAM.inp.deepks_bandgap == 1 || PARAM.inp.deepks_bandgap == 3)
                {
                    o_tot(iks, ib) = ekb(iks, nocc + PARAM.inp.deepks_band_range[1])
                                    - ekb(iks, nocc + PARAM.inp.deepks_band_range[0]);
                }
                else if (PARAM.inp.deepks_bandgap == 2)
                {
                    for (int ir = PARAM.inp.deepks_band_range[0]; ir <= PARAM.inp.deepks_band_range[1]; ++ir)
                    {
                        if (ir != -1)
                        {
                            o_tot(iks, ib) = ekb(iks, nocc + ir) - ekb(iks, nocc - 1);
                            ib++;
                        }
                    }
                    assert(ib == range); // ensure that we have filled all the bandgap values
                }
            }

            const std::string file_otot = get_filename("otot", PARAM.inp.deepks_out_labels, iter);
            LCAO_deepks_io::save_matrix2npy(file_otot, o_tot, rank); // Unit: Hartree

            // don't need these when deepks_out_labels == 2
            // not consider out_base now, because bandgap is not supported in deepks_out_freq_elec now
            if ( output_precalc ) 
            {
                std::vector<ModuleBase::matrix> wg_hl_range(range);
                for (int ir = 0; ir < range; ++ir)
                {
                    wg_hl_range[ir].create(nks, PARAM.inp.nbands);
                    wg_hl_range[ir].zero_out();
                }

                // Calculate O_delta
                for (int iks = 0; iks < nks; ++iks)
                {
                    int ib = 0;
                    if (PARAM.inp.deepks_bandgap == 1 || PARAM.inp.deepks_bandgap == 3)
                    {
                        wg_hl_range[ib](iks, nocc + PARAM.inp.deepks_band_range[0]) = -1.0;
                        wg_hl_range[ib](iks, nocc + PARAM.inp.deepks_band_range[1]) = 1.0;
                    }
                    else if (PARAM.inp.deepks_bandgap == 2)
                    {
                        for (int ir = PARAM.inp.deepks_band_range[0]; ir <= PARAM.inp.deepks_band_range[1]; ++ir)
                        {
                            if (ir != -1)
                            {
                                wg_hl_range[ib](iks, nocc - 1) = -1.0;
                                wg_hl_range[ib](iks, nocc + ir) = 1.0;
                                ib++;
                            }
                        }
                    }
                }

                ModuleBase::matrix o_delta(nks, range);
                torch::Tensor orbital_precalc;
                for (int ir = 0; ir < range; ++ir)
                {
                    std::vector<TH> dm_bandgap(nks);
                    elecstate::cal_dm(ParaV, wg_hl_range[ir], psi, dm_bandgap);

                    torch::Tensor orbital_precalc_temp;
                    ModuleBase::matrix o_delta_temp(nks, 1);
                    DeePKS_domain::cal_orbital_precalc<TK, TH>(dm_bandgap,
                                                                lmaxd,
                                                                inlmax,
                                                                nat,
                                                                nks,
                                                                inl2l,
                                                                kvec_d,
                                                                phialpha,
                                                                gevdm,
                                                                inl_index,
                                                                ucell,
                                                                orb,
                                                                *ParaV,
                                                                GridD,
                                                                orbital_precalc_temp);
                    if (ir == 0)
                    {
                        orbital_precalc = orbital_precalc_temp;
                    }
                    else
                    {
                        orbital_precalc = torch::cat({orbital_precalc, orbital_precalc_temp}, 0);
                    }

                    if (PARAM.inp.deepks_scf)
                    {
                        DeePKS_domain::cal_o_delta<TK, TH>(dm_bandgap, *h_delta, o_delta_temp, *ParaV, nks, nspin);
                        for (int iks = 0; iks < nks; ++iks)
                        {
                            o_delta(iks, ir) = o_delta_temp(iks, 0);
                        }                            
                    }
                }
                // save obase and orbital_precalc
                const std::string file_orbpre = get_filename("orbpre", PARAM.inp.deepks_out_labels, iter);
                LCAO_deepks_io::save_tensor2npy<double>(file_orbpre, orbital_precalc, rank);                        

                if (PARAM.inp.deepks_scf)
                {
                    const std::string file_obase = get_filename("obase", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_matrix2npy(file_obase, o_tot - o_delta, rank); // Unit: Hartree                            
                }
                else
                {
                    const std::string file_obase = get_filename("obase", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_matrix2npy(file_obase, o_tot, rank); // no scf, o_tot=o_base
                }                          
            }                                                                 // end out_precalc
        }                                                                     // end deepks_bandgap > 0

//================================================================================
// 5. HR
//================================================================================
        if ( is_after_scf )
        {
            // not add deepks_out_labels = 2 and deepks_out_freq_elec for HR yet
            // H(R) matrix part, for HR, base will not be calculated since they are HContainer objects
            if (PARAM.inp.deepks_v_delta < 0)
            {
                // set the output
                const double sparse_threshold = 1e-10;
                const int precision = 8;
                const std::string file_hrtot
                    = PARAM.globalv.global_out_dir
                    + (PARAM.inp.deepks_out_labels == 1 ? "deepks_hrtot.csr" : "deepks_hamiltonian_r.csr");
                hamilt::HContainer<TR>* hR_tot = (p_ham->getHR());

                if (rank == 0)
                {
                    std::ofstream ofs_hr(file_hrtot, std::ios::out);
                    ofs_hr << "Matrix Dimension of H(R): " << hR_tot->get_nbasis() << std::endl;
                    ofs_hr << "Matrix number of H(R): " << hR_tot->size_R_loop() << std::endl;
                    hamilt::Output_HContainer<TR> out_hr(hR_tot, ofs_hr, sparse_threshold, precision);
                    out_hr.write(true); // write all the matrices, including empty ones
                    ofs_hr.close();
                }

                if (PARAM.inp.deepks_scf)
                {
                    if (PARAM.inp.deepks_out_labels == 1)
                    {
                        const std::string file_vdeltar = PARAM.globalv.global_out_dir + "deepks_hrdelta.csr";
                        hamilt::HContainer<TR>* h_deltaR = p_ham->get_V_delta_R();

                        if (rank == 0)
                        {
                            std::ofstream ofs_hr(file_vdeltar, std::ios::out);
                            ofs_hr << "Matrix Dimension of H_delta(R): " << h_deltaR->get_nbasis() << std::endl;
                            ofs_hr << "Matrix number of H_delta(R): " << h_deltaR->size_R_loop() << std::endl;
                            hamilt::Output_HContainer<TR> out_hr(h_deltaR, ofs_hr, sparse_threshold, precision);
                            out_hr.write(true); // write all the matrices, including empty ones
                            ofs_hr.close();
                        }

                        if (PARAM.inp.deepks_v_delta == -1)
                        {
                            int R_size = DeePKS_domain::get_R_size(*h_deltaR);
                            torch::Tensor vdr_precalc;
                            DeePKS_domain::cal_vdr_precalc(nlocal,
                                                        lmaxd,
                                                        inlmax,
                                                        nat,
                                                        nks,
                                                        R_size,
                                                        inl2l,
                                                        kvec_d,
                                                        phialpha,
                                                        gevdm,
                                                        inl_index,
                                                        ucell,
                                                        orb,
                                                        *ParaV,
                                                        GridD,
                                                        vdr_precalc);

                            const std::string file_vdrpre = PARAM.globalv.global_out_dir + "deepks_vdrpre.npy";
                            LCAO_deepks_io::save_tensor2npy<double>(file_vdrpre, vdr_precalc, rank);
                        }
                        else if (PARAM.inp.deepks_v_delta == -2)
                        {
                            int R_size = DeePKS_domain::get_R_size(*h_deltaR);
                            torch::Tensor phialpha_r_out;
                            DeePKS_domain::prepare_phialpha_r(nlocal,
                                                            lmaxd,
                                                            inlmax,
                                                            nat,
                                                            R_size,
                                                            phialpha,
                                                            ucell,
                                                            orb,
                                                            *ParaV,
                                                            GridD,
                                                            phialpha_r_out);
                            const std::string file_phialpha_r = PARAM.globalv.global_out_dir + "deepks_phialpha_r.npy";
                            LCAO_deepks_io::save_tensor2npy<double>(file_phialpha_r, phialpha_r_out, rank);

                            torch::Tensor gevdm_out;
                            DeePKS_domain::prepare_gevdm(nat, lmaxd, inlmax, orb, gevdm, gevdm_out);
                            const std::string file_gevdm = PARAM.globalv.global_out_dir + "deepks_gevdm.npy";
                            LCAO_deepks_io::save_tensor2npy<double>(file_gevdm, gevdm_out, rank);
                        }
                    }
                }
            }           
        }

//================================================================================
// 6. Hk
//================================================================================

        if (PARAM.inp.deepks_v_delta > 0)
        {
            std::vector<TH> h_tot(nks);
            DeePKS_domain::get_h_tot<TK, TH, TR>(*ParaV, p_ham, h_tot, nlocal, nks, 'H');

            const std::string file_htot = get_filename("htot", PARAM.inp.deepks_out_labels, iter);
            LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_htot, nlocal, nks, rank);

            if ( output_base )
            {
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
                    const std::string file_hbase = get_filename("hbase", PARAM.inp.deepks_out_labels, iter);
                    for (int ik = 0; ik < nks; ik++)
                    {
                        h_base[ik] = h_tot[ik] - v_delta[ik];
                    }
                    LCAO_deepks_io::save_npy_h<TK, TH>(h_base, file_hbase, nlocal, nks, rank);

                    const std::string file_vdelta = get_filename("vdelta", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_npy_h<TK, TH>(v_delta, file_vdelta, nlocal, nks, rank);                            
                }
                else // deepks_scf == 0
                {
                    const std::string file_hbase = get_filename("hbase", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_hbase, nlocal, nks, rank);
                }                        
            }

            if ( output_precalc )
            {
                if (PARAM.inp.deepks_v_delta == 1) // v_delta_precalc storage method 1
                {
                    torch::Tensor v_delta_precalc;
                    DeePKS_domain::cal_v_delta_precalc<TK>(nlocal,
                                                        lmaxd,
                                                        inlmax,
                                                        nat,
                                                        nks,
                                                        inl2l,
                                                        kvec_d,
                                                        phialpha,
                                                        gevdm,
                                                        inl_index,
                                                        ucell,
                                                        orb,
                                                        *ParaV,
                                                        GridD,
                                                        v_delta_precalc);

                    const std::string file_vdpre = get_filename("vdpre", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_tensor2npy<TK>(file_vdpre, v_delta_precalc, rank);
                }
                else if (PARAM.inp.deepks_v_delta == 2) // v_delta_precalc storage method 2
                {
                    torch::Tensor phialpha_out;
                    DeePKS_domain::prepare_phialpha<TK>(nlocal,
                                                        lmaxd,
                                                        inlmax,
                                                        nat,
                                                        nks,
                                                        kvec_d,
                                                        phialpha,
                                                        ucell,
                                                        orb,
                                                        *ParaV,
                                                        GridD,
                                                        phialpha_out);
                    const std::string file_phialpha = get_filename("phialpha", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_tensor2npy<TK>(file_phialpha, phialpha_out, rank);

                    torch::Tensor gevdm_out;
                    DeePKS_domain::prepare_gevdm(nat, lmaxd, inlmax, orb, gevdm, gevdm_out);
                    const std::string file_gevdm = get_filename("gevdm", PARAM.inp.deepks_out_labels, iter);
                    LCAO_deepks_io::save_tensor2npy<double>(file_gevdm, gevdm_out, rank);
                }
            }
        }     // end v_delta label            

    } // end deepks_out_labels

//================================================================================
// 7. atom.npy, box.npy, overlap.npy
//================================================================================

    if ( is_after_scf ) // don't need to output in multiple electronic steps
    {
        if (PARAM.inp.deepks_out_labels == 2)
        {
            // output atom.npy and box.npy
            torch::Tensor atom_out;
            DeePKS_domain::prepare_atom(ucell, atom_out);
            const std::string file_atom = PARAM.globalv.global_out_dir + "deepks_atom.npy";
            LCAO_deepks_io::save_tensor2npy<double>(file_atom, atom_out, rank);

            torch::Tensor box_out;
            DeePKS_domain::prepare_box(ucell, box_out);
            const std::string file_box = PARAM.globalv.global_out_dir + "deepks_box.npy";
            LCAO_deepks_io::save_tensor2npy<double>(file_box, box_out, rank);

            if (PARAM.inp.deepks_v_delta > 0)
            {
                // prepare for overlap.npy, very much like h_tot except for p_ham->getSk()
                std::vector<TH> s_tot(nks);
                DeePKS_domain::get_h_tot<TK, TH, TR>(*ParaV, p_ham, s_tot, nlocal, nks, 'S');
                const std::string file_stot = PARAM.globalv.global_out_dir + "deepks_overlap.npy";
                LCAO_deepks_io::save_npy_h<TK, TH>(s_tot,
                                                file_stot,
                                                nlocal,
                                                nks,
                                                rank,
                                                1.0); // don't need unit_scale for overlap
            }
        }

//================================================================================
// 8. print, unittest
//================================================================================
        /// print out deepks information to the screen
        if (PARAM.inp.deepks_scf)
        {
            DeePKS_domain::cal_e_delta_band(dm->get_DMK_vector(), *h_delta, nks, nspin, ParaV, e_delta_band);
            if (rank == 0)
            {
                ofs_running << " DeePKS Energy Correction" << std::endl;
                ofs_running << " -----------------------------------------------" << std::endl;
                ofs_running << "  E_delta_band = " << std::setprecision(8) << e_delta_band << " Ry"
                            << " = " << std::setprecision(8) << e_delta_band * ModuleBase::Ry_to_eV << " eV" << std::endl;
                ofs_running << "  E_delta_NN = " << std::setprecision(8) << E_delta << " Ry"
                            << " = " << std::setprecision(8) << E_delta * ModuleBase::Ry_to_eV << " eV" << std::endl;
                ofs_running << " -----------------------------------------------" << std::endl;
            }
            if (PARAM.inp.deepks_out_unittest)
            {
                LCAO_deepks_io::print_dm(nks, PARAM.globalv.nlocal, ParaV->nrow, dm->get_DMK_vector());

                DeePKS_domain::check_gedm(inlmax, inl2l, ld->gedm);

                std::ofstream ofs("E_delta_bands.dat");
                ofs << std::setprecision(10) << e_delta_band;

                std::ofstream ofs1("E_delta.dat");
                ofs1 << std::setprecision(10) << E_delta;
            }
        }        
    }   
    ModuleBase::timer::tick("LCAO_Deepks_Interface", "out_deepks_labels");
}

template class LCAO_Deepks_Interface<double, double>;
template class LCAO_Deepks_Interface<std::complex<double>, double>;
template class LCAO_Deepks_Interface<std::complex<double>, std::complex<double>>;

#endif

#include "source_lcao/setup_deepks.h"

#include "source_io/module_parameter/parameter.h" // use parameter
#include "source_lcao/LCAO_domain.h"

template <typename TK>
Setup_DeePKS<TK>::Setup_DeePKS()
{
}

template <typename TK>
Setup_DeePKS<TK>::~Setup_DeePKS()
{
}

template <typename TK>
void Setup_DeePKS<TK>::build_overlap(const UnitCell& ucell,
                                     const LCAO_Orbitals& orb,
                                     const Parallel_Orbitals& pv,
                                     const Grid_Driver& gd,
                                     TwoCenterIntegrator& overlap_orb_alpha,
                                     const Input_para& inp)
{
#ifdef __MLALGO
    // 9) for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions
    if (PARAM.globalv.deepks_setorb)
    {
        // allocate <phi(0)|alpha(R)>, phialpha is different every ion step, so it is allocated here
        DeePKS_domain::allocate_phialpha(inp.cal_force, ucell, orb, gd, &pv, this->ld.phialpha);

        // build and save <phi(0)|alpha(R)> at beginning
        DeePKS_domain::build_phialpha(inp.cal_force, ucell, orb, gd, &pv, overlap_orb_alpha, this->ld.phialpha);

        if (inp.deepks_out_unittest)
        {
            DeePKS_domain::check_phialpha(inp.cal_force, ucell, orb, gd, &pv, this->ld.phialpha, GlobalV::MY_RANK);
        }
    }
#endif
}

template <typename TK>
void Setup_DeePKS<TK>::before_runner(const UnitCell& ucell,    // unitcell
                                     const int nks,            // number of k points
                                     const LCAO_Orbitals& orb, // orbital info
                                     Parallel_Orbitals& pv,    // parallel orbitals
                                     const Input_para& inp)
{
#ifdef __MLALGO
    LCAO_domain::DeePKS_init(ucell, pv, nks, orb, this->ld, GlobalV::ofs_running);
    if (inp.deepks_scf)
    {
        // load the DeePKS model from deep neural network
        DeePKS_domain::load_model(inp.deepks_model, this->ld.model_deepks);
        // read pdm from file for NSCF or SCF-restart, do it only once in whole calculation
        DeePKS_domain::read_pdm((inp.init_chg == "file"),
                                inp.deepks_equiv,
                                this->ld.init_pdm,
                                ucell.nat,
                                this->ld.deepks_param,
                                *orb.Alpha,
                                this->ld.pdm);
    }
#endif
}

template <typename TK>
void Setup_DeePKS<TK>::delta_e(const UnitCell& ucell,
                               const K_Vectors& kv,
                               const LCAO_Orbitals& orb,
                               const Parallel_Orbitals& pv, // parallel orbitals
                               const Grid_Driver& gd,
                               const std::vector<std::vector<TK>>& dm_vec,
                               elecstate::fenergy& f_en,
                               const Input_para& inp)
{
#ifdef __MLALGO
    if (inp.deepks_scf)
    {
        this->ld.dpks_cal_e_delta_band(dm_vec, kv.get_nks());
        DeePKS_domain::update_dmr(kv.kvec_d, dm_vec, ucell, orb, pv, gd, this->ld.dm_r);
        f_en.edeepks_scf = this->ld.E_delta - this->ld.e_delta_band;
        f_en.edeepks_delta = this->ld.E_delta;
    }
#endif
}

template <typename TK>
void Setup_DeePKS<TK>::write_forces(const ModuleBase::matrix& fcs,
                                    const ModuleBase::matrix& fvnl_dalpha,
                                    const Input_para& inp)
{
#ifdef __MLALGO
    // DeePKS force
    if (inp.deepks_out_labels) // not parallelized yet
    {
        if (inp.deepks_out_base == "none" || (inp.deepks_out_base != "none" && this->dpks_out_type == "tot"))
        {
            const std::string file_ftot
                = PARAM.globalv.global_out_dir + (inp.deepks_out_labels == 1 ? "deepks_ftot.npy" : "deepks_force.npy");
            LCAO_deepks_io::save_matrix2npy(file_ftot, fcs, GlobalV::MY_RANK); // Hartree/Bohr, F_tot

            if (inp.deepks_out_labels == 1)
            {
                // this base only considers subtracting the deepks_scf part
                const std::string file_fbase = PARAM.globalv.global_out_dir + "deepks_fbase.npy";
                if (inp.deepks_scf)
                {
                    LCAO_deepks_io::save_matrix2npy(file_fbase,
                                                    fcs - fvnl_dalpha,
                                                    GlobalV::MY_RANK); // Hartree/Bohr, F_base
                }
                else
                {
                    LCAO_deepks_io::save_matrix2npy(file_fbase, fcs, GlobalV::MY_RANK); // no scf, F_base=F_tot
                }
            }
        }
        if (inp.deepks_out_base != "none")
        {
            // output fcs as tot or base in another dir
            // this base considers changing xc functional to base functional
            const std::string file_f
                = PARAM.globalv.global_deepks_label_elec_dir + (dpks_out_type == "tot" ? "ftot.npy" : "fbase.npy");
            LCAO_deepks_io::save_matrix2npy(file_f, fcs, GlobalV::MY_RANK);
        }
    }
#endif
}

template <typename TK>
void Setup_DeePKS<TK>::write_stress(const ModuleBase::matrix& scs,
                                    const ModuleBase::matrix& svnl_dalpha,
                                    const double& omega,
                                    const Input_para& inp)
{
#ifdef __MLALGO
    if (inp.deepks_out_labels == 1)
    {
        assert(omega > 0.0);

        if (inp.deepks_out_base == "none" || (inp.deepks_out_base != "none" && this->dpks_out_type == "tot"))
        {
            const std::string file_stot = PARAM.globalv.global_out_dir + "deepks_stot.npy";
            LCAO_deepks_io::save_matrix2npy(file_stot,
                                            scs,
                                            GlobalV::MY_RANK,
                                            omega,
                                            'U'); // change to energy unit Ry when printing, S_tot;

            // this base only considers subtracting the deepks_scf part
            const std::string file_sbase = PARAM.globalv.global_out_dir + "deepks_sbase.npy";
            if (inp.deepks_scf)
            {
                LCAO_deepks_io::save_matrix2npy(file_sbase,
                                                scs - svnl_dalpha,
                                                GlobalV::MY_RANK,
                                                omega,
                                                'U'); // change to energy unit Ry when printing, S_base;
            }
            else
            {
                LCAO_deepks_io::save_matrix2npy(file_sbase, scs, GlobalV::MY_RANK, omega,
                                                'U'); // sbase = stot
            }
        }
        if (inp.deepks_out_base != "none")
        {
            // output scs as tot or base in another dir
            // this base considers changing xc functional to base functional
            const std::string file_s = PARAM.globalv.global_deepks_label_elec_dir
                                       + (this->dpks_out_type == "tot" ? "stot.npy" : "sbase.npy");
            LCAO_deepks_io::save_matrix2npy(file_s,
                                            scs,
                                            GlobalV::MY_RANK,
                                            omega,
                                            'U'); // change to energy unit Ry when printing, S_tot;
        }
    }
    else if (inp.deepks_out_labels == 2)
    {
        const std::string file_stot = PARAM.globalv.global_out_dir + "deepks_stress.npy";
        LCAO_deepks_io::save_matrix2npy(file_stot, scs, GlobalV::MY_RANK, omega,
                                        'F'); // flat mode
    }
#endif
}

template class Setup_DeePKS<double>;
template class Setup_DeePKS<std::complex<double>>;

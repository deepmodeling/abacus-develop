#include "esolver_gets.h"

#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/read_pp_ucell.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/param_update.h"
#include "source_io/module_hs/cal_r_overlap_r.h"
#include "source_io/module_hs/write_hs_r.h"
#include "source_io/module_output/print_info.h"
#include "source_lcao/lcao_domain.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"

namespace ModuleESolver
{

ESolver_GetS::ESolver_GetS()
{
    this->classname = "ESolver_GetS";
    this->basisname = "LCAO";
}

ESolver_GetS::~ESolver_GetS()
{
}

void ESolver_GetS::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    this->inp_ = &inp;

    ModuleBase::TITLE("ESolver_GetS", "before_all_runners");
    ModuleBase::timer::start("ESolver_GetS", "before_all_runners");

    // 1.1) read pseudopotentials
    const std::string pseudo_dir = this->inp_->pseudo_dir;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool out_element_info = this->inp_->out_element_info;
    const std::string dft_functional = this->inp_->dft_functional;
    const bool lspinorb = this->inp_->lspinorb;
    const double pseudo_rcut = this->inp_->pseudo_rcut;
    const double soc_lambda = this->inp_->soc_lambda;
    const int nspin = this->inp_->nspin;
    const int npol = PARAM.globalv.npol;
    const std::string basis_type = this->inp_->basis_type;
    const std::string esolver_type = this->inp_->esolver_type;
    const std::string init_wfc = this->inp_->init_wfc;
    const int nbands = this->inp_->nbands;
    const bool two_fermi = PARAM.globalv.two_fermi;
    const double nelec_delta = this->inp_->nelec_delta;
    const std::string smearing_method = this->inp_->smearing_method;
    const std::string ks_solver = this->inp_->ks_solver;
    const int bndpar = this->inp_->bndpar;
    const double nelec = this->inp_->nelec;
    const double nupdown = this->inp_->nupdown;
    // nlocal is calculated inside read_pseudo() via CalAtomsInfo::cal_atoms_info()
    auto atoms_info = unitcell::read_pseudo(GlobalV::ofs_running,
                                            ucell,
                                            pseudo_dir,
                                            global_out_dir,
                                            out_element_info,
                                            dft_functional,
                                            lspinorb,
                                            pseudo_rcut,
                                            soc_lambda,
                                            nspin,
                                            npol,
                                            basis_type,
                                            esolver_type,
                                            init_wfc,
                                            nbands,
                                            two_fermi,
                                            nelec_delta,
                                            smearing_method,
                                            ks_solver,
                                            bndpar,
                                            nelec,
                                            nupdown);
    elecstate::ParamUpdater::update_from_atoms_info(atoms_info);

    // 1.2) symmetrize things
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        const int cal_symm_repr[2] = {this->inp_->cal_symm_repr[0], this->inp_->cal_symm_repr[1]};
        ucell.symm.analy_sys(ucell.lat,
                             ucell.st,
                             ucell.atoms,
                             GlobalV::ofs_running,
                             this->inp_->symmetry_prec,
                             inp.nspin,
                             this->inp_->calculation,
                             cal_symm_repr);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // 1.3) Setup k-points according to symmetry.
    const bool use_ibz = !inp.berry_phase && ModuleSymmetry::Symmetry::symm_flag != -1;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {this->inp_->kspacing[0], this->inp_->kspacing[1], this->inp_->kspacing[2]};
    const std::string kmesh_type = this->inp_->kmesh_type;
    const double koffset[3] = {this->inp_->koffset[0], this->inp_->koffset[1], this->inp_->koffset[2]};
    this->kv.set(ucell,
                 ucell.symm,
                 inp.kpoint_file,
                 inp.nspin,
                 ucell.G,
                 ucell.latvec,
                 GlobalV::ofs_running,
                 use_ibz,
                 global_out_dir,
                 gamma_only_local,
                 kspacing,
                 kmesh_type,
                 koffset);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    ModuleIO::print_parameters(ucell, this->kv, inp);

    // 2) init ElecState
    // autoset nbands in ElecState, it should before basis_init (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and std::complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<std::complex<double>>(&(this->chr), // use which parameter?
                                                                         &(this->kv),
                                                                         this->kv.get_nks(),
                                                                         this->pw_big);
    }

    // 3) init LCAO basis
    // reading the localized orbitals/projectors
    // construct the interpolation tables.
    LCAO_domain::init_basis_lcao(this->pv,
                                 inp.onsite_radius,
                                 inp.lcao_ecut,
                                 inp.lcao_dk,
                                 inp.lcao_dr,
                                 inp.lcao_rmax,
                                 ucell,
                                 two_center_bundle_,
                                 orb_);

    ModuleBase::timer::end("ESolver_GetS", "before_all_runners");
}

void ESolver_GetS::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_GetS", "runner");
    ModuleBase::timer::start("ESolver_GetS", "runner");

    // (1) Find adjacent atoms for each atom.
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                            this->inp_->out_level,
                                            orb_.get_rcutmax_Phi(),
                                            ucell.infoNL->get_rcutmax_Beta(),
                                            PARAM.globalv.gamma_only_local);

    Grid_Driver gd;

    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         gd,
                         ucell,
                         search_radius,
                         this->inp_->test_atom_input);

    Record_adj RA;
    RA.for_2d(ucell, gd, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    if (this->p_hamilt == nullptr)
    {
        if (this->inp_->nspin == 4)
        {
            this->p_hamilt
                = new hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>(ucell,
                                                                                     gd,
                                                                                     &this->pv,
                                                                                     this->kv,
                                                                                     *(two_center_bundle_.overlap_orb),
                                                                                     orb_.cutoffs());
            auto* hamilt_ptr = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
            auto* ops_ptr
                = dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>*>(hamilt_ptr->ops);
            ops_ptr->contributeHR();
        }
        else
        {
            this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>, double>(ucell,
                                                                                  gd,
                                                                                  &this->pv,
                                                                                  this->kv,
                                                                                  *(two_center_bundle_.overlap_orb),
                                                                                  orb_.cutoffs());
            auto* hamilt_ptr = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
            auto* ops_ptr = dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, double>*>(hamilt_ptr->ops);
            ops_ptr->contributeHR();
        }
    }

    const std::string fn = PARAM.globalv.global_out_dir + "sr_nao.csr";

    auto* hamilt_ptr = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
    ModuleIO::output_SR(pv, gd, hamilt_ptr, fn);

    if (this->inp_->out_mat_r[0])
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(ucell, pv, orb_);
        r_matrix.out_rR(ucell, gd, istep, this->inp_->out_mat_r[1]);
    }

    if (this->inp_->out_mat_ds[0])
    {
        LCAO_HS_Arrays HS_Arrays; // store sparse arrays
        //! Print out sparse matrix
        ModuleIO::output_dSR(istep,
                             ucell,
                             this->pv,
                             HS_Arrays,
                             gd, // mohan add 2024-04-06
                             two_center_bundle_,
                             orb_,
                             kv,
                             false,
                             1e-10,
                             this->inp_->out_mat_ds[1]);
    }

    ModuleBase::timer::end("ESolver_GetS", "runner");
}

void ESolver_GetS::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
};
double ESolver_GetS::cal_energy()
{
    return 0.0;
};
void ESolver_GetS::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
};
void ESolver_GetS::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
};

} // namespace ModuleESolver

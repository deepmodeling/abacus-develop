#include "esolver_gets.h"

#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/read_pseudo.h"
#include "source_lcao/LCAO_domain.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_io/module_hs/cal_r_overlap_R.h"
#include "source_io/module_output/print_info.h"
#include "source_io/module_hs/write_HS_R.h"

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

void ESolver_GetS::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_GetS", "before_all_runners");
    ModuleBase::timer::start("ESolver_GetS", "before_all_runners");

    // 1.1) read pseudopotentials
    const std::string pseudo_dir = PARAM.inp.pseudo_dir;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool out_element_info = PARAM.inp.out_element_info;
    const std::string dft_functional = PARAM.inp.dft_functional;
    const bool lspinorb = PARAM.inp.lspinorb;
    const double pseudo_rcut = PARAM.inp.pseudo_rcut;
    const double soc_lambda = PARAM.inp.soc_lambda;
    const int nspin = PARAM.inp.nspin;
    const int nlocal = PARAM.globalv.nlocal;
    const int npol = PARAM.globalv.npol;
    const std::string basis_type = PARAM.inp.basis_type;
    const std::string esolver_type = PARAM.inp.esolver_type;
    const std::string init_wfc = PARAM.inp.init_wfc;
    const int nbands = PARAM.inp.nbands;
    elecstate::read_pseudo(GlobalV::ofs_running, ucell, pseudo_dir, global_out_dir, out_element_info, dft_functional, lspinorb, pseudo_rcut, soc_lambda, nspin, nlocal, npol, basis_type, esolver_type, init_wfc, nbands);

    // 1.2) symmetrize things
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // 1.3) Setup k-points according to symmetry.
    const bool use_ibz = !inp.berry_phase && ModuleSymmetry::Symmetry::symm_flag != -1;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {PARAM.inp.kspacing[0], PARAM.inp.kspacing[1], PARAM.inp.kspacing[2]};
    const std::string kmesh_type = PARAM.inp.kmesh_type;
    const double koffset[3] = {PARAM.inp.koffset[0], PARAM.inp.koffset[1], PARAM.inp.koffset[2]};
    this->kv.set(ucell, ucell.symm, inp.kpoint_file, inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running, use_ibz, global_out_dir, gamma_only_local, kspacing, kmesh_type, koffset);
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

void ESolver_GetS::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_GetS", "runner");
    ModuleBase::timer::start("ESolver_GetS", "runner");

    // (1) Find adjacent atoms for each atom.
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                            PARAM.inp.out_level,
                                            orb_.get_rcutmax_Phi(),
                                            ucell.infoNL->get_rcutmax_Beta(),
                                            PARAM.globalv.gamma_only_local);

    Grid_Driver gd;

    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         gd,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    Record_adj RA;
    RA.for_2d(ucell, gd, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    if (this->p_hamilt == nullptr)
    {
        if (PARAM.inp.nspin == 4)
        {
            this->p_hamilt
                = new hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>(ucell,
                                                                                     gd,
                                                                                     &this->pv,
                                                                                     this->kv,
                                                                                     *(two_center_bundle_.overlap_orb),
                                                                                     orb_.cutoffs());
            auto* hamilt_ptr = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
            auto* ops_ptr = dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>*>(hamilt_ptr->ops);
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

    if (PARAM.inp.out_mat_r[0])
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(ucell, pv, orb_);
        r_matrix.out_rR(ucell, gd, istep, PARAM.inp.out_mat_r[1]);
    }

    if (PARAM.inp.out_mat_ds[0])
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
                             PARAM.inp.out_mat_ds[1]);
    }

    ModuleBase::timer::end("ESolver_GetS", "runner");
}

void ESolver_GetS::after_all_runners(UnitCell& ucell) {};
double ESolver_GetS::cal_energy()
{
    return 0.0;
};
void ESolver_GetS::cal_force(UnitCell& ucell, ModuleBase::matrix& force) {};
void ESolver_GetS::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) {};

} // namespace ModuleESolver

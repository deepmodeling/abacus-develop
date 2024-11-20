#include "esolver_gets.h"

#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_io/print_info.h"
#include "module_io/write_HS_R.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_GetS<TK, TR>::ESolver_GetS()
{
    this->classname = "ESolver_GetS";
    this->basisname = "LCAO";
}

template <typename TK, typename TR>
ESolver_GetS<TK, TR>::~ESolver_GetS()
{
}

template <typename TK, typename TR>
void ESolver_GetS<TK, TR>::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_GetS", "before_all_runners");
    ModuleBase::timer::tick("ESolver_GetS", "before_all_runners");

    // 1.1) read pseudopotentials
    ucell.read_pseudo(GlobalV::ofs_running);

    // 1.2) symmetrize things
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // 1.3) Setup k-points according to symmetry.
    this->kv.set(ucell.symm, inp.kpoint_file, inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    ModuleIO::setup_parameters(ucell, this->kv);

    // 2) init ElecState
    // autoset nbands in ElecState, it should before basis_init (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), // use which parameter?
                                                       &(this->kv),
                                                       this->kv.get_nks(),
                                                       &(this->GG), // mohan add 2024-04-01
                                                       &(this->GK), // mohan add 2024-04-01
                                                       this->pw_rho,
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

    // 4) initialize the density matrix
    // DensityMatrix is allocated here, DMK is also initialized here
    // DMR is not initialized here, it will be constructed in each before_scf
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->init_DM(&this->kv, &(this->pv), inp.nspin);

    ModuleBase::timer::tick("ESolver_GetS", "before_all_runners");
}

template <>
void ESolver_GetS<double, double>::runner(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_GetS", "runner");
    ModuleBase::WARNING_QUIT("ESolver_GetS<double, double>::runner", "not implemented");
}

template <>
void ESolver_GetS<std::complex<double>, std::complex<double>>::runner(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_GetS", "runner");
    ModuleBase::timer::tick("ESolver_GetS", "runner");

    // (1) Find adjacent atoms for each atom.
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                            PARAM.inp.out_level,
                                            orb_.get_rcutmax_Phi(),
                                            ucell.infoNL.get_rcutmax_Beta(),
                                            PARAM.globalv.gamma_only_local);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    this->RA.for_2d(this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    if (this->p_hamilt == nullptr)
    {
        this->p_hamilt
            = new hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>(&this->pv,
                                                                                 this->kv,
                                                                                 *(two_center_bundle_.overlap_orb),
                                                                                 orb_.cutoffs());
        dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt->ops)
            ->contributeHR();
    }

    const std::string fn = PARAM.globalv.global_out_dir + "SR.csr";
    std::cout << " The file is saved in " << fn << std::endl;
    ModuleIO::output_SR(pv, GlobalC::GridD, this->p_hamilt, fn);

    ModuleBase::timer::tick("ESolver_GetS", "runner");
}

template <>
void ESolver_GetS<std::complex<double>, double>::runner(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_GetS", "runner");
    ModuleBase::timer::tick("ESolver_GetS", "runner");

    // (1) Find adjacent atoms for each atom.
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                            PARAM.inp.out_level,
                                            orb_.get_rcutmax_Phi(),
                                            ucell.infoNL.get_rcutmax_Beta(),
                                            PARAM.globalv.gamma_only_local);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    this->RA.for_2d(this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    if (this->p_hamilt == nullptr)
    {
        this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>, double>(&this->pv,
                                                                              this->kv,
                                                                              *(two_center_bundle_.overlap_orb),
                                                                              orb_.cutoffs());
        dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, double>*>(this->p_hamilt->ops)->contributeHR();
    }

    const std::string fn = PARAM.globalv.global_out_dir + "SR.csr";
    std::cout << " The file is saved in " << fn << std::endl;
    ModuleIO::output_SR(pv, GlobalC::GridD, this->p_hamilt, fn);

    ModuleBase::timer::tick("ESolver_GetS", "runner");
}

template class ESolver_GetS<double, double>;
template class ESolver_GetS<std::complex<double>, double>;
template class ESolver_GetS<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver

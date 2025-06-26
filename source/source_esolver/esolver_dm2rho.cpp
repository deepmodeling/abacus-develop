#include "esolver_dm2rho.h"

#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/read_pseudo.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_io/cube_io.h"
#include "module_io/io_npz.h"
#include "module_io/print_info.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_DM2rho<TK, TR>::ESolver_DM2rho()
{
    this->classname = "ESolver_DM2rho";
    this->basisname = "LCAO";
}

template <typename TK, typename TR>
ESolver_DM2rho<TK, TR>::~ESolver_DM2rho()
{
}

template <typename TK, typename TR>
void ESolver_DM2rho<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_DM2rho", "before_all_runners");
    ModuleBase::timer::tick("ESolver_DM2rho", "before_all_runners");

    ESolver_KS_LCAO<TK, TR>::before_all_runners(ucell, inp);

    ModuleBase::timer::tick("ESolver_DM2rho", "before_all_runners");
}

template <typename TK, typename TR>
void ESolver_DM2rho<TK, TR>::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_DM2rho", "runner");
    ModuleBase::timer::tick("ESolver_DM2rho", "runner");

    ESolver_KS_LCAO<TK, TR>::before_scf(ucell, istep);

    // file name of DM
    std::string zipname = "output_DM0.npz";
    elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();

    // read DM from file
    ModuleIO::read_mat_npz(&(this->pv), ucell, zipname, *(dm->get_DMR_pointer(1)));

    // if nspin=2, need extra reading
    if (PARAM.inp.nspin == 2)
    {
        zipname = "output_DM1.npz";
        ModuleIO::read_mat_npz(&(this->pv), ucell, zipname, *(dm->get_DMR_pointer(2)));
    }

    this->pelec->psiToRho(*this->psi);

    int nspin0 = PARAM.inp.nspin == 2 ? 2 : 1;

    for (int is = 0; is < nspin0; is++)
    {
        std::string fn = PARAM.globalv.global_out_dir + "/SPIN" + std::to_string(is + 1) + "_CHG.cube";

        // write electron density
        ModuleIO::write_vdata_palgrid(this->Pgrid,
                                      this->chr.rho[is],
                                      is,
                                      PARAM.inp.nspin,
                                      istep,
                                      fn,
                                      this->pelec->eferm.get_efval(is),
                                      &(ucell),
                                      3,
                                      1);
    }

    ModuleBase::timer::tick("ESolver_DM2rho", "runner");
}

template <typename TK, typename TR>
void ESolver_DM2rho<TK, TR>::after_all_runners(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DM2rho", "after_all_runners");
    ModuleBase::timer::tick("ESolver_DM2rho", "after_all_runners");

    ESolver_KS_LCAO<TK, TR>::after_all_runners(ucell);

    ModuleBase::timer::tick("ESolver_DM2rho", "after_all_runners");
};

template class ESolver_DM2rho<std::complex<double>, double>;
template class ESolver_DM2rho<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver

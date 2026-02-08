#include "ctrl_output_td.h"

#include "source_base/parallel_global.h"
#include "source_io/dipole_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_current/td_current_io.h"

namespace ModuleIO
{

template <typename TR>
void ctrl_output_td(const UnitCell& ucell,
                    double** rho_save,
                    const ModulePW::PW_Basis* rhopw,
                    const int istep,
                    const psi::Psi<std::complex<double>>* psi,
                    const elecstate::ElecState* pelec,
                    const K_Vectors& kv,
                    const TwoCenterIntegrator* intor,
                    const Parallel_Orbitals* pv,
                    const LCAO_Orbitals& orb,
                    const Velocity_op<TR>* velocity_mat,
                    const Grid_Driver& grid,
                    hamilt::HamiltLCAO<std::complex<double>, TR>* p_hamilt,
                    Record_adj& RA,
                    TD_info* td_p,
                    const Exx_NAO<std::complex<double>>& exx_nao
                    )
{
    ModuleBase::TITLE("ModuleIO", "ctrl_output_td");

#ifdef __LCAO
    // (1) Write dipole information
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        if (PARAM.inp.out_dipole == 1)
        {
            std::stringstream ss_dipole;
            ss_dipole << PARAM.globalv.global_out_dir << "dipole_s" << is + 1 << ".txt";
            ModuleIO::write_dipole(ucell, rho_save[is], rhopw, is, istep, ss_dipole.str());
        }
    }

    // (2) Write current information
    const elecstate::ElecStateLCAO<std::complex<double>>* pelec_lcao
        = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec);

    if (!pelec_lcao)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::ctrl_output_td", "Failed to cast ElecState to ElecStateLCAO");
    }

    if (TD_info::out_current == 1)
    {
        if (TD_info::out_current_k)
        {
            ModuleIO::write_current_eachk<TR>(ucell, istep, psi, pelec, kv, intor, pv, orb, velocity_mat, RA);
        }
        else
        {
            ModuleIO::write_current<TR>(ucell, istep, psi, pelec, kv, intor, pv, orb, velocity_mat, RA);
        }
    }
    else if(TD_info::out_current==2)
    {
        ModuleIO::write_current(ucell, grid, istep, psi, pelec, kv, pv, orb, td_p->r_calculator, p_hamilt->getSR(), p_hamilt->getHR(), exx_nao);
    }

    // (3) Output file for restart
    if (PARAM.inp.out_freq_td > 0) // default value of out_freq_td is 0
    {
        if (istep % PARAM.inp.out_freq_td == 0)
        {
            if (td_p != nullptr)
            {
                td_p->out_restart_info(istep, elecstate::H_TDDFT_pw::At, elecstate::H_TDDFT_pw::At_laststep);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ModuleIO::ctrl_output_td",
                                         "TD_info pointer is null, cannot output restart info.");
            }
        }
    }
#endif // __LCAO
}

template void ctrl_output_td<double>(const UnitCell&,
                                     double**,
                                     const ModulePW::PW_Basis*,
                                     const int,
                                     const psi::Psi<std::complex<double>>*,
                                     const elecstate::ElecState*,
                                     const K_Vectors&,
                                     const TwoCenterIntegrator*,
                                     const Parallel_Orbitals*,
                                     const LCAO_Orbitals&,
                                     const Velocity_op<double>*,
                                     const Grid_Driver&,
                                     hamilt::HamiltLCAO<std::complex<double>, double>*,
                                     Record_adj&,
                                     TD_info*,
                                     const Exx_NAO<std::complex<double>>&
                                     );

template void ctrl_output_td<std::complex<double>>(const UnitCell&,
                                                   double**,
                                                   const ModulePW::PW_Basis*,
                                                   const int,
                                                   const psi::Psi<std::complex<double>>*,
                                                   const elecstate::ElecState*,
                                                   const K_Vectors&,
                                                   const TwoCenterIntegrator*,
                                                   const Parallel_Orbitals*,
                                                   const LCAO_Orbitals&,
                                                   const Velocity_op<std::complex<double>>*,
                                                   const Grid_Driver&,
                                                   hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*,
                                                   Record_adj&,
                                                   TD_info*,
                                                   const Exx_NAO<std::complex<double>>&
                                                   );

} // namespace ModuleIO
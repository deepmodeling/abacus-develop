#include "evolve_elec.h"

#include "evolve_psi.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace module_tddft
{
Evolve_elec::Evolve_elec() {};
Evolve_elec::~Evolve_elec() {};

double Evolve_elec::td_force_dt;
bool Evolve_elec::td_vext;
std::vector<int> Evolve_elec::td_vext_dire_case;
bool Evolve_elec::out_dipole;
bool Evolve_elec::out_efield;
double Evolve_elec::td_print_eij; // the threshold to output Eij elements
int Evolve_elec::td_edm;          // 0: new edm method   1: old edm method

// this routine only serves for TDDFT using LCAO basis set
void Evolve_elec::solve_psi(const int& istep,
                            const int nband,
                            const int nlocal,
                            hamilt::Hamilt<std::complex<double>>* phm,
                            Parallel_Orbitals& para_orb,
                            psi::Psi<std::complex<double>>* psi,
                            psi::Psi<std::complex<double>>* psi_laststep,
                            std::complex<double>** Hk_laststep,
                            std::complex<double>** Sk_laststep,
                            ModuleBase::matrix& ekb,
                            int htype,
                            int propagator,
                            const int& nks)
{
    ModuleBase::TITLE("Evolve_elec", "solve_psi");
    ModuleBase::timer::tick("Evolve_elec", "solve_psi");

    for (int ik = 0; ik < nks; ik++)
    {
        phm->updateHk(ik);

        ModuleBase::timer::tick("Efficience", "evolve_k");
        psi->fix_k(ik);
        psi_laststep->fix_k(ik);
        if (htype == 0)
        {
            evolve_psi(nband,
                       nlocal,
                       &(para_orb),
                       phm,
                       psi[0].get_pointer(),
                       psi_laststep[0].get_pointer(),
                       nullptr,
                       nullptr,
                       &(ekb(ik, 0)),
                       htype,
                       propagator);
        }
        else if (htype == 1)
        {
            evolve_psi(nband,
                       nlocal,
                       &(para_orb),
                       phm,
                       psi[0].get_pointer(),
                       psi_laststep[0].get_pointer(),
                       Hk_laststep[ik],
                       Sk_laststep[ik],
                       &(ekb(ik, 0)),
                       htype,
                       propagator);

            const bool use_tensor = false;
            if (use_tensor)
            {
                std::cout << "Print ekb: " << std::endl;
                ekb.print(std::cout);
                std::cout << "nband = " << nband << std::endl;
                std::cout << "psi->get_nbands() = " << psi->get_nbands() << std::endl;
                std::cout << "nlocal = " << nlocal << std::endl;
                std::cout << "psi->get_nbasis() = " << psi->get_nbasis() << std::endl;
                std::cout << "ekb.nr = " << ekb.nr << std::endl;
                std::cout << "ekb.nc = " << ekb.nc << std::endl;

                // Create TensorMap for psi_k, psi_k_laststep, H_laststep, S_laststep, ekb
                container::TensorMap psi_k_tensor(psi[0].get_pointer(),
                                                  container::DataType::DT_COMPLEX_DOUBLE,
                                                  container::DeviceType::CpuDevice,
                                                  container::TensorShape({psi->get_nbands(), psi->get_nbasis()}));
                container::TensorMap psi_k_laststep_tensor(
                    psi_laststep[0].get_pointer(),
                    container::DataType::DT_COMPLEX_DOUBLE,
                    container::DeviceType::CpuDevice,
                    container::TensorShape({psi->get_nbands(), psi->get_nbasis()}));
                container::TensorMap H_laststep_tensor(Hk_laststep[ik],
                                                       container::DataType::DT_COMPLEX_DOUBLE,
                                                       container::DeviceType::CpuDevice,
                                                       container::TensorShape({para_orb.nloc}));
                container::TensorMap S_laststep_tensor(Sk_laststep[ik],
                                                       container::DataType::DT_COMPLEX_DOUBLE,
                                                       container::DeviceType::CpuDevice,
                                                       container::TensorShape({para_orb.nloc}));
                container::TensorMap ekb_tensor(&(ekb(ik, 0)),
                                                container::DataType::DT_DOUBLE,
                                                container::DeviceType::CpuDevice,
                                                container::TensorShape({nband}));

                evolve_psi_tensor(nband,
                                  nlocal,
                                  &(para_orb),
                                  phm,
                                  psi_k_tensor,
                                  psi_k_laststep_tensor,
                                  H_laststep_tensor,
                                  S_laststep_tensor,
                                  ekb_tensor,
                                  htype,
                                  propagator);
                // evolve_psi_tensor(nband,
                //                   nlocal,
                //                   &(para_orb),
                //                   phm,
                //                   psi[0].get_pointer(),
                //                   psi_laststep[0].get_pointer(),
                //                   Hk_laststep[ik],
                //                   Sk_laststep[ik],
                //                   &(ekb(ik, 0)),
                //                   htype,
                //                   propagator);
                std::cout << "Print ekb tensor: " << std::endl;
                ekb.print(std::cout);
            }
        }
        else
        {
            std::cout << "method of htype is wrong" << std::endl;
        }

        ModuleBase::timer::tick("Efficience", "evolve_k");
    } // end k

    ModuleBase::timer::tick("Evolve_elec", "solve_psi");
    return;
}
} // namespace module_tddft
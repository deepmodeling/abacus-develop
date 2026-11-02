#include "source_pw/module_rt/td_endpoint_pw.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_hamilt/hamilt_hs_adapter.h"

namespace pw
{
template <typename T, typename Device>
void orth_init_pw(hsolver::TDOrbOps<T, Device>& orbitals,
                  hsolver::HSOperator<T, Device>& op,
                  psi::Psi<T, Device>& current,
                  elecstate::ElecStatePW<T, Device>& estate,
                  Charge& charge,
                  UnitCell& cell,
                  const ModulePW::PW_Basis& rhod,
                  int nspin,
                  int nbands,
                  int iter,
                  bool do_orth)
{
    ModuleBase::timer::start("pw", "orth_init_pw");
    // Padding remains zero even when physical orbital corrections are disabled.
    orbitals.zero_padding(&current);
    if (do_orth)
    {
        int changed = orbitals.orth_initial(&current, iter);
        // Occupations and density construction contain collectives across pools.
        Parallel_Reduce::reduce_all(changed);
        if (changed > 0)
        {
            orbitals.cal_band_energy(op, current, &estate.ekb);
            elecstate::calculate_weights(estate.ekb,
                                         estate.wg,
                                         estate.klist,
                                         estate.eferm,
                                         estate.f_en,
                                         estate.nelec_spin,
                                         nbands,
                                         estate.skip_weights);
            elecstate::calEBand(estate.ekb, estate.wg, estate.f_en);

            estate.psiToRho(current);
            module_charge::symmetrize_rho(nspin, charge, &rhod, cell.symm);
        }
    }
    ModuleBase::timer::end("pw", "orth_init_pw");
}

template <typename T, typename Device>
void rebuild_td_end(psi::Psi<T, Device>& current,
                    elecstate::ElecStatePW<T, Device>& estate,
                    Charge& charge,
                    UnitCell& cell,
                    hamilt::HamiltPW<T, Device>& hamilt,
                    hsolver::TDOrbOps<T, Device>& orbitals,
                    const ModulePW::PW_Basis_K& basis,
                    const ModulePW::PW_Basis& rhod,
                    const ModuleBase::Vector3<double>& vecpot,
                    int nspin,
                    int imp_sol,
                    bool sc_mag,
                    int dft_u,
                    const std::string& assume_iso)
{
    ModuleBase::timer::start("pw", "rebuild_td_end");
    // Discard a mixed trial density while keeping the final orbitals and occupations.
    estate.psiToRho(current);
    module_charge::symmetrize_rho(nspin, charge, &rhod, cell.symm);
    charge.save_rho_before_sum_band();
    estate.pot->update_from_charge(&charge, &cell);

    using Real = typename GetTypeReal<T>::type;
    const Real* veff = estate.pot->template get_veff_smooth_data<Real>();
    hamilt.bind_td_state(veff, vecpot);
    hamilt::HamiltHSOperator<T, Device> op(&hamilt, &basis);
    orbitals.cal_band_energy(op, current, &estate.ekb);
    elecstate::calEBand(estate.ekb, estate.wg, estate.f_en);
    estate.f_en.deband = estate.cal_delta_eband(cell);
    estate.f_en.descf = 0.0;

    // Do not reuse a force correction from an earlier self-consistent endpoint.
    estate.vnew_exist = false;
    estate.cal_energies(2, imp_sol, sc_mag, dft_u, assume_iso);
    ModuleBase::timer::end("pw", "rebuild_td_end");
}

template void orth_init_pw(hsolver::TDOrbOps<std::complex<float>, base_device::DEVICE_CPU>&,
                           hsolver::HSOperator<std::complex<float>, base_device::DEVICE_CPU>&,
                           psi::Psi<std::complex<float>, base_device::DEVICE_CPU>&,
                           elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>&,
                           Charge&,
                           UnitCell&,
                           const ModulePW::PW_Basis&,
                           int,
                           int,
                           int,
                           bool);
template void rebuild_td_end(psi::Psi<std::complex<float>, base_device::DEVICE_CPU>&,
                             elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>&,
                             Charge&,
                             UnitCell&,
                             hamilt::HamiltPW<std::complex<float>, base_device::DEVICE_CPU>&,
                             hsolver::TDOrbOps<std::complex<float>, base_device::DEVICE_CPU>&,
                             const ModulePW::PW_Basis_K&,
                             const ModulePW::PW_Basis&,
                             const ModuleBase::Vector3<double>&,
                             int,
                             int,
                             bool,
                             int,
                             const std::string&);
template void orth_init_pw(hsolver::TDOrbOps<std::complex<double>, base_device::DEVICE_CPU>&,
                           hsolver::HSOperator<std::complex<double>, base_device::DEVICE_CPU>&,
                           psi::Psi<std::complex<double>, base_device::DEVICE_CPU>&,
                           elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>&,
                           Charge&,
                           UnitCell&,
                           const ModulePW::PW_Basis&,
                           int,
                           int,
                           int,
                           bool);
template void rebuild_td_end(psi::Psi<std::complex<double>, base_device::DEVICE_CPU>&,
                             elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>&,
                             Charge&,
                             UnitCell&,
                             hamilt::HamiltPW<std::complex<double>, base_device::DEVICE_CPU>&,
                             hsolver::TDOrbOps<std::complex<double>, base_device::DEVICE_CPU>&,
                             const ModulePW::PW_Basis_K&,
                             const ModulePW::PW_Basis&,
                             const ModuleBase::Vector3<double>&,
                             int,
                             int,
                             bool,
                             int,
                             const std::string&);
#if defined(__CUDA) || defined(__ROCM)
template void orth_init_pw(hsolver::TDOrbOps<std::complex<float>, base_device::DEVICE_GPU>&,
                           hsolver::HSOperator<std::complex<float>, base_device::DEVICE_GPU>&,
                           psi::Psi<std::complex<float>, base_device::DEVICE_GPU>&,
                           elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>&,
                           Charge&,
                           UnitCell&,
                           const ModulePW::PW_Basis&,
                           int,
                           int,
                           int,
                           bool);
template void rebuild_td_end(psi::Psi<std::complex<float>, base_device::DEVICE_GPU>&,
                             elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>&,
                             Charge&,
                             UnitCell&,
                             hamilt::HamiltPW<std::complex<float>, base_device::DEVICE_GPU>&,
                             hsolver::TDOrbOps<std::complex<float>, base_device::DEVICE_GPU>&,
                             const ModulePW::PW_Basis_K&,
                             const ModulePW::PW_Basis&,
                             const ModuleBase::Vector3<double>&,
                             int,
                             int,
                             bool,
                             int,
                             const std::string&);
template void orth_init_pw(hsolver::TDOrbOps<std::complex<double>, base_device::DEVICE_GPU>&,
                           hsolver::HSOperator<std::complex<double>, base_device::DEVICE_GPU>&,
                           psi::Psi<std::complex<double>, base_device::DEVICE_GPU>&,
                           elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>&,
                           Charge&,
                           UnitCell&,
                           const ModulePW::PW_Basis&,
                           int,
                           int,
                           int,
                           bool);
template void rebuild_td_end(psi::Psi<std::complex<double>, base_device::DEVICE_GPU>&,
                             elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>&,
                             Charge&,
                             UnitCell&,
                             hamilt::HamiltPW<std::complex<double>, base_device::DEVICE_GPU>&,
                             hsolver::TDOrbOps<std::complex<double>, base_device::DEVICE_GPU>&,
                             const ModulePW::PW_Basis_K&,
                             const ModulePW::PW_Basis&,
                             const ModuleBase::Vector3<double>&,
                             int,
                             int,
                             bool,
                             int,
                             const std::string&);
#endif
} // namespace pw

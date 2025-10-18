#include "source_io/ctrl_runner_lcao.h" // use ctrl_runner_lcao() 

#include "source_estate/elecstate_lcao.h" // use elecstate::ElecState
#include "source_lcao/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>

#include "source_io/write_proj_band_lcao.h" // projcted band structure
#include "source_io/cal_ldos.h" // cal LDOS
#include "source_io/write_eband_terms.hpp"
#include "source_io/write_vxc.hpp"
#include "source_io/write_vxc_r.hpp"

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_runner_lcao(UnitCell& ucell,      // unitcell
        const Input_para &inp,              // input
		K_Vectors &kv,                      // k-point
		elecstate::ElecStateLCAO<TK>* pelec,// electronic info
		Parallel_Orbitals &pv,              // orbital info
        Parallel_Grid &pgrid,               // grid info
		Grid_Driver &gd,                    // search for adjacent atoms
		psi::Psi<TK>* psi,                  // wave function
        Charge &chr,                  // charge density
		hamilt::HamiltLCAO<TK, TR>* p_hamilt, // hamiltonian
		TwoCenterBundle &two_center_bundle,   // use two-center integration
        Gint_Gamma &gg,                     // gint for Gamma-only
		Gint_k &gk,                         // gint for multi k-points
		LCAO_Orbitals &orb,                 // LCAO orbitals
		ModulePW::PW_Basis* pw_rho,   // charge density
		ModulePW::PW_Basis* pw_rhod,  // dense charge density 
		Structure_Factor &sf,         // structure factor
        ModuleBase::matrix &vloc,     // local pseudopotential 
		Exx_NAO<TK> &exx_nao,
        surchem &solvent)             // solvent model
{
    ModuleBase::TITLE("ModuleIO", "ctrl_runner_lcao");
    ModuleBase::timer::tick("ModuleIO", "ctrl_runner_lcao");

    // 1) write projected band structure
    if (inp.out_proj_band)
    {
        ModuleIO::write_proj_band_lcao(psi, pv, pelec, kv, ucell, p_hamilt);
    }

	// 2) out ldos
	if (inp.out_ldos[0])
    {
        ModuleIO::Cal_ldos<TK>::cal_ldos_lcao(pelec, psi[0], pgrid, ucell);
    }

    // 3) print out exchange-correlation potential
    if (inp.out_mat_xc)
    {
        ModuleIO::write_Vxc<TK, TR>(inp.nspin,
                                    PARAM.globalv.nlocal,
                                    GlobalV::DRANK,
                                    &pv,
                                    *psi,
                                    ucell,
                                    sf,
                                    solvent,
                                    *pw_rho,
                                    *pw_rhod,
                                    vloc,
                                    chr,
                                    gg,
                                    gk,
                                    kv,
                                    orb.cutoffs(),
                                    pelec->wg,
                                    gd
#ifdef __EXX
                                    ,
                                    exx_nao.exd ? &exx_nao.exd->get_Hexxs() : nullptr,
                                    exx_nao.exc ? &exx_nao.exc->get_Hexxs() : nullptr
#endif
        );
    }

    if (inp.out_mat_xc2)
    {
        ModuleIO::write_Vxc_R<TK, TR>(inp.nspin,
                                      &pv,
                                      ucell,
                                      sf,
                                      solvent,
                                      *pw_rho,
                                      *pw_rhod,
                                      vloc,
                                      chr,
                                      gg,
                                      gk,
                                      kv,
                                      orb.cutoffs(),
                                      gd
#ifdef __EXX
                                      ,
                                      exx_nao.exd ? &exx_nao.exd->get_Hexxs() : nullptr,
                                      exx_nao.exc ? &exx_nao.exc->get_Hexxs() : nullptr
#endif
        );
    }


    // write eband terms
    if (inp.out_eband_terms)
    {
        ModuleIO::write_eband_terms<TK, TR>(inp.nspin,
                                            PARAM.globalv.nlocal,
                                            GlobalV::DRANK,
                                            &pv,
                                            *psi,
                                            ucell,
                                            sf,
                                            solvent,
                                            *pw_rho,
                                            *pw_rhod,
                                            vloc,
                                            chr,
                                            gg,
                                            gk,
                                            kv,
                                            pelec->wg,
                                            gd,
                                            orb.cutoffs(),
                                            two_center_bundle
#ifdef __EXX
                                            ,
                                            exx_nao.exd ? &exx_nao.exd->get_Hexxs() : nullptr,
                                            exx_nao.exc ? &exx_nao.exc->get_Hexxs() : nullptr
#endif
       );
    }

}




// TK: double  TR: double 
template void ModuleIO::ctrl_runner_lcao<double, double>(UnitCell& ucell,      // unitcell
        const Input_para &inp,              // input
		K_Vectors &kv,                      // k-point
		elecstate::ElecStateLCAO<double>* pelec,// electronic info
		Parallel_Orbitals &pv,              // orbital info
        Parallel_Grid &pgrid,               // grid info
		Grid_Driver &gd,                    // search for adjacent atoms
		psi::Psi<double>* psi,                  // wave function
        Charge &chr,                  // charge density
		hamilt::HamiltLCAO<double, double>* p_hamilt, // hamiltonian
		TwoCenterBundle &two_center_bundle,   // use two-center integration
        Gint_Gamma &gg,                     // gint for Gamma-only
		Gint_k &gk,                         // gint for multi k-points
		LCAO_Orbitals &orb,                 // LCAO orbitals
		ModulePW::PW_Basis* pw_rho,   // charge density
		ModulePW::PW_Basis* pw_rhod,  // dense charge density 
		Structure_Factor &sf,         // structure factor
        ModuleBase::matrix &vloc,     // local pseudopotential 
        Exx_NAO<double> &exx_nao,
        surchem &solvent);             // solvent model

// TK: complex<double>  TR: double 
template void ctrl_runner_lcao<std::complex<double>, double>(UnitCell& ucell,      // unitcell
        const Input_para &inp,              // input
		K_Vectors &kv,                      // k-point
		elecstate::ElecStateLCAO<std::complex<double>>* pelec,// electronic info
		Parallel_Orbitals &pv,              // orbital info
        Parallel_Grid &pgrid,               // grid info
		Grid_Driver &gd,                    // search for adjacent atoms
		psi::Psi<std::complex<double>>* psi,                  // wave function
        Charge &chr,                  // charge density
		hamilt::HamiltLCAO<std::complex<double>, double>* p_hamilt, // hamiltonian
		TwoCenterBundle &two_center_bundle,   // use two-center integration
        Gint_Gamma &gg,                     // gint for Gamma-only
		Gint_k &gk,                         // gint for multi k-points
		LCAO_Orbitals &orb,                 // LCAO orbitals
		ModulePW::PW_Basis* pw_rho,   // charge density
		ModulePW::PW_Basis* pw_rhod,  // dense charge density 
		Structure_Factor &sf,         // structure factor
        ModuleBase::matrix &vloc,     // local pseudopotential 
        Exx_NAO<std::complex<double>> &exx_nao,
        surchem &solvent);             // solvent model

// TK: complex<double>  TR: complex<double>
template void ctrl_runner_lcao<std::complex<double>, std::complex<double>>(UnitCell& ucell,      // unitcell
        const Input_para &inp,              // input
		K_Vectors &kv,                      // k-point
		elecstate::ElecStateLCAO<std::complex<double>>* pelec,// electronic info
		Parallel_Orbitals &pv,              // orbital info
        Parallel_Grid &pgrid,               // grid info
		Grid_Driver &gd,                    // search for adjacent atoms
		psi::Psi<std::complex<double>>* psi,                  // wave function
        Charge &chr,                  // charge density
		hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_hamilt, // hamiltonian
		TwoCenterBundle &two_center_bundle,   // use two-center integration
        Gint_Gamma &gg,                     // gint for Gamma-only
		Gint_k &gk,                         // gint for multi k-points
		LCAO_Orbitals &orb,                 // LCAO orbitals
		ModulePW::PW_Basis* pw_rho,   // charge density
		ModulePW::PW_Basis* pw_rhod,  // dense charge density 
		Structure_Factor &sf,         // structure factor
        ModuleBase::matrix &vloc,     // local pseudopotential 
        Exx_NAO<std::complex<double>> &exx_nao,
        surchem &solvent);             // solvent model

} // end namespace

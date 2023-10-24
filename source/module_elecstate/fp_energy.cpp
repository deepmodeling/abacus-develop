#include "fp_energy.h"
#include "module_base/global_variable.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

#include <iomanip>
#include <iostream>

#include "module_base/tool_quit.h"
namespace elecstate
{

/// @brief calculate etot
double fenergy::calculate_etot()
{
    if(GlobalV::use_paw)
    {
        etot = eband + deband + etxc + ewald_energy - hartree_energy + demet + descf + exx + efield + gatefield
               + evdw + esol_el + esol_cav + edftu + edeepks_scf;
    }
    else
    {
        etot = eband + deband + (etxc - etxcc) + ewald_energy + hartree_energy + demet + descf + exx + efield + gatefield
               + evdw + esol_el + esol_cav + edftu + edeepks_scf;
    }

#ifdef USE_PAW
    if(GlobalV::use_paw)
    {
        double ecore = GlobalC::paw_cell.calculate_ecore();
        double epawdc = GlobalC::paw_cell.get_epawdc();
        etot += ( ecore + epawdc );
    }
#endif
    return etot;
}

/// @brief calculate etot_harris
double fenergy::calculate_harris()
{
    if(GlobalV::use_paw)
    {
        etot_harris = eband + deband_harris + etxc + ewald_energy - hartree_energy + demet + descf + exx + efield
                  + gatefield + evdw + esol_el + esol_cav + edftu + edeepks_scf;
    }
    else
    {
        etot_harris = eband + deband_harris + (etxc - etxcc) + ewald_energy + hartree_energy + demet + descf + exx + efield
                  + gatefield + evdw + esol_el + esol_cav + edftu + edeepks_scf;
    }
#ifdef USE_PAW
    if(GlobalV::use_paw)
    {
        double ecore = GlobalC::paw_cell.calculate_ecore();
        double epawdc = GlobalC::paw_cell.get_epawdc();
        etot_harris += ( ecore + epawdc );
    }
#endif
    return etot_harris;
}

/// @brief set all energies to zero
void fenergy::clear_all()
{
    etot = etot_old = eband = deband = etxc = etxcc = vtxc = ewald_energy = hartree_energy = demet = descf = exx
        = efield = gatefield = evdw = etot_harris = deband_harris = esol_el = esol_cav = edftu = edeepks_scf = 0.0;
}

/// @brief print all energies
void fenergy::print_all() const
{
    std::cout << std::resetiosflags(std::ios::scientific) << std::endl;
    std::cout << std::setprecision(16) << std::endl;
    std::cout << " eband=" << eband << std::endl;
    std::cout << " deband=" << deband << std::endl;
    std::cout << " etxc-etxcc=" << etxc - etxcc << std::endl;
    std::cout << " ewld=" << ewald_energy << std::endl;
    std::cout << " ehart=" << hartree_energy << std::endl;
    std::cout << " entropy(-TS)=" << demet << std::endl;
    std::cout << " descf=" << descf << std::endl;
    std::cout << " exx=" << exx << std::endl;
    std::cout << " efiled=" << efield << std::endl;
    std::cout << " gatefiled=" << gatefield << std::endl;
    std::cout << " evdw=" << evdw << std::endl;
    std::cout << " esol_el=" << esol_el << std::endl;
    std::cout << " esol_cav=" << esol_cav << std::endl;
    std::cout << " edftu=" << edftu << std::endl;
    std::cout << " edeepks_scf=" << edeepks_scf << std::endl;
    std::cout << std::endl;
    std::cout << " total= " << etot << std::endl;
}

/// @brief get the reference of fermi of a specific spin
/// @param is SPIN
/// @return a reference of fermi(is)
double& efermi::get_ef(const int& is)
{
    if (!two_efermi)
    {
        return this->ef;
    }
    else if (is == 0)
    {
        return this->ef_up;
    }
    else if (is == 1)
    {
        return this->ef_dw;
    }
    else
    {
        ModuleBase::WARNING_QUIT("energy", "Please check NSPIN when TWO_EFERMI is true");
        __builtin_unreachable();
    }
}

/// @brief get the value of fermi of a specific spin
/// @param is SPIN
/// @return value of fermi(is)
double efermi::get_efval(const int& is) const
{
    if (!two_efermi)
    {
        return this->ef;
    }
    else if (is == 0)
    {
        return this->ef_up;
    }
    else if (is == 1)
    {
        return this->ef_dw;
    }
    else
    {
        ModuleBase::WARNING_QUIT("energy", "Please check NSPIN when TWO_EFERMI is true");
        __builtin_unreachable();
    }
}

} // namespace elecstate

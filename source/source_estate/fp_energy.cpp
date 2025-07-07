#include "fp_energy.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_variable.h"


#include "source_base/tool_quit.h"

#include <iomanip>
#include <iostream>

namespace elecstate
{

/// @brief calculate etot
double fenergy::calculate_etot()
{
    etot = eband + deband + (etxc - etxcc) + ewald_energy + hartree_energy + demet + descf + exx + efield
            + gatefield + evdw + esol_el + esol_cav + edftu + edeepks_scf + escon;
    return etot;
}

/// @brief calculate etot_harris
double fenergy::calculate_harris()
{
    etot_harris = eband + deband_harris + (etxc - etxcc) + ewald_energy + hartree_energy + demet + descf + exx
                    + efield + gatefield + evdw + esol_el + esol_cav + edftu + edeepks_scf + escon;
    return etot_harris;
}

/// @brief set all energies to zero
void fenergy::clear_all()
{
    etot = etot_old = eband = deband = etxc = etxcc = vtxc = ewald_energy = hartree_energy = demet = descf = exx
        = efield = gatefield = evdw = etot_harris = deband_harris = esol_el = esol_cav = edftu = edeepks_scf = escon
        = 0.0;
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
    std::cout << " escon=" << escon << std::endl;
    std::cout << std::endl;
    std::cout << " total= " << etot << std::endl;
}

/// @brief set efermi of a specific spin
/// @param is SPIN
/// @param ef_in fermi(is)
void efermi::set_efval(const int& is, const double& ef_in)
{
    if (!two_efermi)
    {
        this->ef = ef_in;
    }
    else if (is == 0)
    {
        this->ef_up = ef_in;
    }
    else if (is == 1)
    {
        this->ef_dw = ef_in;
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

/// @brief get all fermi energies for all spins
/// @return all fermi energies for all spins
std::vector<double> efermi::get_all_ef() const
{
    if (two_efermi)
    {
        return {ef_up, ef_dw};
    }
    else
    {
        return {ef, ef}; // For NSPIN=1, ef_up=ef_dw=ef
    }
}

} // namespace elecstate

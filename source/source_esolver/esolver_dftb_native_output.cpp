#include "esolver_dftb_native.h"

#include "source_base/parallel_common.h"
#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace ModuleESolver
{
namespace
{
constexpr double dftb_hartree_to_ev = 27.211386245988;
} // namespace

void ESolver_DFTBNative::write_result_outputs(const UnitCell& ucell,
                                               const ModuleDFTB::DftbPeriodicInput& input,
                                               const std::string& output_dir,
                                               std::ostream& running_log,
                                               std::ostream& dftb_log) const
{
    const double total_free_energy_ev = this->result_.total_free_energy_hartree * dftb_hartree_to_ev;
    running_log << std::setprecision(16)
                         << "\n Native DFTB SCC converged: " << std::boolalpha << this->result_.converged
                         << ", iterations=" << this->result_.scc_iterations
                         << ", max |delta q|=" << this->result_.maximum_charge_residual << " e\n"
                         << " Native DFTB Fermi level: " << this->result_.fermi_energy_hartree << " Ha\n"
                         << " Native DFTB band free energy: " << this->result_.band_free_energy_hartree << " Ha\n"
                         << " Native DFTB H0 energy: " << this->result_.h0_energy_hartree << " Ha\n"
                         << " Native DFTB SCC energy: " << this->result_.scc_energy_hartree << " Ha\n"
                         << " Native DFTB third-order energy: " << this->result_.third_order_energy_hartree << " Ha\n"
                         << " Native DFTB repulsive energy: " << this->result_.repulsive_energy_hartree << " Ha\n"
                         << " Native DFTB total free energy: " << this->result_.total_free_energy_hartree << " Ha\n"
                         << " Native DFTB total free energy: " << total_free_energy_ev << " eV\n"
                         << " Native DFTB component sum: "
                         << this->result_.h0_energy_hartree + this->result_.scc_energy_hartree
                            + this->result_.third_order_energy_hartree + this->result_.repulsive_energy_hartree
                         << " Ha\n";
    dftb_log << "# SCC converged: " << std::boolalpha << this->result_.converged
             << ", iterations: " << this->result_.scc_iterations
             << ", maximum charge residual: " << this->result_.maximum_charge_residual << " e\n"
             << "# Fermi energy: " << this->result_.fermi_energy_hartree << " Ha\n"
             << "# Band energy: " << this->result_.band_energy_hartree << " Ha\n"
             << "# Band free energy: " << this->result_.band_free_energy_hartree << " Ha\n"
             << "# H0 energy: " << this->result_.h0_energy_hartree << " Ha\n"
             << "# SCC energy: " << this->result_.scc_energy_hartree << " Ha\n"
             << "# Third-order energy: " << this->result_.third_order_energy_hartree << " Ha\n"
             << "# Repulsive energy: " << this->result_.repulsive_energy_hartree << " Ha\n"
             << "# Total free energy: " << this->result_.total_free_energy_hartree << " Ha / "
             << total_free_energy_ev << " eV\n";
    double total_electron_excess = 0.0;
    std::size_t atom_index = 0;
    running_log << " Native DFTB Mulliken populations and electron excess charges:\n"
                << "  atom  element  population(e)  electron_excess(e)\n";
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const auto& species = ucell.atoms[it];
        double species_excess = 0.0;
        for (int ia = 0; ia < species.na; ++ia, ++atom_index)
        {
            const double charge = this->result_.electron_excess_charges[atom_index];
            const double population = this->result_.atomic_electron_populations[atom_index];
            species_excess += charge;
            total_electron_excess += charge;
            running_log << "  " << std::setw(4) << atom_index + 1 << "  " << std::setw(7) << species.label
                        << "  " << std::setw(18) << population << "  " << charge << "\n";
        }
        if (species.na > 0)
            running_log << " Native DFTB mean electron excess " << species.label << ": "
                        << species_excess / species.na << " e\n";
    }
    running_log << " Native DFTB net electron excess: " << total_electron_excess << " e\n"
                << " #TOTAL ENERGY# " << total_free_energy_ev << " eV (native DFTB)\n";
    dftb_log << "# Net electron excess: " << total_electron_excess << " e\n"
             << "# Atomic Mulliken populations and excess charges\n"
             << "# atom element population(e) electron_excess(e)\n";
    atom_index = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const auto& species = ucell.atoms[it];
        for (int ia = 0; ia < species.na; ++ia, ++atom_index)
            dftb_log << atom_index + 1 << ' ' << species.label << ' '
                     << this->result_.atomic_electron_populations[atom_index] << ' '
                     << this->result_.electron_excess_charges[atom_index] << '\n';
    }
    dftb_log.flush();

    const std::string eigen_file = output_dir + "eig_occ.txt";
    std::ofstream eigenvalues(eigen_file.c_str());
    if (!eigenvalues)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB eigenvalue file: " + eigen_file);
    eigenvalues << std::setprecision(this->output_precision_)
                << "# K-point eigenvalues and occupations from the converged native DFTB SCC solution\n"
                << "# Module: Native DFTB eigensolver\n"
                << "# Units: energy in eV; occupations are electrons per spin-degenerate state; k-point weights sum to 1\n"
                << "# k_index weight kx(Direct) ky(Direct) kz(Direct) band_index energy(eV) occupation(e)\n";
    for (std::size_t ik = 0; ik < this->result_.kpoint_eigenvalues.size(); ++ik)
    {
        const auto& point = this->result_.kpoint_eigenvalues[ik];
        for (std::size_t band = 0; band < point.eigenvalues_hartree.size(); ++band)
            eigenvalues << ik + 1 << ' ' << point.weight << ' '
                        << point.fractional[0] << ' ' << point.fractional[1] << ' ' << point.fractional[2] << ' '
                        << band + 1 << ' ' << point.eigenvalues_hartree[band] * dftb_hartree_to_ev << ' '
                        << point.occupations[band] << '\n';
    }
    eigenvalues.close();
    running_log << " Native DFTB eigenvalues and occupations: wrote " << eigen_file << "\n";

    const std::string mulliken_file = output_dir + "mulliken.txt";
    std::ofstream mulliken(mulliken_file.c_str());
    if (!mulliken)
        ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB Mulliken file: " + mulliken_file);
    mulliken << std::setprecision(this->output_precision_)
             << "# Atomic Mulliken populations from the converged native DFTB SCC solution\n"
             << "# Module: Native DFTB Mulliken analysis\n"
             << "# Units: population and electron excess in e; excess = population - SKF neutral valence\n"
             << "# atom element population(e) electron_excess(e)\n";
    atom_index = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const auto& species = ucell.atoms[it];
        for (int ia = 0; ia < species.na; ++ia, ++atom_index)
            mulliken << atom_index + 1 << ' ' << species.label << ' '
                     << this->result_.atomic_electron_populations[atom_index] << ' '
                     << this->result_.electron_excess_charges[atom_index] << '\n';
    }
    mulliken.close();
    running_log << " Native DFTB Mulliken populations: wrote " << mulliken_file << "\n";

    if (!this->result_.band_structure.empty())
    {
        const std::string band_file = output_dir + "band.txt";
        std::ofstream bands(band_file.c_str());
        if (!bands)
            ModuleBase::WARNING_QUIT("ESolver_DFTBNative", "Cannot write native DFTB band file: " + band_file);
        bands << std::setprecision(this->output_precision_)
              << "# Band energies from the converged native DFTB SCC potential (frozen-potential path solve)\n"
              << "# Module: Native DFTB band structure\n"
              << "# Units: k_distance in bohr^-1, eigenvalues and E-Ef in eV; coordinates are fractional reciprocal\n"
              << "# Fermi energy: " << this->result_.fermi_energy_hartree * dftb_hartree_to_ev << " eV\n"
              << "# k_index k_distance(bohr^-1) kx ky kz label band_index energy(eV) E-Ef(eV)\n";
        for (std::size_t ik = 0; ik < this->result_.band_structure.size(); ++ik)
        {
            const auto& point = this->result_.band_structure[ik];
            for (std::size_t band = 0; band < point.eigenvalues_hartree.size(); ++band)
            {
                const double energy = point.eigenvalues_hartree[band];
                bands << ik + 1 << ' ' << point.distance_inverse_bohr << ' '
                      << point.fractional[0] << ' ' << point.fractional[1] << ' ' << point.fractional[2] << ' '
                      << (point.label.empty() ? "-" : point.label) << ' ' << band + 1 << ' '
                      << energy * dftb_hartree_to_ev << ' '
                      << (energy - this->result_.fermi_energy_hartree) * dftb_hartree_to_ev << '\n';
            }
        }
        bands.close();
        running_log << " Native DFTB frozen-potential band structure: "
                    << this->result_.band_structure.size() << " k-points, "
                    << (this->result_.band_structure.front().eigenvalues_hartree.size())
                    << " bands; wrote " << band_file << "\n";
    }
    dftb_log << "#TOTAL ENERGY# " << total_free_energy_ev << " eV (native DFTB)\n"
             << "!FINAL_ETOT_IS " << total_free_energy_ev << " eV (native DFTB)\n"
             << "# Output files: dftb.log, eig_occ.txt, mulliken.txt"
             << (this->result_.band_structure.empty() ? "\n" : ", band.txt\n");
}

void ESolver_DFTBNative::after_all_runners(BaseCell& cell)
{
    static_cast<void>(cell);
    if (Parallel_Common::get_rank() != 0) return;
    std::cout << std::setprecision(16)
              << "\n --------------------------------------------" << std::endl
              << " !FINAL_ETOT_IS " << this->result_.total_free_energy_hartree * dftb_hartree_to_ev
              << " eV (native DFTB)" << std::endl;
}
} // namespace ModuleESolver

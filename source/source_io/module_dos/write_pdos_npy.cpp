#include "write_pdos_npy.h"
#include "pdos_orbitals.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/npy_record.h"

#include <algorithm>
#include <locale>
#include <ostream>
#include <sstream>

std::string ModuleIO::write_pdos_npy(const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                                     int nspin, double emin, double de, const std::string& out_dir,
                                     int istep, int out_app_flag)
{
    const int step = istep < 0 ? 1 : istep + 1;
    const bool numbered = out_app_flag == 2 && istep >= 0;
    const std::string suffix = numbered ? "g" + std::to_string(step) : "";
    const std::string path = out_dir + "pdos" + suffix + ".npy";
    const int ne = pdos.front().nc;
    std::vector<std::vector<ModuleIO::PdosOrbital>> by_type(cell.ntype);
    size_t width = 1;
    for (int it = 0; it < cell.ntype; ++it)
    {
        by_type[it] = ModuleIO::pdos_orbitals(cell.atoms[it]);
        width = std::max(width, cell.atoms[it].label.size());
    }
    width = (width + 7) / 8 * 8; // keep the float64 spectrum aligned in every row
    size_t norb = 0;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        norb += by_type[cell.iat2it[iat]].size();
    }

    std::ostringstream fields;
    fields.imbue(std::locale::classic());
    fields << "('nspin','=i4'),('atom_index','=i4'),('l','=i4'),('m','=i4'),"
           << "('energy_min_eV','=f8'),('energy_step_eV','=f8'),('species','|S" << width
           << "'),('pdos','=f8',(" << pdos.size() << "," << ne << "))";
    const uint64_t row_bytes = 16 + 4 * 4 + 2 * 8 + width + 8 * static_cast<uint64_t>(ne) * pdos.size();
    const bool append = out_app_flag == 1 && istep > 0;
    const double grid[] = {emin, de};
    const bool spinor = nspin == 4;
    const std::string descr = fields.str();
    write_npy_record(path, append, step, descr, "pdos", norb, row_bytes, [&](std::ostream& out, const NpyPrefix& prefix) {
        std::vector<double> spectrum(ne);
        for (int iat = 0; iat < cell.nat; ++iat)
        {
            const int it = cell.iat2it[iat];
            std::string label = cell.atoms[it].label;
            label.resize(width, '\0');
            const int first = cell.get_iat2iwt()[iat];
            for (const auto& orbital : by_type[it])
            {
                const int32_t attributes[] = {nspin, iat + 1, orbital.l, orbital.m};
                out.write(prefix.data(), prefix.size());
                out.write(reinterpret_cast<const char*>(attributes), sizeof(attributes));
                out.write(reinterpret_cast<const char*>(grid), sizeof(grid));
                out.write(label.data(), width);
                for (const auto& spin : pdos)
                {
                    for (int n = 0; n < ne; ++n)
                    {
                        spectrum[n] = orbital.value(spin, first, n, spinor);
                    }
                    out.write(reinterpret_cast<const char*>(spectrum.data()), spectrum.size() * sizeof(double));
                }
            }
        }
    });
    return path;
}

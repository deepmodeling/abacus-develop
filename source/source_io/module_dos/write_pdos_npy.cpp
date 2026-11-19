#include "write_pdos_npy.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/npy_record.h"

#include <algorithm>
#include <cmath>
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
    size_t width = 1;
    size_t norb = 0;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        const Atom& atom = cell.atoms[cell.iat2it[iat]];
        width = std::max(width, atom.label.size());
        norb += std::count(atom.iw2n.begin(), atom.iw2n.end(), 0);
    }
    width = (width + 7) / 8 * 8; // align the binary64 spectrum in every row

    std::ostringstream fields;
    fields.imbue(std::locale::classic());
    fields << "('nspin','=i4'),('atom_index','=i4'),('l','=i4'),('m','=i4'),"
           << "('energy_min_eV','=f8'),('energy_step_eV','=f8'),('species','|S" << width
           << "'),('pdos','=f8',(" << pdos.size() << "," << ne << "))";
    const uint64_t row_bytes = 16 + 4 * 4 + 2 * 8 + width + 8 * static_cast<uint64_t>(ne) * pdos.size();
    const bool append = out_app_flag == 1 && istep > 0;
    const double grid[] = {emin, de};
    const int npol = nspin == 4 ? 2 : 1;
    const std::string descr = fields.str();
    write_npy_record(path, append, step, descr, "pdos", norb, row_bytes, [&](std::ostream& out, const NpyPrefix& prefix) {
        std::vector<double> spectrum(ne);
        for (int iat = 0; iat < cell.nat; ++iat)
        {
            const Atom& atom = cell.atoms[cell.iat2it[iat]];
            std::string label = atom.label;
            label.resize(width, '\0');
            const int first = cell.get_iat2iwt()[iat];
            // Atom::set_index orders orbitals by l, then zeta, then real-harmonic column.
            for (int iw = 0; iw < atom.nw; ++iw)
            {
                if (atom.iw2n[iw] != 0) { continue; }
                const int l = atom.iw2l[iw];
                const int col = atom.iw2m[iw];
                const int m = col % 2 == 0 ? -col / 2 : (col + 1) / 2;
                const int32_t attributes[] = {nspin, iat + 1, l, m};
                out.write(prefix.data(), prefix.size());
                out.write(reinterpret_cast<const char*>(attributes), sizeof(attributes));
                out.write(reinterpret_cast<const char*>(grid), sizeof(grid));
                out.write(label.data(), width);
                for (const auto& spin : pdos)
                {
                    std::fill(spectrum.begin(), spectrum.end(), 0.0);
                    for (int zeta = 0; zeta < atom.l_nchi[l]; ++zeta)
                    {
                        const int row = first + npol * (iw + zeta * (2 * l + 1));
                        for (int n = 0; n < ne; ++n)
                        {
                            spectrum[n] += npol == 2 ? spin(row, n) + spin(row + 1, n) : spin(row, n);
                        }
                    }
                    for (double& value : spectrum)
                    {
                        if (std::abs(value) < 1e-6) { value = 0.0; }
                    }
                    out.write(reinterpret_cast<const char*>(spectrum.data()), spectrum.size() * sizeof(double));
                }
            }
        }
    });
    return path;
}

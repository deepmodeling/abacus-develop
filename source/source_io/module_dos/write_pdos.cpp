#include "write_pdos.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/npy_record.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>

ModuleIO::PdosBlock ModuleIO::atom_pdos(const UnitCell& cell, int iat, const ModuleBase::matrix* pdos,
                                        int spins, bool spinor)
{
    const Atom& atom = cell.atoms[cell.iat2it[iat]];
    const int ne = pdos[0].nc;
    const int first = cell.get_iat2iwt()[iat];
    const int npol = spinor ? 2 : 1;
    PdosBlock block;
    for (int iw = 0; iw < atom.nw; ++iw)
    {
        auto& channel = block[{atom.iw2l[iw], atom.iw2m[iw]}];
        if (channel.empty()) { channel.resize(static_cast<size_t>(spins) * ne); }
        const int row = first + npol * iw;
        for (int is = 0; is < spins; ++is)
        {
            double* values = channel.data() + static_cast<size_t>(is) * ne;
            for (int n = 0; n < ne; ++n)
            {
                values[n] += spinor ? pdos[is](row, n) + pdos[is](row + 1, n) : pdos[is](row, n);
            }
        }
    }
    for (auto& channel : block)
    {
        for (double& value : channel.second)
        {
            if (std::abs(value) < 1e-6) { value = 0.0; }
        }
    }
    return block;
}

std::string ModuleIO::write_pdos_npy(const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                                     int nspin, double emin, double de, const std::string& out_dir,
                                     int istep, int out_app_flag)
{
    const int step = std::max(1, istep + 1);
    const bool numbered = out_app_flag == 2 && istep >= 0;
    const std::string suffix = numbered ? "g" + std::to_string(step) : "";
    const std::string path = out_dir + "pdos" + suffix + ".npy";
    const int ne = pdos.front().nc;
    size_t width = 1;
    size_t count = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        const Atom& atom = cell.atoms[it];
        width = std::max(width, atom.label.size());
        for (int l = 0; l <= atom.nwl; ++l)
        {
            if (atom.l_nchi[l] > 0) { count += static_cast<size_t>(atom.na) * (2 * l + 1); }
        }
    }
    width = (width + 7) / 8 * 8;
    std::ostringstream fields;
    fields.imbue(std::locale::classic());
    fields << "('nspin','=i4'),('atom_index','=i4'),('l','=i4'),('m','=i4'),"
           << "('energy_min_eV','=f8'),('energy_step_eV','=f8'),('species','|S" << width
           << "'),('pdos','=f8',(" << pdos.size() << "," << ne << "))";
    const std::uint64_t row_bytes = 48 + width + sizeof(double) * static_cast<std::uint64_t>(ne) * pdos.size();
    const bool append = out_app_flag == 1 && istep > 0;
    const std::string descr = fields.str();
    std::ofstream out;
    const auto prefix = open_npy_record(out, path, append, step, descr, "pdos", count, row_bytes);
    const double grid[] = {emin, de};
    const bool spinor = nspin == 4;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        std::string label = cell.atoms[cell.iat2it[iat]].label;
        label.resize(width, '\0');
        const auto channels = atom_pdos(cell, iat, pdos.data(), pdos.size(), spinor);
        for (const auto& channel : channels)
        {
            const int col = channel.first.second;
            const int m = col % 2 == 0 ? -col / 2 : (col + 1) / 2;
            const std::int32_t attrs[] = {nspin, iat + 1, channel.first.first, m};
            out.write(prefix.data(), prefix.size());
            out.write(reinterpret_cast<const char*>(attrs), sizeof(attrs));
            out.write(reinterpret_cast<const char*>(grid), sizeof(grid));
            out.write(label.data(), label.size());
            out.write(reinterpret_cast<const char*>(channel.second.data()), channel.second.size() * sizeof(double));
        }
    }
    out.close();
    return path;
}

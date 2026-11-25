#include "write_pdos_text.h"
#include "dos_npy_header.h"

#include "source_base/global_variable.h"
#include "source_io/module_output/provenance_header.h"
#include "source_io/module_output/write_notice.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>

namespace
{
// Each (l, real-harmonic column) contains the zeta sum in spin-major order.
using PdosBlock = std::map<std::pair<int, int>, std::vector<double>>;
PdosBlock atom_pdos(const UnitCell& cell, int iat, const ModuleBase::matrix* pdos,
                                        int spins, bool spinor)
{
    const Atom& atom = cell.atoms[cell.iat2it[iat]];
    const int ne = pdos[0].nc;
    const int first = cell.get_iat2iwt()[iat];
    const int npol = spinor ? 2 : 1;
    PdosBlock block;
    for (int iw = 0; iw < atom.nw; ++iw)
    {
        std::vector<double>& channel = block[{atom.iw2l[iw], atom.iw2m[iw]}];
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
    for (PdosBlock::value_type& channel : block)
    {
        for (double& value : channel.second)
        {
            if (std::abs(value) < 1e-6) { value = 0.0; }
        }
    }
    return block;
}
}

void ModuleIO::write_pdos_text(
        const UnitCell& ucell,
        const ModuleBase::matrix* pdos,
        const int nspin,
        const int nlocal,
        const int npoints,
        const double& emin,
        const double& dos_edelta_ev,
        const std::string& out_dir,
        const int istep,
        const int out_app_flag,
        const int precision,
        std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_pdos_text");

    // number of spin channels to write
    const int nspin0 = (nspin == 2) ? 2 : 1;

    // Modes 0 and 1 reuse the fixed filename (mode 0 overwrites, mode 1 appends
    // each ionic step); only mode 2 embeds the step number in the filename.
    // Single-point SCF passes istep = -1 and truncates.
    // The first written step (istep == 0) truncates so a fresh run does not
    // keep records from a previous run; only later steps (istep > 0) append.
    const int istep_file = (out_app_flag == 2) ? istep : -1;
    const bool append_step = (out_app_flag == 1) && (istep > 0);
    const std::string time_stamp = ModuleIO::local_timestamp();
    const int istep_one_based = istep + 1;

    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;
        ss << out_dir << "pdoss" << is + 1;
        if (istep_file >= 0)
        {
            ss << "g" << istep_file + 1;
        }
        ss << ".txt";

        std::ios_base::openmode mode = std::ios::out;
        if (append_step)
        {
            mode |= std::ios::app;
        }
        std::ofstream ofs(ss.str().c_str(), mode);

        // iw2m stores m as 0..2l, mapping to physical m values:
        // 0->0, 1->+1, 2->-1, 3->+2, 4->-2, 5->+3, 6->-3
        // So columns are ordered: s, p(m=0,+1,-1), d(m=0,+1,-1,+2,-2), f(...)
        ofs << ModuleIO::prov_step(istep_one_based, time_stamp);
        ofs << "# npoints: " << npoints << std::endl;
        ofs << "# energy(eV)  atom  species  pdos(1/eV), columns: s(m=0) p(m=0,+1,-1) d(m=0,+1,-1,+2,-2) f(m=0,+1,-1,+2,-2,+3,-3)" << std::endl;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int it = ucell.iat2it[iat];
            const bool spinor = nspin == 4;
            const ModuleBase::matrix* spin_pdos = pdos + is;
            const PdosBlock channels = atom_pdos(ucell, iat, spin_pdos, 1, spinor);

            for (int n = 0; n < npoints; ++n)
            {
                const double en = emin + n * dos_edelta_ev;

                ofs << std::setw(12) << std::fixed << std::setprecision(precision) << en
                    << "  " << std::setw(3) << iat + 1
                    << "  " << std::setw(4) << ucell.atoms[it].label;

                for (const PdosBlock::value_type& channel : channels)
                {
                    ofs << "  " << std::setw(8) << std::fixed << std::setprecision(precision) << channel.second[n];
                }
                ofs << std::endl;
            }
        }

        ofs.close();

        ModuleIO::write_notice(ofs_running, "PDOS", is, nspin, ss.str());
    }
}

void ModuleIO::write_pdos_npy(const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                           int nspin, double emin, double de, const std::string& path)
{
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
    const std::uint64_t row_bytes = 32 + width + sizeof(double) * static_cast<std::uint64_t>(ne) * pdos.size();
    if (row_bytes > std::numeric_limits<std::int32_t>::max())
    {
        throw std::runtime_error("PDOS array exceeds NumPy v1 dtype limits");
    }
    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.open(path, std::ios::binary | std::ios::trunc);
    write_dos_npy_header(out, fields.str(), count);
    const double grid[] = {emin, de};
    const bool spinor = nspin == 4;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        std::string label = cell.atoms[cell.iat2it[iat]].label;
        label.resize(width, '\0');
        const PdosBlock channels = atom_pdos(cell, iat, pdos.data(), pdos.size(), spinor);
        for (const PdosBlock::value_type& channel : channels)
        {
            const int col = channel.first.second;
            const int m = col % 2 == 0 ? -col / 2 : (col + 1) / 2;
            const std::int32_t attrs[] = {nspin, iat + 1, channel.first.first, m};
            out.write(reinterpret_cast<const char*>(attrs), sizeof(attrs));
            out.write(reinterpret_cast<const char*>(grid), sizeof(grid));
            out.write(label.data(), label.size());
            const std::size_t spectrum_bytes = channel.second.size() * sizeof(double);
            out.write(reinterpret_cast<const char*>(channel.second.data()), spectrum_bytes);
        }
    }
    out.close();
}

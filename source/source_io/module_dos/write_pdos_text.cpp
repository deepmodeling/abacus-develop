#include "write_pdos_text.h"
#include "write_pdos.h"

#include "source_base/global_variable.h"

#include <fstream>
#include <iomanip>
#include <sstream>

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
        const int out_app_flag)
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
        ofs << "# istep: " << istep + 1 << std::endl;
        ofs << "# npoints: " << npoints << std::endl;
        ofs << "# energy(eV)  atom  species  pdos(1/eV), columns: s(m=0) p(m=0,+1,-1) d(m=0,+1,-1,+2,-2) f(m=0,+1,-1,+2,-2,+3,-3)" << std::endl;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int it = ucell.iat2it[iat];
            const bool spinor = nspin == 4;
            const auto channels = atom_pdos(ucell, iat, pdos + is, 1, spinor);

            for (int n = 0; n < npoints; ++n)
            {
                const double en = emin + n * dos_edelta_ev;

                ofs << std::setw(12) << std::fixed << std::setprecision(6) << en
                    << "  " << std::setw(3) << iat + 1
                    << "  " << std::setw(4) << ucell.atoms[it].label;

                for (const auto& channel : channels)
                {
                    ofs << "  " << std::setw(8) << std::fixed << std::setprecision(6) << channel.second[n];
                }
                ofs << std::endl;
            }
        }

        ofs.close();
    }
}

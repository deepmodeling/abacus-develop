#include "relax_stru_io.h"

#include "source_base/formatter.h"
#include "source_cell/cif_io.h"
#include "source_cell/print_cell.h"
#include "source_main/version.h"

#include <ctime>

namespace relax_stru_io
{

std::string build_stru_header(const int istep,
                              const double etot,
                              const ModuleBase::matrix& stress,
                              const bool is_final)
{
    std::time_t now = std::time(nullptr);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    const char* step_label = is_final ? "# RELAX STEP %d (FINAL), Energy: %.8f eV\n"
                                      : "# RELAX STEP %d, Energy: %.8f eV\n";
    std::string header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n",
                                         VERSION,
                                         time_buf);
    header += FmtCore::format(step_label, istep + 1, etot * ModuleBase::Ry_to_eV);

    // stress in kbar: Ry/Bohr^3 -> kbar, 3 rows
    const double stress_transform = ModuleBase::RYDBERG_SI
                                    / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI
                                       * ModuleBase::BOHR_RADIUS_SI)
                                    * 1.0e-8;
    for (int i = 0; i < 3; i++)
    {
        header += FmtCore::format("# Stress (kbar): %.6f %.6f %.6f\n",
                                  stress(i, 0) * stress_transform,
                                  stress(i, 1) * stress_transform,
                                  stress(i, 2) * stress_transform);
    }
    return header;
}

bool need_orbital(const Input_para& inp)
{
    bool need_orb = inp.basis_type == "pw";
    need_orb = need_orb && inp.init_wfc.substr(0, 3) == "nao";
    need_orb = need_orb || inp.basis_type == "lcao";
    need_orb = need_orb || inp.basis_type == "lcao_in_pw";
    return need_orb;
}

void write_stru(UnitCell& ucell,
                const Input_para& inp,
                const std::string& filename,
                const std::string& header,
                const ModuleBase::matrix& force,
                const bool need_orb,
                const bool deepks_setorb,
                const int my_rank)
{
    if (inp.out_stru == 1)
    {
        unitcell::print_stru_file(ucell,
                                  ucell.atoms,
                                  ucell.latvec,
                                  filename,
                                  header,
                                  inp.nspin,
                                  true,
                                  inp.calculation == "md",
                                  inp.out_mul,
                                  need_orb,
                                  deepks_setorb,
                                  my_rank,
                                  force);
    }
    else if (inp.out_stru == 2)
    {
        ModuleIO::CifParser::write(filename, ucell, header, "data_?", my_rank);
    }
}

} // namespace relax_stru_io

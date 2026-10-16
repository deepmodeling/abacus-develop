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
                              const Input_para& inp,
                              const bool is_final,
                              const bool geometry_evaluated)
{
    std::time_t now = std::time(nullptr);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    const char* final_tag = is_final ? " (FINAL)" : "";
    const double etot_ev = etot * ModuleBase::Ry_to_eV;
    // Line 3 is labeled with an energy value only when that energy belongs
    // to the geometry written below. When the optimizer proposed the final
    // geometry but it was never evaluated, etot is stale, so print N/A to
    // keep the energy label and the coordinates on the same frame.
    std::string header;
    if (geometry_evaluated)
    {
        header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n# RELAX STEP %d%s, Energy: %.8f eV\n",
                                 VERSION,
                                 time_buf,
                                 istep + 1,
                                 final_tag,
                                 etot_ev);
    }
    else
    {
        header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n# RELAX STEP %d%s, Energy: N/A\n",
                                 VERSION,
                                 time_buf,
                                 istep + 1,
                                 final_tag);
    }

    // stress in kbar: Ry/Bohr^3 -> kbar, always exactly 3 lines for
    // downstream parsers. N/A marks uncomputed stress; no inline comment is
    // appended so that every stress line has the same token count.
    const double stress_transform = ModuleBase::RYDBERG_SI
                                    / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI
                                       * ModuleBase::BOHR_RADIUS_SI)
                                    * 1.0e-8;
    const bool stress_valid = inp.cal_stress && geometry_evaluated;
    if (stress_valid)
    {
        for (int i = 0; i < 3; i++)
        {
            header += FmtCore::format("# Stress (kbar): %.6f %.6f %.6f\n",
                                      stress(i, 0) * stress_transform,
                                      stress(i, 1) * stress_transform,
                                      stress(i, 2) * stress_transform);
        }
    }
    else
    {
        header += "# Stress (kbar): N/A N/A N/A\n"
                  "# Stress (kbar): N/A N/A N/A\n"
                  "# Stress (kbar): N/A N/A N/A\n";
    }

    // Line 7: single NOTE line describing the state of stress/forces/geometry.
    // Keeping the header at exactly 7 lines makes it easy for downstream
    // parsers to skip a fixed-size comment block.
    std::string note;
    if (!geometry_evaluated)
    {
        note = "# NOTE: geometry proposed by optimizer but not evaluated; energy N/A; stress N/A; forces omitted";
    }
    else if (!inp.cal_stress && !inp.cal_force)
    {
        note = "# NOTE: stress not computed (cal_stress=0); forces not computed (cal_force=0); per-atom f fields omitted intentionally";
    }
    else if (!inp.cal_stress)
    {
        note = "# NOTE: stress not computed (cal_stress=0)";
    }
    else if (!inp.cal_force)
    {
        note = "# NOTE: forces not computed (cal_force=0); per-atom f fields omitted intentionally";
    }
    else
    {
        note = "# NOTE: stress and forces computed for this geometry";
    }
    header += note + "\n";

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
                const int my_rank,
                const bool has_force)
{
    if (inp.out_stru == 1)
    {
        unitcell::print_stru_file(ucell,
                                  ucell.atoms,
                                  ucell.latvec,
                                  filename,
                                  header,
                                  inp.nspin,
                                  inp.calculation == "md",
                                  inp.out_mul,
                                  need_orb,
                                  deepks_setorb,
                                  my_rank,
                                  force,
                                  has_force);
    }
    else if (inp.out_stru == 2)
    {
        ModuleIO::CifParser::write(filename, ucell, header, "data_?", my_rank);
    }
}

} // namespace relax_stru_io

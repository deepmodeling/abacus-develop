#include "relax_driver.h"

#include "source_base/formatter.h"
#include "source_main/version.h"

#include <ctime>

std::string Relax_Driver::build_stru_header(const int istep,
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
    std::string header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n# RELAX STEP %d%s, Energy: %.8f eV\n",
                                          VERSION,
                                          time_buf,
                                          istep + 1,
                                          final_tag,
                                          etot * ModuleBase::Ry_to_eV);

    // stress in kbar: Ry/Bohr^3 -> kbar, always exactly 3 lines for
    // downstream parsers. N/A marks uncomputed stress; no inline comment is
    // appended so that every stress line has the same token count.
    const double stress_transform = ModuleBase::RYDBERG_SI
                                    / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI)
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
        note = "# NOTE: geometry proposed by optimizer but not evaluated; stress N/A; forces omitted; energy above belongs to the last evaluated geometry";
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

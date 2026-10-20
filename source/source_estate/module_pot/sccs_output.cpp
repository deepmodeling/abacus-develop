#include "sccs_output.h"
#include "source_base/constants.h"
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
// The report is a grid: every cell (tag, key or value) has one width and is
// left-aligned, so every cell starts on a column edge. Energies keep 9
// significant digits, other reals 8.
const int cell_width = 16;
const int energy_digits = 9;
const int value_digits = 8;

std::string real_text(double value, int digits)
{
    std::ostringstream text;
    const int mantissa_digits = digits - 1;
    text << std::scientific << std::setprecision(mantissa_digits) << value;
    return text.str();
}
std::string int_text(int value)
{
    return std::to_string(value);
}
// One line of cells, without trailing blanks. A cell wider than the column
// keeps one blank after it.
void write_row(std::ostream& output, const std::vector<std::string>& cells)
{
    const std::size_t width = cell_width;
    std::string line = " ";
    for (std::size_t i = 0; i + 1 < cells.size(); ++i)
    {
        line += cells[i];
        const std::size_t padding = cells[i].size() < width ? width - cells[i].size() : 1;
        line.append(padding, ' ');
    }
    line += cells.back();
    output << line << '\n';
}
void write_origin(std::ostream& output, const elecstate::PccOutput& result)
{
    if (result.slab)
    {
        const std::vector<std::string> cells = {"PCC2D_ORIGIN", "AXIS", int_text(result.axis),
                                                "COORDINATE/Bohr", real_text(result.coordinate, value_digits)};
        write_row(output, cells);
        return;
    }
    const std::vector<std::string> cells = {"PCC0D_ORIGIN",
                                            "X/Bohr", real_text(result.origin.x, value_digits),
                                            "Y/Bohr", real_text(result.origin.y, value_digits),
                                            "Z/Bohr", real_text(result.origin.z, value_digits)};
    write_row(output, cells);
}
// Smooth and point solute, polarization and screened charge: one header row,
// then one row per kind. A slab reports the dipole and second moment along the normal.
void write_moments(std::ostream& output, const elecstate::PccOutput& result)
{
    const std::string tag = result.slab ? "PCC2D_MOMENTS" : "PCC0D_MOMENTS";
    std::vector<std::string> header = {tag, "KIND", "Q/e"};
    if (result.slab)
    {
        header.push_back("PN/eBohr");
        header.push_back("QNN/eBohr2");
    }
    else
    {
        header.push_back("PX/eBohr");
        header.push_back("PY/eBohr");
        header.push_back("PZ/eBohr");
        header.push_back("QRR/eBohr2");
    }
    write_row(output, header);
    const char* labels[] = {"SMOOTH", "POINT", "POLARIZATION", "SCREENED"};
    for (int i = 0; i < 4; ++i)
    {
        const elecstate::ChargeMoments& moments = result.moments[i];
        std::vector<std::string> row = {tag, labels[i], real_text(moments.charge, value_digits)};
        if (result.slab)
        {
            const double dipole = moments.dipole * result.normal;
            row.push_back(real_text(dipole, value_digits));
        }
        else
        {
            row.push_back(real_text(moments.dipole.x, value_digits));
            row.push_back(real_text(moments.dipole.y, value_digits));
            row.push_back(real_text(moments.dipole.z, value_digits));
        }
        row.push_back(real_text(moments.second_moment, value_digits));
        write_row(output, row);
    }
}
}
namespace elecstate
{
void write_sccs_wait(std::ostream& output, double residual, double start_drho, int start_nmax)
{
    const std::vector<std::string> cells = {"SCCS_DEFERRED", "DRHO", real_text(residual, value_digits),
                                            "START_DRHO", real_text(start_drho, value_digits),
                                            "START_NMAX", int_text(start_nmax)};
    write_row(output, cells);
}

void write_sccs_out(std::ostream& output, const SccsOutput& result, int level,
                    double solvation_energy, bool pcc, double pcc_energy)
{
    if (level == 0 || !result.valid)
    {
        return;
    }
    std::vector<std::string> iteration = {"SCCS_ITER", "CG_STEPS", int_text(result.iterations),
                                          "E_SOL/Ry", real_text(solvation_energy, energy_digits)};
    if (pcc)
    {
        iteration.push_back("E_PCC/Ry");
        iteration.push_back(real_text(pcc_energy, energy_digits));
    }
    write_row(output, iteration);
    if (level < 2)
    {
        return;
    }
    const std::vector<std::string> residual = {"SCCS_RESIDUAL", "RMS", real_text(result.residual_rms, value_digits),
                                               "MAX", real_text(result.residual_max, value_digits),
                                               "WARM_START", int_text(result.warm_started)};
    write_row(output, residual);
    if (pcc)
    {
        const std::vector<std::string> gauss = {"SCCS_GAUSS",
                                                "Q_FAR_FIELD/e", real_text(result.far_field_charge, value_digits),
                                                "Q_POL/e", real_text(result.pcc.moments[2].charge, value_digits),
                                                "Q_EXPECTED/e", real_text(result.expected_charge, value_digits)};
        write_row(output, gauss);
    }
    std::vector<std::string> cavity = {"SCCS_CAVITY", "VOLUME/Bohr3", real_text(result.volume, value_digits),
                                       "SURFACE/Bohr2", real_text(result.surface, value_digits)};
    if (result.filled)
    {
        cavity.push_back("V_FILLED/Bohr3");
        cavity.push_back(real_text(result.filled_volume, value_digits));
    }
    write_row(output, cavity);
    const std::vector<std::string> transforms = {"SCCS_FFT",
                                                 "R2G_CALLS", int_text(result.transforms.forward_calls),
                                                 "G2R_CALLS", int_text(result.transforms.inverse_calls),
                                                 "CACHED_SOURCES", int_text(result.reused_fixed_sources)};
    write_row(output, transforms);
    if (!pcc)
    {
        return;
    }
    const PccOutput& data = result.pcc;
    write_origin(output, data);
    write_moments(output, data);
    // PCC energies are in Ry; ModuleBase::e2 converts Hartree to Rydberg.
    const double point_energy = pcc_energy / ModuleBase::e2;
    const std::string tag = data.slab ? "PCC2D_ENERGY" : "PCC0D_ENERGY";
    const std::vector<std::string> energies = {tag, "REACTION/Ha", real_text(result.reaction_energy, energy_digits),
                                               "PCC_SMOOTH/Ha", real_text(data.smooth_energy, energy_digits),
                                               "PCC_POINT/Ha", real_text(point_energy, energy_digits)};
    write_row(output, energies);
}
} // namespace elecstate

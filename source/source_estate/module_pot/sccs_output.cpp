#include "sccs_output.h"
#include "source_base/constants.h"
#include <algorithm>
#include <iomanip>
#include <ostream>

namespace
{
class StreamState
{
public:
    explicit StreamState(std::ostream& stream) : stream_(stream)
    {
        flags_ = stream.flags();
        precision_ = stream.precision();
    }
    ~StreamState()
    {
        stream_.flags(flags_);
        stream_.precision(precision_);
    }
private:
    std::ostream& stream_;
    std::ios_base::fmtflags flags_;
    std::streamsize precision_;
};
void write_origin(std::ostream& output, const elecstate::PccOutput& result)
{
    if (result.slab)
    {
        output << " PCC2D_ORIGIN AXIS " << result.axis
               << " COORDINATE/Bohr " << result.coordinate << '\n';
    }
    else
    {
        output << " PCC0D_ORIGIN X/Bohr " << result.origin.x
               << " Y/Bohr " << result.origin.y << " Z/Bohr " << result.origin.z << '\n';
    }
}
void write_moments(std::ostream& output, const char* label, const elecstate::ChargeMoments& moments)
{
    output << " PCC0D_MOMENTS " << label << " Q/e " << moments.charge
           << " PX/eBohr " << moments.dipole.x << " PY/eBohr " << moments.dipole.y
           << " PZ/eBohr " << moments.dipole.z << " QRR/eBohr2 " << moments.second_moment << '\n';
}
// Slab moments along the normal.
void write_slab_moments(std::ostream& output,
                        const char* label,
                        const elecstate::ChargeMoments& moments,
                        const ModuleBase::Vector3<double>& normal)
{
    const double dipole = moments.dipole * normal;
    output << " PCC2D_MOMENTS " << label << " Q/e " << moments.charge << " PN/eBohr " << dipole
           << " QNN/eBohr2 " << moments.second_moment << '\n';
}
}
namespace elecstate
{
void write_sccs_out(std::ostream& output, const SccsOutput& result, int level,
                    double solvation_energy, bool pcc, double pcc_energy)
{
    if (level == 0 || !result.valid)
    {
        return;
    }
    StreamState state(output);
    output << " SCCS_ITER " << result.iterations << " E_SOL/Ry " << std::setprecision(8) << solvation_energy;
    if (pcc)
    {
        output << " E_PCC/Ry " << pcc_energy;
    }
    output << '\n';
    if (level < 2)
    {
        return;
    }
    output << " SCCS_RESIDUAL RMS " << result.residual_rms << " MAX " << result.residual_max
           << " WARM_START " << result.warm_started << '\n';
    if (pcc)
    {
        output << " SCCS_GAUSS Q_POL_FAR_FIELD/e " << result.far_field_charge
               << " Q_POL_DENSITY/e " << result.pcc.moments[2].charge
               << " Q_POL_EXPECTED/e " << result.expected_charge << '\n';
    }
    output << " SCCS_CAVITY VOLUME/Bohr3 " << result.volume << " SURFACE/Bohr2 " << result.surface << '\n'
           << " SCCS_FFT R2G_CALLS " << result.transforms.forward_calls
           << " G2R_CALLS " << result.transforms.inverse_calls
           << " CACHED_SOURCES " << result.reused_fixed_sources << '\n';
    if (!pcc)
    {
        return;
    }
    output << std::setprecision(12);
    const PccOutput& data = result.pcc;
    write_origin(output, data);
    const char* labels[] = {"SMOOTH", "POINT", "POLARIZATION", "SCREENED"};
    for (int i = 0; i < 4; ++i)
    {
        if (data.slab)
        {
            write_slab_moments(output, labels[i], data.moments[i], data.normal);
        }
        else
        {
            write_moments(output, labels[i], data.moments[i]);
        }
    }
    const char* label = data.slab ? " PCC2D_ENERGY" : " PCC0D_ENERGY";
    const double point_energy = pcc_energy / ModuleBase::e2;
    output << label << " REACTION/Ha " << result.reaction_energy << " PCC_SMOOTH/Ha " << data.smooth_energy
           << " PCC_POINT/Ha " << point_energy << " PCC_USED/Ry " << pcc_energy << '\n';
}
} // namespace elecstate

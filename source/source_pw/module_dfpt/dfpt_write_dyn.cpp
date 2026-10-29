// ============================================================
// QE-compatible dynamical-matrix file writer (dfpt_out_dyn).
//
// One file per q point (OUT.<suffix>/dyn.q<i>), format aligned with
// QE ph.x's fildyn output so files can be diffed element-by-element
// against QE references or post-processed by external tools. The
// matrix content is the raw (pre-ASR, pre-LO-TO) UNREDUCED force
// constants Phi = D * sqrt(m_a * m_b) in Ry/bohr^2, matching the
// QE dyn convention. At Gamma with compute_q0, the Born effective
// charges and dielectric tensor are appended.
// ============================================================

#include "dfpt_pw_data.h"

#include "source_cell/unitcell.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace ModuleDFPT
{

void DFPT_Write_Dyn(const UnitCell& ucell,
                    const DFPT_PW_Data& data,
                    int nq,
                    const std::string& out_dir)
{
    ModuleBase::TITLE("DFPT_Write_Dyn", "write");
    ModuleBase::timer::start("DFPT_Write_Dyn", "write");
    const int nat = ucell.nat;
    const int nat3 = 3 * nat;
    for (int q = 0; q < nq; ++q)
    {
        if (!data.has_dynmat_raw(q) && data.get_dynmat(q).nr != nat3)
        {
            continue;
        }
        std::ostringstream name;
        name << out_dir << "/dyn.q" << q;
        std::ofstream out(name.str());
        if (!out.is_open())
        {
            continue;
        }
        out << std::setprecision(10) << std::scientific;
        const ModuleBase::ComplexMatrix dm = data.has_dynmat_raw(q) ? data.get_dynmat_raw(q)
                                                                    : data.get_dynmat(q);
        const ModuleBase::Vector3<double> q_frac = data.get_qvec(q);

        // QE dyn header
        out << "Dynamical matrix file\n";
        out << "    " << nq << "    " << nat << "\n";
        // unit-cell lattice vectors (QE writes them in alat units; here
        // we use lat0-scaled bohr, which for cubic cells matches alat)
        const double inv_scale = 1.0 / ucell.lat0;
        const ModuleBase::Matrix3 lv = ucell.latvec * inv_scale;
        out << " " << lv.e11 << " " << lv.e12 << " " << lv.e13 << "\n";
        out << " " << lv.e21 << " " << lv.e22 << " " << lv.e23 << "\n";
        out << " " << lv.e31 << " " << lv.e32 << " " << lv.e33 << "\n";
        // species and masses
        for (int it = 0; it < ucell.ntype; ++it)
        {
            out << " " << ucell.atoms[it].label << " " << ucell.atoms[it].mass;
        }
        out << "\n";
        // atom types and positions (fractional)
        for (int ia = 0; ia < nat; ++ia)
        {
            const int it = ucell.iat2it[ia];
            const int ial = ucell.iat2ia[ia];
            out << " " << it + 1 << " (" << ucell.atoms[it].tau[ial].x << ", "
                << ucell.atoms[it].tau[ial].y << ", " << ucell.atoms[it].tau[ial].z << ")\n";
        }
        // q vector
        out << "q = ( " << q_frac.x << " " << q_frac.y << " " << q_frac.z << " )\n";
        // dynamical matrix blocks (unreduced: Phi = D * sqrt(m_a m_b))
        for (int ia = 0; ia < nat; ++ia)
        {
            const double ma = ucell.atoms[ucell.iat2it[ia]].mass;
            for (int ja = 0; ja < nat; ++ja)
            {
                const double mb = ucell.atoms[ucell.iat2it[ja]].mass;
                const double mm = std::sqrt(ma * mb);
                out << "    " << ia + 1 << "    " << ja + 1 << "\n";
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        const std::complex<double> v = dm(3 * ia + i, 3 * ja + j) * mm;
                        out << " " << v.real() << " " << v.imag();
                    }
                    out << "\n";
                }
            }
        }
        // Gamma-point Born charges and dielectric tensor
        if (q == 0 && data.get_loto())
        {
            out << "\n     Dielectric constant in cartesian axis\n\n";
            const ModuleBase::matrix& eps = data.get_dielectric();
            for (int i = 0; i < 3; ++i)
            {
                out << "          (";
                for (int j = 0; j < 3; ++j)
                {
                    out << " " << std::setw(15) << eps(i, j);
                }
                out << " )\n";
            }
            out << "\n     Effective charges (d Force / dE) in cartesian axis\n\n";
            for (int ia = 0; ia < nat; ++ia)
            {
                out << "           atom      " << ia + 1 << "   " << ucell.atoms[ucell.iat2it[ia]].label << "\n";
                const ModuleBase::matrix& zstar = data.get_born(ia);
                for (int b = 0; b < 3; ++b)
                {
                    out << "      E" << static_cast<char>('x' + b) << "  (";
                    for (int a = 0; a < 3; ++a)
                    {
                        out << " " << std::setw(15) << zstar(b, a);
                    }
                    out << " )\n";
                }
            }
        }
    }
    ModuleBase::timer::end("DFPT_Write_Dyn", "write");
}

} // namespace ModuleDFPT

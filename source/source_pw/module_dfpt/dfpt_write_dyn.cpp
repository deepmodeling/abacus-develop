// ============================================================
// QE-compatible dynamical-matrix file writer (dfpt_out_dyn).
//
// One file per q point (OUT.<suffix>/dyn.q<i>), format matching
// Quantum ESPRESSO ph.x's fildyn output (io_dyn_mat_old.f90):
// the file can be parsed by QE's read routine or diffed against
// QE references. The matrix content is the raw (pre-ASR, pre-LO-TO)
// unreduced force constants Phi = D * sqrt(m_a * m_b) in Ry/bohr^2.
// At Gamma with dfpt_compute_q0, the Born effective charges and
// dielectric tensor are appended in the QE convention.
//
// Only the ionode (rank 0 of Comm_Comm) writes, to avoid competing
// truncations of the same relative path on a shared directory.
// ============================================================

#include "dfpt_write_dyn.h"

#include "dfpt_pw_data.h"

#include "source_cell/unitcell.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace ModuleDFPT
{

/// QE's mass unit conversion: dyn files store mass[amu] * AMU_RY,
/// where AMU_RY = AMU_AU / 2 = 911.444... (constants.f90 line 57)
static double amu_to_ry(double mass_amu)
{
    return mass_amu * 911.4442431045;
}

void DFPT_Write_Dyn(const UnitCell& ucell,
                    const DFPT_PW_Data& data,
                    int nq,
                    const std::string& out_dir,
                    int my_rank,
                    std::ostream* running_log)
{
    ModuleBase::TITLE("ModuleDFPT", "DFPT_Write_Dyn");
    ModuleBase::timer::start("ModuleDFPT", "DFPT_Write_Dyn");

    // only the output rank writes, to avoid competing truncations;
    // the rank is passed by the caller (no global-state read here)
    if (my_rank != 0)
    {
        ModuleBase::timer::end("ModuleDFPT", "DFPT_Write_Dyn");
        return;
    }

    const int nat = ucell.nat;
    const int nat3 = 3 * nat;
    const double alat = ucell.lat0; // bohr, the QE celldm(1)

    for (int q = 0; q < nq; ++q)
    {
        const ModuleBase::ComplexMatrix dm = data.has_dynmat_raw(q) ? data.get_dynmat_raw(q)
                                                                    : data.get_dynmat(q);
        if (dm.nr != nat3)
        {
            continue;
        }

        std::ostringstream name;
        name << out_dir << "/dyn.q" << q;
        std::ofstream out(name.str());
        if (!out.is_open())
        {
            // report the full failed path and skip this q
            std::ostringstream err;
            err << "DFPT_Write_Dyn: cannot open " << name.str() << " for writing; skipping q="
                << q;
            ModuleBase::WARNING("DFPT_Write_Dyn", err.str());
            continue;
        }
        out << std::setprecision(10) << std::scientific;

        // --- QE dyn header (io_dyn_mat_old.f90 write_dyn_mat_header) ---
        out << "Dynamical matrix file\n";
        out << "ABACUS DFPT output                                                          \n";
        // ntyp nat ibrav celldm(1..6); ibrav=0 (basis vectors follow)
        out << "  " << ucell.ntype << "  " << nat << "  0  " << alat
            << "  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000" << "\n";
        out << "Basis vectors\n";
        // latvec is already in lat0 (= alat) units (unitcell_setup.cpp:
        // "latvec have the unit of lat0"); write directly, no division
        for (int i = 0; i < 3; ++i)
        {
            const ModuleBase::Matrix3& lv = ucell.latvec;
            out << "  ";
            if (i == 0)
            {
                out << lv.e11 << "  " << lv.e12 << "  " << lv.e13;
            }
            else if (i == 1)
            {
                out << lv.e21 << "  " << lv.e22 << "  " << lv.e23;
            }
            else
            {
                out << lv.e31 << "  " << lv.e32 << "  " << lv.e33;
            }
            out << "\n";
        }
        // species: index, label, mass * AMU_RY (io_dyn_mat_old.f90 line 78)
        for (int it = 0; it < ucell.ntype; ++it)
        {
            const double mass_file = amu_to_ry(ucell.atoms[it].mass);
            out << "  " << it + 1 << "  '" << ucell.atoms[it].label << "'  " << mass_file << "\n";
        }
        // atoms: global index, type index, tau in alat units
        for (int ia = 0; ia < nat; ++ia)
        {
            const int it = ucell.iat2it[ia];
            const int ial = ucell.iat2ia[ia];
            const ModuleBase::Vector3<double> tau_alat = ucell.atoms[it].tau[ial];
            out << "  " << ia + 1 << "  " << it + 1 << "  " << tau_alat.x << "  " << tau_alat.y
                << "  " << tau_alat.z << "\n";
        }

        // --- q vector in cartesian 2*pi/alat units (QE convention) ---
        const ModuleBase::Vector3<double> q_frac = data.get_qvec(q);
        // q_cart = q_frac * G, already in 2*pi/lat0 = 2*pi/alat units
        const ModuleBase::Vector3<double> q_cart = q_frac * ucell.G;
        out << "\n     Dynamical  Matrix in cartesian axes\n\n";
        out << "     q = ( " << q_cart.x << "  " << q_cart.y << "  " << q_cart.z << " )\n\n";

        // --- force-constant matrix blocks (unreduced) ---
        for (int ia = 0; ia < nat; ++ia)
        {
            const double ma = ucell.atoms[ucell.iat2it[ia]].mass;
            for (int ja = 0; ja < nat; ++ja)
            {
                const double mb = ucell.atoms[ucell.iat2it[ja]].mass;
                const double mass_product = ma * mb;
                const double mm = std::sqrt(mass_product);
                out << "  " << ia + 1 << "  " << ja + 1 << "\n";
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        const int row = 3 * ia + i;
                        const int col = 3 * ja + j;
                        const std::complex<double> phi = dm(row, col) * mm;
                        out << " " << phi.real() << "  " << phi.imag();
                    }
                    out << "\n";
                }
            }
        }

        // --- dielectric tensor and Born charges at Gamma (QE order:
        // Dielectric Tensor first, then Effective Charges) ---
        // gate on the computed q0 response (compute_q0) at the actual
        // Gamma entry -- a q-vector check, not a list-index check, so a
        // reordered explicit q file cannot attach the tensors to a
        // nonzero-q entry (setup_q_list additionally rejects non-Gamma-
        // first orderings when compute_q0 is on, keeping the q0 solver
        // and this export consistent); independent of get_loto(): with
        // compute_q0=true and loto=false the tensors are still computed
        // and should be exported
        const double gamma_tol = 1.0e-8; ///< empirical parameter: fractional-q Gamma tolerance (dfpt_phon.cpp value)
        const bool q_is_gamma = std::abs(q_frac.x) <= gamma_tol && std::abs(q_frac.y) <= gamma_tol
                                && std::abs(q_frac.z) <= gamma_tol;
        if (q_is_gamma && data.get_compute_q0())
        {
            const ModuleBase::matrix& eps = data.get_dielectric();
            const ModuleBase::matrix& born0 = data.get_born(0);
            const bool has_eps = eps.nr == 3 && eps.nc == 3;
            const bool has_born = born0.nr == 3 && born0.nc == 3;
            if (has_eps)
            {
                out << "\n     Dielectric Tensor:\n\n";
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        out << "  " << std::setw(22) << eps(i, j);
                    }
                    out << "\n";
                }
            }
            if (has_born)
            {
                out << "\n     Effective Charges E-U: Z_{alpha}{s,beta}\n\n";
                for (int ia = 0; ia < nat; ++ia)
                {
                    const ModuleBase::matrix& zstar = data.get_born(ia);
                    out << "     atom #  " << ia + 1 << "\n";
                    for (int b = 0; b < 3; ++b)
                    {
                        for (int a = 0; a < 3; ++a)
                        {
                            out << "  " << std::setw(22) << zstar(b, a);
                        }
                        out << "\n";
                    }
                }
            }
        }

        out.close();
        if (out.fail())
        {
            std::ostringstream err;
            err << "DFPT_Write_Dyn: write error on " << name.str();
            ModuleBase::WARNING("DFPT_Write_Dyn", err.str());
        }
        else
        {
            // log the created file path (AGENTS rule 16): the running log
            // receives the record through the explicit running_log
            // dependency (serial runs open a separate running log stream
            // in global_file.cpp, so stdout alone does not reach it);
            // stdout keeps a terminal echo
            std::ostringstream ok;
            ok << "DFPT_Write_Dyn: wrote " << name.str() << std::endl;
            if (running_log != nullptr)
            {
                (*running_log) << " " << ok.str();
            }
            std::cout << " " << ok.str();
        }
    }
    ModuleBase::timer::end("ModuleDFPT", "DFPT_Write_Dyn");
}

} // namespace ModuleDFPT

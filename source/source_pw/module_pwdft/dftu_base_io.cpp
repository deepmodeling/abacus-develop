#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/timer.h"
#include "source_main/version.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// local helpers for eigenvalue calculation
// migrated from dftu_base.cpp, mohan 2025-11-08
namespace
{

inline void JacobiRotate(std::vector<std::vector<double>>& A, int p, int q, int n)
{
    if (std::abs(A[p][q]) > 1e-10)
    {
        double r = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t = 0.0;
        if (r >= 0)
        {
            t = 1.0 / (r + sqrt(1.0 + r * r));
        }
        else
        {
            t = -1.0 / (-r + sqrt(1.0 + r * r));
        }
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        A[p][p] -= t * A[p][q];
        A[q][q] += t * A[p][q];
        A[p][q] = A[q][p] = 0.0;

        for (int k = 0; k < n; k++)
        {
            if (k != p && k != q)
            {
                double Akp = c * A[k][p] - s * A[k][q];
                double Akq = s * A[k][p] + c * A[k][q];
                A[k][p] = A[p][k] = Akp;
                A[k][q] = A[q][k] = Akq;
            }
        }
    }
}

inline std::vector<double> CalculateEigenvalues(std::vector<std::vector<double>>& A, int n)
{
    std::vector<double> eigenvalues(n);
    while (true)
    {
        int p = 0, q = 1;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (std::abs(A[i][j]) > std::abs(A[p][q]))
                {
                    p = i;
                    q = j;
                }
            }
        }

        if (std::abs(A[p][q]) < 1e-10)
        {
            for (int i = 0; i < n; i++)
            {
                eigenvalues[i] = A[i][i];
            }
            break;
        }

        JacobiRotate(A, p, q, n);
    }
    return eigenvalues;
}

/// @brief Extract the four Pauli-component blocks of one SOC correlated shell.
///
/// The blocks, each (2l+1)x(2l+1) stored row-major with index m0*m+m1, are
/// the charge channel and the sigma_x/y/z spin channels. The PW storage keeps
/// the four blocks contiguous; the LCAO storage is a real symmetric matrix in
/// spin basis, whose spin-off-diagonal imaginary channel is not stored.
///
/// @param dftu DFT+U object holding the occupation matrix
/// @param iat atom index
/// @param l angular-momentum channel
/// @param layout storage layout of the nspin == 4 occupation matrix
/// @param blocks [out] four Pauli-component blocks
inline void extract_soc_pauli_blocks(const Plus_U_Base& dftu,
                                     const int iat,
                                     const int l,
                                     const DFTU_BASE::OccmatSocLayout layout,
                                     std::vector<std::vector<double>>& blocks)
{
    const int m = 2 * l + 1;
    const int m2 = m * m;
    blocks.assign(4, std::vector<double>(m2, 0.0));

    const OccupationMatrix& occmat = dftu.occmat();
    if (layout == DFTU_BASE::SOC_LAYOUT_PAULI)
    {
        const ModuleBase::matrix& occ = occmat.mat(iat, l, 0);
        for (int is = 0; is < 4; ++is)
        {
            for (int k = 0; k < m2; ++k)
            {
                blocks[is][k] = occ.c[is * m2 + k];
            }
        }
    }
    else // SOC_LAYOUT_SPIN_BASIS_REAL
    {
        for (int m0 = 0; m0 < m; ++m0)
        {
            for (int m1 = 0; m1 < m; ++m1)
            {
                const int k = m0 * m + m1;
                const double suu = occmat.get(iat, l, 0, m0, m1);
                const double sud = occmat.get(iat, l, 0, m0, m1 + m);
                const double sdu = occmat.get(iat, l, 0, m0 + m, m1);
                const double sdd = occmat.get(iat, l, 0, m0 + m, m1 + m);
                blocks[0][k] = suu + sdd;
                blocks[1][k] = sud + sdu;
                // The spin-basis storage drops imaginary parts, so the
                // sigma_y channel has no contribution.
                blocks[2][k] = 0.0;
                blocks[3][k] = suu - sdd;
            }
        }
    }
}

/// @brief Diagonalize the SOC occupation matrix reconstructed in spinor space.
///
/// The Hermitian spinor matrix is built from the Pauli blocks as
///   n_upup   = (b0 + b3)/2, n_dndn = (b0 - b3)/2,
///   n_updown = (b1 + i*b2)/2.
/// Its 2*(2l+1) real eigenvalues are computed with LAPACK zheev.
///
/// @param blocks four Pauli-component blocks
/// @param m 2l+1
/// @return eigenvalues sorted ascending by zheev
inline std::vector<double> calculate_spinor_eigenvalues(
    const std::vector<std::vector<double>>& blocks,
    const int m)
{
    const int n = 2 * m;
    ModuleBase::ComplexMatrix spinor_n(n, n);
    for (int m0 = 0; m0 < m; ++m0)
    {
        for (int m1 = 0; m1 < m; ++m1)
        {
            const int k = m0 * m + m1;
            const double b0 = blocks[0][k];
            const double b1 = blocks[1][k];
            const double b2 = blocks[2][k];
            const double b3 = blocks[3][k];

            spinor_n(m0, m1) = std::complex<double>(0.5 * (b0 + b3), 0.0);
            spinor_n(m0, m1 + m) = std::complex<double>(0.5 * b1, 0.5 * b2);
            spinor_n(m0 + m, m1) = std::complex<double>(0.5 * b1, -0.5 * b2);
            spinor_n(m0 + m, m1 + m) = std::complex<double>(0.5 * (b0 - b3), 0.0);
        }
    }

    std::vector<double> eigenvalues(n, 0.0);
    std::vector<std::complex<double>> work(1, std::complex<double>(0.0, 0.0));
    std::vector<double> rwork(std::max(1, 3 * n - 2), 0.0);
    int info = 0;
    int lwork = -1;

    LapackConnector::zheev('N', 'U', n, spinor_n, n, eigenvalues.data(),
                           work.data(), lwork, rwork.data(), &info);
    if (info != 0)
    {
        ModuleBase::WARNING_QUIT("calculate_spinor_eigenvalues",
                                 "LAPACK zheev workspace query failed");
    }
    lwork = std::max(1, static_cast<int>(work[0].real()));
    work.resize(lwork);
    LapackConnector::zheev('N', 'U', n, spinor_n, n, eigenvalues.data(),
                           work.data(), lwork, rwork.data(), &info);
    if (info != 0)
    {
        ModuleBase::WARNING_QUIT("calculate_spinor_eigenvalues",
                                 "LAPACK zheev diagonalization failed");
    }
    return eigenvalues;
}

} // namespace


namespace
{

// Write the STRU-style provenance header shared by occ_mat.txt and the
// per-ionic-step snapshot files: ABACUS version, local timestamp and the
// 1-based relaxation step this file belongs to.
void write_provenance_header(std::ostream& os, const int istep)
{
    std::time_t now = std::time(nullptr);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    os << "# ABACUS version: " << VERSION << "\n";
    os << "# Written at " << time_buf << "\n";
    os << "# RELAX STEP " << istep + 1 << "\n";
}

// Consume the label of one SOC block in the compact layout, e.g.
// "spin 1 nelec 1.07922327" or "spin 12 re". The optional
// "nelec <value>" token pair is skipped when present.
void skip_compact_soc_label(std::ifstream& ifs, char* word)
{
    ifs >> word; // "spin"
    ifs >> word; // block index: "1", "2" or "12"
    ifs >> word; // "nelec", "re" or "im"
    if (strcmp(word, "nelec") == 0)
    {
        double nelec = 0.0;
        ifs >> nelec;
    }
}

// Read n*n doubles of one compact matrix block row by row. When out is
// non-null the values are stored, otherwise they are read and discarded.
void read_matrix_block(std::ifstream& ifs, const int n, double* out)
{
    for (int i = 0; i < n * n; ++i)
    {
        double value = 0.0;
        ifs >> value;
        if (out != nullptr)
        {
            out[i] = value;
        }
    }
}

} // namespace


namespace DFTU_BASE
{

bool is_ion_step_output_step(const int istep, const OccmatOutputCfg& cfg)
{
    if (cfg.out_freq_ion <= 0)
    {
        return false;
    }
    return (istep % cfg.out_freq_ion == 0);
}

bool is_elec_snapshot_trigger(const int iter,
                              const bool conv_esolver,
                              const OccmatOutputCfg& cfg)
{
    const bool periodic = (iter % cfg.out_freq_elec == 0);
    const bool last_step = (iter == cfg.scf_nmax);
    return periodic || last_step || conv_esolver;
}

std::string gen_ion_step_occ_mat_filename(const std::string& out_dir, const int istep)
{
    std::stringstream ss;
    ss << out_dir << "occ_matg" << (istep + 1) << ".txt";
    return ss.str();
}

std::string find_first_existing_file(const std::string& dir,
                                     const std::vector<std::string>& candidates)
{
    for (const std::string& name : candidates)
    {
        const std::string full = dir + name;
        std::ifstream probe(full.c_str());
        if (probe.is_open())
        {
            return full;
        }
    }
    return std::string();
}

void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const int init_occ_mat,
                  const std::string& fn,
                  const std::string& init_chg,
                  int nspin,
                  int npol,
                  OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "read_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ifstream ifdftu(fn.c_str(), std::ios::in);

    if (!ifdftu)
    {
        if (init_occ_mat > 0)
        {
            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "Can not find the file dm_onsite_ini.txt. Please check your dm_onsite_ini.txt");
        }
        else
        {
            if (init_chg == "file")
            {
                const std::string not_found_msg = "Can not find the file " + fn
                                                  + ". Please do scf calculation first";
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", not_found_msg);
            }
        }
        const std::string open_failed_msg = "Can not open " + fn + " file";
        ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", open_failed_msg);
    }

    GlobalV::ofs_running << " DFT+U: read occupation matrix from " << fn << std::endl;
    std::cout << " DFT+U: read occupation matrix from " << fn << std::endl;

    ifdftu.clear();
    ifdftu.seekg(0);

    char word[20];

    int T = 0;
    int iat = 0;
    int spin = 0;
    int L = 0;
    int zeta = 0;

    while (true)
    {
        if (!(ifdftu >> word))
        {
            break;
        }

        // Comment lines ("# ...") carry the provenance header and the
        // electronic-step marker in the compact layout; skip the line.
        if (word[0] == '#')
        {
            ifdftu.ignore(150, '\n');
            continue;
        }

        // The legacy "Atom=" token and the compact "Atom" label both
        // start one atom block; every other token (element prefix, log
        // sentences) is ignored.
        const bool legacy_atom = (strcmp(word, "Atom=") == 0);
        const bool compact_atom = (strcmp(word, "Atom") == 0);
        if (!legacy_atom && !compact_atom)
        {
            continue;
        }
        const bool legacy_layout = legacy_atom;

        ifdftu >> iat;
        iat -= 1;

        if (legacy_layout)
        {
            ifdftu >> word;

            // Accept both the compact "L=2" and the legacy split "L= 2".
            if (strncmp("L=", word, 2) != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            if (strlen(word) > 2)
            {
                L = atoi(word + 2);
            }
            else
            {
                ifdftu >> L;
            }
            ifdftu >> word;

            // The ORBITAL= token is optional in newly written files; when
            // absent the radial channel defaults to 0.
            if (strncmp("ORBITAL=", word, 8) == 0)
            {
                if (strlen(word) > 8)
                {
                    zeta = atoi(word + 8);
                }
                else
                {
                    ifdftu >> zeta;
                }
            }
            else
            {
                zeta = 0;
            }
            ifdftu.ignore(150, '\n');

            if (zeta != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m",
                                         "only the first radial channel (ORBITAL=0) is supported");
            }
        }
        else
        {
            // Compact header, e.g. " Fe Atom 1 L 2 mag 0.03..."
            ifdftu >> word;
            if (strcmp(word, "L") != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> L;
            // The optional "mag ..." tail is redundant for a restart.
            ifdftu.ignore(150, '\n');
            zeta = 0;
        }

        T = ucell.iat2it[iat];
        const int NL = ucell.atoms[T].nwl + 1;

        for (int l = 0; l < NL; l++)
        {
            if (l != l_channel[T])
            {
                continue;
            }

            const int nm = 2 * L + 1;

            if (nspin == 1 || nspin == 2)
            {
                for (int is = 0; is < 2; is++)
                {
                    if (legacy_layout)
                    {
                        ifdftu >> word;
                        if (strcmp("spin=", word) != 0)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                        ifdftu >> spin;
                        spin -= 1;
                        ifdftu.ignore(150, '\n');

                        double value = 0.0;
                        for (int m0 = 0; m0 < nm; m0++)
                        {
                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                ifdftu >> value;
                                occ.set(iat, L, spin, m0, m1, value);
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                    else
                    {
                        // Compact label: "spin 1 nelec 0.468..."
                        ifdftu >> word;
                        if (strcmp("spin", word) != 0)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                        int spin_idx = 0;
                        ifdftu >> spin_idx;
                        if (spin_idx != is + 1)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG SPIN INDEX IN Plus_U FILE");
                        }
                        ifdftu >> word;
                        if (strcmp("nelec", word) != 0)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                        double nelec = 0.0;
                        ifdftu >> nelec;

                        double value = 0.0;
                        for (int m0 = 0; m0 < nm; m0++)
                        {
                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                ifdftu >> value;
                                occ.set(iat, L, is, m0, m1, value);
                            }
                        }
                    }
                }
            }
            else if (nspin == 4) // SOC
            {
                if (legacy_layout)
                {
                    double value = 0.0;
                    for (int m0 = 0; m0 < nm; m0++)
                    {
                        for (int ipol0 = 0; ipol0 < npol; ipol0++)
                        {
                            const int m0_all = m0 + nm * ipol0;

                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                {
                                    int m1_all = m1 + nm * ipol1;
                                    ifdftu >> value;
                                    occ.set(iat, L, 0, m0_all, m1_all, value);
                                }
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                }
                else
                {
                    // Four labeled blocks, in the same order they are
                    // written: up/up, up/down Re, up/down Im, down/down.
                    std::vector<double> uu(nm * nm);
                    std::vector<double> re(nm * nm);
                    std::vector<double> im(nm * nm);
                    std::vector<double> dd(nm * nm);

                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, uu.data());

                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, re.data());

                    // The imaginary block is needed by the PW (Pauli)
                    // path to reconstruct b2 = 2*Im(n_ud); the LCAO
                    // (real spin-basis) path discards it because its
                    // storage can not represent Im(n_ud).
                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, im.data());

                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, dd.data());

                    if (soc_layout == SOC_LAYOUT_PAULI)
                    {
                        // Reconstruct the 4 contiguous Pauli blocks
                        // [b0, b1, b2, b3] in the 2m x 2m flat buffer
                        // from the file's (n_uu, Re(n_ud), Im(n_ud),
                        // n_dd) representation. The writer prints
                        //   n_uu     = (b0 + b3)/2
                        //   Re(n_ud) = b1/2
                        //   Im(n_ud) = b2/2
                        //   n_dd     = (b0 - b3)/2
                        // so the inverse is
                        //   b0 = n_uu + n_dd
                        //   b1 = 2 * Re(n_ud)
                        //   b2 = 2 * Im(n_ud)
                        //   b3 = n_uu - n_dd.
                        ModuleBase::matrix& occ0 = occ.mat(iat, L, 0);
                        const int m2 = nm * nm;
                        for (int k = 0; k < m2; ++k)
                        {
                            occ0.c[0 * m2 + k] = uu[k] + dd[k];
                            occ0.c[1 * m2 + k] = 2.0 * re[k];
                            occ0.c[2 * m2 + k] = 2.0 * im[k];
                            occ0.c[3 * m2 + k] = uu[k] - dd[k];
                        }
                    }
                    else // SOC_LAYOUT_SPIN_BASIS_REAL
                    {
                        // Reconstruct the real 2m x 2m spin-basis matrix.
                        // Mirror the legacy writer: the Re block fills both
                        // (up, down) and (down, up) at the same (m0, m1).
                        // Im(n_ud) is dropped because the storage can not
                        // represent it.
                        for (int m0 = 0; m0 < nm; m0++)
                        {
                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                const int k = m0 * nm + m1;
                                occ.set(iat, L, 0, m0, m1, uu[k]);
                                occ.set(iat, L, 0, m0, nm + m1, re[k]);
                                occ.set(iat, L, 0, nm + m0, m1, re[k]);
                                occ.set(iat, L, 0, nm + m0, nm + m1, dd[k]);
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

#ifdef __MPI
/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Each occupation matrix is broadcast as one contiguous block
/// (matrix::c stores nr * nc consecutive doubles) instead of element by
/// element.
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& l_channel,
                       int nspin,
                       int npol)
{
    ModuleBase::TITLE("DFTU_BASE", "local_occup_bcast");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (l_channel[T] == -1)
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            const int L = l_channel[T];

            for (int l = 0; l <= ucell.atoms[T].nwl; l++)
            {
                if (l != l_channel[T])
                {
                    continue;
                }

                if (nspin == 1 || nspin == 2)
                {
                    for (int spin = 0; spin < 2; spin++)
                    {
                        Parallel_Common::bcast_double(occ.mat(iat, l, spin).c,
                                                      occ.mat(iat, l, spin).nr * occ.mat(iat, l, spin).nc);
                    }
                }
                else if (nspin == 4) // SOC
                {
                    Parallel_Common::bcast_double(occ.mat(iat, l, 0).c,
                                                  occ.mat(iat, l, 0).nr * occ.mat(iat, l, 0).nc);
                }
            }
        }
    }
    return;
}
#endif


void prepare_ion_step_file(const std::string& global_out_dir,
                           const int istep,
                           const OccmatOutputCfg& cfg)
{
    // Only the root process creates the per-ionic-step file.
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    // trunc makes a rerun start from a clean file instead of appending
    // snapshots left over from a previous calculation.
    const std::string ion_step_fn = gen_ion_step_occ_mat_filename(global_out_dir, istep);
    std::ofstream ofs_ion_step;
    ofs_ion_step.open(ion_step_fn, std::ios::out | std::ios::trunc);
    if (!ofs_ion_step)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::prepare_ion_step_file",
                                 "Can't create per-ionic-step occupation-matrix file");
    }
    // Provenance header shared with occ_mat.txt.
    write_provenance_header(ofs_ion_step, istep);
    ofs_ion_step.close();
    return;
}


void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            const std::string& global_out_dir,
            int nspin,
            int npol,
            int istep,
            int iter,
            const OccmatOutputCfg& cfg,
            OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "output");

    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output", "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output", "iter must be >= 1");
    }
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output",
                                "invalid occupation-matrix output frequency configuration");
    }

    // The per-type U values and the per-atom occupation matrices are not
    // dumped into the running log here: for large cells (thousands of
    // atoms) this floods the log at every electronic step. The dedicated
    // occupation-matrix files (occ_mat.txt, occ_matg{#}.txt) gated by
    // out_occ_mat carry the same information instead.

    // occ_mat.txt is not written here: its section records drho, which
    // is only available after the electronic solve. write_latest_occmat()
    // overwrites it at the iter_finish stage.

    // At the first electronic step of an output ionic step, (re)create the
    // per-ionic-step file and write its header. The LCAO path enters here
    // directly; the PW path additionally prepares the file from
    // iter_init_dftu_pw(), because it skips output() at istep 0 / iter 1.
    // Electronic-step sections are appended later by append_ion_step_snapshot().
    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    if (cfg.out_occ_mat && ion_step_output && iter == 1)
    {
        prepare_ion_step_file(global_out_dir, istep, cfg);
    }

    return;
}

// Write one electronic-step section shared by occ_mat.txt and the
// per-ionic-step snapshot files: the step marker, the configured
// charge-density convergence threshold, the actual charge-density
// residual and the full occupation matrices. When occmat_ready is false,
// the matrix body is an "N/A" placeholder.
void write_snapshot_section(std::ostream& os,
                            const Plus_U_Base& dftu,
                            const UnitCell& ucell,
                            const int iter,
                            const double scf_thr,
                            const double drho,
                            const bool occmat_ready,
                            const int nspin,
                            const int npol,
                            const OccmatSocLayout soc_layout)
{
    os << "# Electronic step " << iter << "\n";
    os << "# scf_thr " << std::scientific << std::setprecision(8) << scf_thr << "\n";
    os << "# drho " << std::scientific << std::setprecision(8) << drho << "\n";

    // The PW path has no occupation matrix at istep 0 / iter 1 unless it was
    // loaded from file; record an explicit placeholder instead of a silent
    // zero matrix. read_occup_m() skips the token like any other non-atom
    // word, so the file stays parseable on restart.
    if (!occmat_ready)
    {
        os << "\n N/A\n";
        return;
    }

    // Full occupation matrices together with per-atom magnetism.
    write_occup_m(dftu, ucell, os, true, nspin, npol,
                  OCMAT_FMT_READABLE, soc_layout);
}

void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              int nspin,
                              int npol,
                              int istep,
                              int iter,
                              bool conv_esolver,
                              bool occmat_ready,
                              double scf_thr,
                              double drho,
                              const OccmatOutputCfg& cfg,
                              OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "append_ion_step_snapshot");

    // dft_plus_u <= 0 means the occupation matrix was never computed
    // (Plus_U_Base::l_channel stays empty); writing here would dereference
    // a null l_channel.data() in write_occup_m. Match the documented contract
    // that out_occ_mat only takes effect for DFT+U calculations. Check this
    // before validating iter/istep/nspin: when DFT+U is off this function is
    // a no-op, so an invalid iter (e.g. from an upstream bug) must not crash
    // here. Mirrors the guard-first pattern in write_latest_occmat.
    if (!cfg.out_occ_mat || cfg.dft_plus_u <= 0)
    {
        return;
    }

    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "nspin must be 1, 2 or 4");
    }
    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "iter must be >= 1");
    }
    // out_freq_ion == 0 is the valid default (no numbered file) and is
    // handled by the gates below; only negative values are an error.
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                "invalid occupation-matrix output frequency configuration");
    }

    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    const bool elec_trigger = is_elec_snapshot_trigger(iter, conv_esolver, cfg);
    if (!ion_step_output || !elec_trigger || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string ion_step_fn = gen_ion_step_occ_mat_filename(global_out_dir, istep);

    // The provenance header written by prepare_ion_step_file() must not be
    // counted as a recorded section: scan the file for an existing section
    // marker instead of testing the stream position.
    bool section_exists = false;
    {
        std::ifstream checker(ion_step_fn.c_str());
        std::string line;
        while (std::getline(checker, line))
        {
            if (line.rfind("# Electronic step", 0) == 0)
            {
                section_exists = true;
                break;
            }
        }
    }

    std::ofstream ofs_ion_step;
    ofs_ion_step.open(ion_step_fn, std::ios::app);
    if (!ofs_ion_step)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                 "Can't open per-ionic-step occupation-matrix file");
    }

    // Separate consecutive sections by a blank line; the first section
    // directly follows the provenance header.
    if (section_exists)
    {
        ofs_ion_step << "\n";
    }
    write_snapshot_section(ofs_ion_step, dftu, ucell, iter, scf_thr, drho,
                           occmat_ready, nspin, npol, soc_layout);

    ofs_ion_step.close();

    GlobalV::ofs_running << " Append local occupation matrices of electronic step " << iter
                         << " to file: " << ion_step_fn << std::endl;
    return;
}


void write_latest_occmat(const Plus_U_Base& dftu,
                         const UnitCell& ucell,
                         const std::string& global_out_dir,
                         int nspin,
                         int npol,
                         int istep,
                         int iter,
                         double scf_thr,
                         double drho,
                         const OccmatOutputCfg& cfg,
                         OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "write_latest_occmat");

    if (!cfg.out_occ_mat || cfg.dft_plus_u <= 0 || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string latest_fn = global_out_dir + "occ_mat.txt";
    std::ofstream ofdftu;
    ofdftu.open(latest_fn, std::ios::out | std::ios::trunc);
    if (!ofdftu)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::write_latest_occmat", "Can't create file occ_mat.txt");
    }

    // occ_mat.txt is a single-section snapshot file: the provenance
    // header shared with the g files followed by the same section writer.
    write_provenance_header(ofdftu, istep);
    // Callers of write_latest_occmat() guarantee the matrix exists: the PW
    // esolver gates on latest_ready, the LCAO esolver computes it every
    // electronic iteration.
    write_snapshot_section(ofdftu, dftu, ucell, iter, scf_thr, drho,
                           true, nspin, npol, soc_layout);

    ofdftu.close();
    return;
}


void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ostream& ofs,
                   bool diag,
                   int nspin,
                   int npol,
                   OccmatTextFormat fmt,
                   OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "write_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!dftu.has_l_channel(T))
        {
            continue;
        }
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = dftu.get_l_channel(T);

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < NL; l++)
            {
                if (l != dftu.get_l_channel(T))
                {
                    continue;
                }

                const std::string& elem_label = ucell.atoms[T].label;
                const int index_in_type = I + 1;
                const ModuleBase::Vector3<double>& tau = ucell.atoms[T].tau[I];
                const bool readable = (fmt == OCMAT_FMT_READABLE);

                // In the compact collinear/SOC layout the header carries the
                // per-atom magnetism, which is only known after the Pauli
                // block traces; that header is written in the spin branches.
                const bool defer_header = readable
                                          && (nspin == 1 || nspin == 2 || nspin == 4);
                if (!defer_header)
                {
                    if (readable)
                    {
                        // Readable header, e.g. " Fe Atom 1 L 2".
                        // The position is not repeated: it is recorded in
                        // STRU and constant within one ionic (g) step.
                        ofs << "\n " << elem_label;
                        ofs << " Atom " << iat + 1;
                        ofs << " L " << l << std::endl;
                    }
                    else
                    {
                        // Legacy token layout parsed by read_occup_m(); the
                        // "<Element> <index>" prefix is skipped when scanning
                        // for "Atom=".
                        ofs << "\n" << elem_label << " " << index_in_type;
                        ofs << " Atom=" << iat + 1;
                        ofs << " L=" << l;
                        ofs << std::setprecision(8) << std::fixed
                            << std::setw(14) << tau.x
                            << std::setw(14) << tau.y
                            << std::setw(14) << tau.z << std::endl;
                    }
                }

                if (nspin == 1 || nspin == 2)
                {
                    const int nm = 2 * l + 1;
                    double sum0[2] = {0.0, 0.0};

                    if (readable)
                    {
                        // Compact snapshot layout, e.g.
                        // " Fe Atom 1 L 2 mag 0.84215429"
                        // " spin 1 nelec 1.76153075"
                        // followed by the matrix rows.
                        if (diag)
                        {
                            for (int is = 0; is < 2; is++)
                            {
                                for (int m0 = 0; m0 < nm; m0++)
                                {
                                    sum0[is] += dftu.occmat().get(iat, l, is, m0, m0);
                                }
                            }
                        }
                        ofs << "\n " << elem_label;
                        ofs << " Atom " << iat + 1;
                        ofs << " L " << l;
                        if (diag)
                        {
                            ofs << " mag " << std::fixed << std::setprecision(8)
                                << sum0[0] - sum0[1];
                        }
                        ofs << std::endl;

                        ofs << std::fixed << std::setprecision(8);
                        for (int is = 0; is < 2; is++)
                        {
                            ofs << " spin " << is + 1;
                            if (diag)
                            {
                                ofs << " nelec " << sum0[is];
                            }
                            ofs << std::endl;
                            for (int m0 = 0; m0 < nm; m0++)
                            {
                                for (int m1 = 0; m1 < nm; m1++)
                                {
                                    // One separating space before every
                                    // value; positive values are 10 chars,
                                    // negative ones 11.
                                    ofs << " " << std::setw(10)
                                        << dftu.occmat().get(iat, l, is, m0, m1);
                                }
                                ofs << std::endl;
                            }
                        }
                    }
                    else
                    {
                        const std::string trace_label = " Trace (electrons) ";
                        const std::string matrix_label = " spin= ";
                        for (int is = 0; is < 2; is++)
                        {
                            if (diag)
                            {
                                std::vector<std::vector<double>> A(nm,
                                                                   std::vector<double>(nm));
                                for (int m0 = 0; m0 < nm; m0++)
                                {
                                    for (int m1 = 0; m1 < nm; m1++)
                                    {
                                        A[m0][m1] = dftu.occmat().get(iat, l, is, m0, m1);
                                    }
                                }
                                std::vector<double> eigenvalues = CalculateEigenvalues(A, nm);
                                sum0[is] = 0.0;
                                ofs << " Eigenvalues for spin=" << is + 1 << std::endl;
                                ofs << std::setprecision(8) << std::fixed;
                                for (int i = 0; i < nm; i++)
                                {
                                    ofs << std::setw(12) << eigenvalues[i];
                                    sum0[is] += eigenvalues[i];
                                }
                                ofs << std::endl;
                                ofs << std::fixed << std::setprecision(8);
                                ofs << trace_label << sum0[is] << std::endl;
                            }
                            ofs << matrix_label << is + 1 << std::endl;
                            // Matrix elements use fixed-point notation with 8
                            // fractional digits (~1e-8 absolute precision).
                            ofs << std::fixed << std::setprecision(8);
                            for (int m0 = 0; m0 < nm; m0++)
                            {
                                for (int m1 = 0; m1 < nm; m1++)
                                {
                                    // Legacy keeps its historical
                                    // 12-character column alignment.
                                    ofs << " " << std::setw(11)
                                        << dftu.occmat().get(iat, l, is, m0, m1);
                                }
                                ofs << std::endl;
                            }
                        }
                        if (diag)
                        {
                            ofs << " mag " << std::fixed << std::setprecision(8)
                                << sum0[0] - sum0[1] << std::endl;
                        }
                    }
                }
                else if (nspin == 4) // SOC
                {
                    const int m = 2 * l + 1;
                    std::vector<std::vector<double>> blocks;
                    extract_soc_pauli_blocks(dftu, iat, l, soc_layout, blocks);

                    if (readable)
                    {
                        // Compact snapshot layout aligned with the
                        // nspin == 2 style, e.g.
                        // " Fe Atom 1 L 2 mag mx my mz"
                        // " spin 1 nelec ..."
                        // followed by the four spin blocks.
                        double mag_x = 0.0;
                        double mag_y = 0.0;
                        double mag_z = 0.0;
                        for (int mm = 0; mm < m; ++mm)
                        {
                            const int k = mm * m + mm;
                            mag_x += blocks[1][k];
                            mag_y += blocks[2][k];
                            mag_z += blocks[3][k];
                        }
                        // Spin-resolved electron counts: traces of the
                        // up/up and down/down blocks.
                        double n_up = 0.0;
                        double n_down = 0.0;
                        for (int mm = 0; mm < m; ++mm)
                        {
                            const int k = mm * m + mm;
                            n_up += 0.5 * (blocks[0][k] + blocks[3][k]);
                            n_down += 0.5 * (blocks[0][k] - blocks[3][k]);
                        }

                        ofs << "\n " << elem_label;
                        ofs << " Atom " << iat + 1;
                        ofs << " L " << l;
                        if (diag)
                        {
                            ofs << " mag " << std::fixed << std::setprecision(8)
                                << mag_x << " " << mag_y << " " << mag_z;
                        }
                        ofs << std::endl;

                        ofs << std::fixed << std::setprecision(8);

                        // spin 1: spin up/up block (Re).
                        ofs << " spin 1";
                        if (diag)
                        {
                            ofs << " nelec " << n_up;
                        }
                        ofs << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * (blocks[0][k] + blocks[3][k]);
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }

                        // Off-diagonal Re block.
                        ofs << " spin 12 re" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * blocks[1][k];
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }

                        // Off-diagonal Im block.
                        ofs << " spin 12 im" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * blocks[2][k];
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }

                        // spin 2: spin down/down block (Re).
                        ofs << " spin 2";
                        if (diag)
                        {
                            ofs << " nelec " << n_down;
                        }
                        ofs << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * (blocks[0][k] - blocks[3][k]);
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }
                    }
                    else if (diag)
                    {
                        // Occupation numbers of the correlated spinor shell.
                        std::vector<double> eigenvalues = calculate_spinor_eigenvalues(blocks, m);
                        ofs << " Eigenvalues in spinor space:" << std::endl;
                        ofs << std::setprecision(8) << std::fixed;
                        for (int i = 0; i < 2 * m; ++i)
                        {
                            ofs << std::setw(12) << eigenvalues[i];
                            if ((i + 1) % m == 0)
                            {
                                ofs << std::endl;
                            }
                        }

                        // Per-atom charge and spin-moment traces of the shell.
                        double n_tot = 0.0;
                        double mag_x = 0.0;
                        double mag_y = 0.0;
                        double mag_z = 0.0;
                        for (int mm = 0; mm < m; ++mm)
                        {
                            const int k = mm * m + mm;
                            n_tot += blocks[0][k];
                            mag_x += blocks[1][k];
                            mag_y += blocks[2][k];
                            mag_z += blocks[3][k];
                        }
                        ofs << std::fixed << std::setprecision(8);
                        ofs << " Trace (electrons) " << n_tot << std::endl;
                        ofs << " Magnetism for atom " << iat + 1 << " (mx, my, mz):"
                            << std::setw(12) << mag_x << std::setw(12) << mag_y
                            << std::setw(12) << mag_z << std::endl;

                        // Matrix elements use fixed-point notation with 8
                        // fractional digits (~1e-8 absolute precision).
                        ofs << std::fixed << std::setprecision(8);

                        // Hermitian occupation matrix shown spin block by spin block.
                        ofs << " Occupation matrix, spin up/up (Re):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * (blocks[0][k] + blocks[3][k]);
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }
                        ofs << " Occupation matrix, spin up/down (Re):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * blocks[1][k];
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }
                        ofs << " Occupation matrix, spin up/down (Im):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * blocks[2][k];
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }
                        ofs << " Occupation matrix, spin down/down (Re):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * (blocks[0][k] - blocks[3][k]);
                                ofs << " " << std::setw(10) << val;
                            }
                            ofs << std::endl;
                        }
                    }
                    else
                    {
                        // Real 2m x 2m matrix in spin basis; this is the
                        // layout parsed by read_occup_m() on restart.
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ++ipol0)
                            {
                                for (int m1 = 0; m1 < m; ++m1)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ++ipol1)
                                    {
                                        const int k = m0 * m + m1;
                                        double val = 0.0;
                                        if (ipol0 == 0 && ipol1 == 0)
                                        {
                                            val = 0.5 * (blocks[0][k] + blocks[3][k]);
                                        }
                                        else if (ipol0 == 1 && ipol1 == 1)
                                        {
                                            val = 0.5 * (blocks[0][k] - blocks[3][k]);
                                        }
                                        else
                                        {
                                            val = 0.5 * blocks[1][k];
                                        }
                                        ofs << " " << std::setw(10)
                                            << std::setprecision(8)
                                            << std::fixed << val;
                                    }
                                }
                                ofs << std::endl;
                            }
                        }
                    }
                }
            }     // l
        }         // I
    }             // T

    return;
}


} // namespace DFTU_BASE

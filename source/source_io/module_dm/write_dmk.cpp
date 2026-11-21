#include "source_io/module_dm/write_dmk.h"
#include "source_io/module_hs/dense_text_body.h"
#include "source_io/module_output/spin_tag.h"
#include "source_io/module_output/write_notice.h"

#include "source_base/parallel_common.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_cell/ucell_io.h"

#include <cerrno>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
/// @brief Read the next line that carries data.
/// @details Empty lines and lines whose first non-blank character is '#' are
///          skipped, so a reader does not depend on the exact blank-line layout
///          of the writer.
/// @return false when the end of the file is reached.
bool read_data_line(std::ifstream& ifs, std::string& line)
{
    while (std::getline(ifs, line))
    {
        const std::size_t first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos)
        {
            continue;
        }
        if (line[first] == '#')
        {
            continue;
        }
        return true;
    }
    return false;
}

/// @brief Parse the leading integer of a header line.
bool parse_leading_int(const std::string& line, int& value)
{
    std::istringstream iss(line);
    return static_cast<bool>(iss >> value);
}

/// @brief Parse the leading real number of a header line.
bool parse_leading_double(const std::string& line, double& value)
{
    std::istringstream iss(line);
    return static_cast<bool>(iss >> value);
}

/// @brief Parse the three leading real numbers of a header line into a vector.
bool parse_leading_vec3(const std::string& line, ModuleBase::Vector3<double>& vec)
{
    std::istringstream iss(line);
    return static_cast<bool>(iss >> vec.x >> vec.y >> vec.z);
}

/// @brief Compare two k-point coordinates in direct space, allowing integer
///        offsets (i.e. differing by a reciprocal lattice vector).
/// @details The density matrix is invariant under k -> k + G because the
///          phase factor e^{ik·R} is unchanged; only a mismatch modulo the
///          reciprocal lattice is a genuine inconsistency.
bool same_kvec_direct(const ModuleBase::Vector3<double>& a,
                      const ModuleBase::Vector3<double>& b,
                      const double tol = 1e-6)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::abs(dx - std::round(dx)) < tol
        && std::abs(dy - std::round(dy)) < tol
        && std::abs(dz - std::round(dz)) < tol;
}

/// @brief Advance the stream to the frame opened by "Ionic_Step target_step".
/// @details Appended multi-step files (out_app_flag == 1) concatenate one
///          record per ionic step, each opened by a marker line, so reading a
///          given step requires locating its marker first.
/// @return false when the requested frame is absent.
bool seek_frame(std::ifstream& ifs, const int target_step)
{
    const std::string marker = "Ionic_Step";
    std::string line;
    while (std::getline(ifs, line))
    {
        const std::size_t pos = line.find(marker);
        if (pos == std::string::npos)
        {
            continue;
        }
        std::istringstream iss(line.substr(pos + marker.size()));
        int step = 0;
        if (iss >> step && step == target_step)
        {
            return true;
        }
    }
    return false;
}
} // namespace

std::string ModuleIO::dmk_gen_fname(const bool gamma_only, const int ispin, const int nspin, const int ik, const int istep)
{
    // set istep = -1 if you don't want the 'g' index appears in the file name
    assert(istep>=-1);

    // ik should be the correct one

    std::string fname = "dm";

    if (!gamma_only)
    {
        fname += "k" + std::to_string(ik + 1);
    }

    if (nspin == 2)
    {
        fname += "s" + std::to_string(ispin + 1);
    }

    if( istep >= 0 )
    {
        fname += "g" + std::to_string(istep + 1);
    }

    fname += "_nao.txt";

    return fname;
}



bool ModuleIO::dmk_readData(std::ifstream& ifs, double& data)
{
    if (!(ifs >> data))
    {
        return false;
    }
    if (ifs.fail())
    {
        return false;
    }
    // The stream extractor accepts the literals "inf"/"nan" without setting
    // failbit; a restart density must stay finite, so reject them here.
    return std::isfinite(data);
}

bool ModuleIO::dmk_readData(std::ifstream& ifs, std::complex<double>& data)
{
    std::string complex_str;
    if (!(ifs >> complex_str))
    {
        return false;
    }

    // A body element has the fixed form "(real,imag)".
    if (complex_str.size() < 5 || complex_str.front() != '(' || complex_str.back() != ')')
    {
        return false;
    }
    const std::size_t comma_pos = complex_str.find(',');
    if (comma_pos == std::string::npos || comma_pos == 1 || comma_pos == complex_str.size() - 2)
    {
        return false;
    }

    const std::string real_str = complex_str.substr(1, comma_pos - 1);
    const std::string imag_str = complex_str.substr(comma_pos + 1, complex_str.size() - comma_pos - 2);

    // strtod consumes the full substring on overflow (e.g. "1e999") and only
    // reports the range error through errno, so the end-pointer check alone
    // would accept an infinite component; errno must be checked explicitly.
    errno = 0;
    char* end_ptr = nullptr;
    const double real_v = std::strtod(real_str.c_str(), &end_ptr);
    if (end_ptr == real_str.c_str() || *end_ptr != '\0' || errno == ERANGE)
    {
        return false;
    }
    errno = 0;
    const double imag_v = std::strtod(imag_str.c_str(), &end_ptr);
    if (end_ptr == imag_str.c_str() || *end_ptr != '\0' || errno == ERANGE)
    {
        return false;
    }
    // Also reject the literal "inf"/"nan" spellings that strtod parses
    // without setting errno: a density matrix component must be finite.
    if (!std::isfinite(real_v) || !std::isfinite(imag_v))
    {
        return false;
    }

    data = std::complex<double>(real_v, imag_v);
    return true;
}

template <typename T>
bool ModuleIO::read_dmk(const int nspin,
        const int nk,
        const K_Vectors &kv,
        const Parallel_Blkcyclic& pv,
        const bool column_major,
        const std::string& dmk_dir,
        const int istep,
        const int out_app_flag,
        std::vector<std::vector<T>>& dmk,
        std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "read_dmk");
    ModuleBase::timer::start("ModuleIO", "read_dmk");

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &my_rank);
#endif

    const int nlocal = pv.nrglb();
    const bool gamma_only = std::is_same<double, T>::value;
    std::vector<T> dmk_global(my_rank == 0 ? nlocal * nlocal : 0);
    std::vector<std::vector<T>> dmk_read(
        nspin * nk, std::vector<T>(pv.nrloc() * pv.ncloc(), T(0)));

    // write_dmk emits an optional "Ionic_Step" frame marker followed by
    // eleven header lines: nspin, spin index, nkstot_nospin, nkstot, ik,
    // kvec_c, kvec_d, weight, Fermi energy, nlocal and the matrix size.
    // Reading them in order also rejects files of a different layout,
    // spin channel or k-mesh instead of silently mismatching them.
    const int n_header = 11;
    std::vector<std::string> header(n_header);
    std::string line;

#ifdef __MPI
    Parallel_Blkcyclic pv_glb;
    pv_glb.set(nlocal, nlocal, nlocal, pv.ictxt());
#endif

    bool read_success = true;
    // Mode 0 rewrites every k/spin file under the same fixed name at each
    // ionic step, so an interrupted write can leave files from different
    // steps on disk. The first parsed marker becomes the reference every
    // other file has to match; files without a marker (older writer) skip
    // the cross-check.
    int common_step = -1;
    for (int ispin = 0; ispin < nspin && read_success; ispin++)
    {
        for (int ik = 0; ik < nk && read_success; ik++)
        {
            const std::size_t index_k = static_cast<std::size_t>(ik) + static_cast<std::size_t>(nk) * ispin;
            if (my_rank == 0)
            {
                // Keep file-local error exits inside this scope so every rank
                // reaches the cooperative status broadcast below.
                do
                {
                    // Mirror the writer: only mode 2 embeds the ionic step in the
                    // filename, modes 0 and 1 reuse a fixed name.
                    const int istep_file = (out_app_flag == 2) ? istep : -1;
                    const std::string fn = dmk_dir + dmk_gen_fname(
                                                         gamma_only,
                                                         ispin,
                                                         nspin,
                                                         ik,
                                                         istep_file);
                    std::ifstream ifs(fn.c_str());

                    if (!ifs)
                    {
                        ofs_running << " Cannot find density matrix file " << fn << " for k-point " << ik + 1
                                    << std::endl;
                        ModuleBase::WARNING("ModuleIO::read_dmk", "Can't open density matrix (k) file < " + fn + " >.");
                        read_success = false;
                        break;
                    }
                    ofs_running << " Read density matrix file " << fn << " for k-point " << ik + 1 << std::endl;

                    // A fixed-name file of a multi-step run (mode 1) holds one
                    // record per ionic step, so the requested frame has to be
                    // located first; a negative step keeps the first frame.
                    if (out_app_flag == 1 && istep >= 0)
                    {
                        if (!seek_frame(ifs, istep + 1))
                        {
                            ofs_running << " Cannot find ionic step " << istep + 1 << " in density matrix file " << fn
                                        << std::endl;
                            read_success = false;
                            break;
                        }
                    }

                    // ---- header ----
                    int frame_step = -1;
                    if (!read_data_line(ifs, line))
                    {
                        ofs_running << " Empty density matrix file " << fn << std::endl;
                        read_success = false;
                        break;
                    }
                    if (line.find("Ionic_Step") != std::string::npos)
                    {
                        std::istringstream iss_step(line.substr(line.find("Ionic_Step") + 10));
                        iss_step >> frame_step;
                        if (!read_data_line(ifs, line))
                        {
                            ofs_running << " Truncated density matrix file " << fn << std::endl;
                            read_success = false;
                            break;
                        }
                    }
                    // A per-step file (mode 2) must describe the step its name
                    // announces, otherwise the restart would silently mix steps.
                    if (out_app_flag == 2 && frame_step >= 0 && frame_step != istep + 1)
                    {
                        ofs_running << " Ionic step mismatch in density matrix file " << fn << ": file says "
                                    << frame_step << ", expected " << istep + 1 << std::endl;
                        read_success = false;
                        break;
                    }
                    // Mode 0 reuses fixed names, so the only cross-file guard
                    // is that every file carries the same step marker.
                    if (out_app_flag == 0 && frame_step >= 0)
                    {
                        if (common_step < 0)
                        {
                            common_step = frame_step;
                        }
                        else if (frame_step != common_step)
                        {
                            ofs_running << " Inconsistent ionic step in density matrix file " << fn
                                        << ": file says " << frame_step << ", expected " << common_step
                                        << " (from the first file)" << std::endl;
                            read_success = false;
                            break;
                        }
                    }
                    header[0] = line;
                    for (int ih = 1; ih < n_header; ++ih)
                    {
                        if (!read_data_line(ifs, header[ih]))
                        {
                            ofs_running << " Truncated density matrix file " << fn << std::endl;
                            read_success = false;
                            break;
                        }
                    }
                    if (!read_success)
                    {
                        break;
                    }

                    int nspin_file = 0;
                    int ispin_file = 0;
                    int nkstot_nospin_file = 0;
                    int nkstot_file = 0;
                    int ik_file = 0;
                    double ef_file = 0.0;
                    int nlocal_file = 0;
                    int nrow_file = 0;
                    int ncol_file = 0;
                    bool header_ok = parse_leading_int(header[0], nspin_file);
                    header_ok = header_ok && parse_leading_int(header[1], ispin_file);
                    header_ok = header_ok && parse_leading_int(header[2], nkstot_nospin_file);
                    header_ok = header_ok && parse_leading_int(header[3], nkstot_file);
                    header_ok = header_ok && parse_leading_int(header[4], ik_file);
                    header_ok = header_ok && parse_leading_double(header[8], ef_file);
                    header_ok = header_ok && parse_leading_int(header[9], nlocal_file);
                    {
                        std::istringstream iss(header[10]);
                        header_ok = header_ok && static_cast<bool>(iss >> nrow_file >> ncol_file);
                    }
                    if (!header_ok)
                    {
                        ofs_running << " Malformed header in density matrix file " << fn << std::endl;
                        read_success = false;
                        break;
                    }
                    if (nspin_file != nspin || ispin_file != ispin + 1)
                    {
                        ofs_running << " Spin mismatch in density matrix file " << fn
                                    << ": file says nspin = " << nspin_file << ", spin index = " << ispin_file
                                    << "; expected nspin = " << nspin << ", spin index = " << ispin + 1 << std::endl;
                        read_success = false;
                        break;
                    }
                    if (nkstot_nospin_file != kv.get_nkstot_nospin() || nkstot_file != kv.get_nkstot() ||
                        ik_file != ik + 1)
                    {
                        ofs_running << " k-point mismatch in density matrix file " << fn
                                    << ": file says nkstot = " << nkstot_file << " (nospin " << nkstot_nospin_file
                                    << "), ik = " << ik_file << "; expected nkstot = " << kv.get_nkstot() << " (nospin "
                                    << kv.get_nkstot_nospin() << "), ik = " << ik + 1 << std::endl;
                        read_success = false;
                        break;
                    }

                    // The counts and the index alone do not pin down the k-point:
                    // a reordered or shifted mesh of the same size would pass the
                    // check above but describe a different density. The stored
                    // direct coordinates must therefore match the current mesh
                    // (modulo a reciprocal lattice vector, which leaves the
                    // phase e^{ik·R} unchanged).
                    ModuleBase::Vector3<double> kvec_d_file;
                    if (!parse_leading_vec3(header[6], kvec_d_file))
                    {
                        ofs_running << " Malformed k-point coordinates in density matrix file " << fn << std::endl;
                        read_success = false;
                        break;
                    }
                    if (!same_kvec_direct(kvec_d_file, kv.kvec_d[ik]))
                    {
                        ofs_running << " k-point coordinate mismatch in density matrix file " << fn
                                    << ": file says kvec_d = (" << kvec_d_file.x << ", " << kvec_d_file.y << ", "
                                    << kvec_d_file.z << "); expected (" << kv.kvec_d[ik].x << ", "
                                    << kv.kvec_d[ik].y << ", " << kv.kvec_d[ik].z << ")" << std::endl;
                        read_success = false;
                        break;
                    }
                    if (nlocal_file != nlocal || nrow_file != nlocal || ncol_file != nlocal)
                    {
                        ofs_running << " Basis size mismatch in density matrix file " << fn << ": file says "
                                    << nlocal_file << " basis and a " << nrow_file << " x " << ncol_file
                                    << " matrix, expected " << nlocal << std::endl;
                        read_success = false;
                        break;
                    }
                    ofs_running << " Fermi energy of " << fn << " is " << ef_file << " Ry" << std::endl;

                    // ---- UnitCell block between the header and the matrix body ----
                    // Consume the blank separator line(s) first, then let the
                    // UnitCell reader skip its fixed set of lines.
                    while (true)
                    {
                        const std::streampos pos = ifs.tellg();
                        if (!std::getline(ifs, line))
                        {
                            break;
                        }
                        if (line.find_first_not_of(" \t\r") != std::string::npos)
                        {
                            ifs.seekg(pos);
                            break;
                        }
                    }
                    if (!ModuleIO::UcellIO::read_ucell(ifs))
                    {
                        ofs_running << " Malformed UnitCell block in density matrix file " << fn
                                    << ": could not parse the atom-count line or ran out of data"
                                    << std::endl;
                        read_success = false;
                        break;
                    }

                    // ---- matrix body, always row-major in the file ----
                    for (int i = 0; i < nlocal && read_success; ++i)
                    {
                        for (int j = 0; j < nlocal; ++j)
                        {
                            T value = T(0);
                            if (!dmk_readData(ifs, value))
                            {
                                read_success = false;
                                break;
                            }
                            dmk_global[ModuleIO::dense_buffer_index(
                                i,
                                j,
                                nlocal,
                                column_major)] = value;
                        }
                    }
                    if (!read_success)
                    {
                        ofs_running << " Failed to read the matrix body of density matrix file " << fn << std::endl;
                        break;
                    }
                    ifs.close();
                } while (false);
            } // rank0

#ifdef __MPI
            // Cooperative failure: every rank distributes the same matrix or
            // stops together, so a rank-0 parse error cannot deadlock the run.
            Parallel_Common::bcast_bool(read_success, pv.comm());
#endif
            if (!read_success)
            {
                break;
            }

#ifdef __MPI
            Cpxgemr2d(
                nlocal,
                nlocal,
                dmk_global.data(),
                1,
                1,
                const_cast<int*>(pv_glb.desc()),
                dmk_read[index_k].data(),
                1,
                1,
                const_cast<int*>(pv.desc()),
                pv_glb.ictxt());
#else
            dmk_read[index_k] = dmk_global;
#endif
        } // ik
    } // ispin

    if (read_success)
    {
        dmk.swap(dmk_read);
    }
    ModuleBase::timer::end("ModuleIO", "read_dmk");
    return read_success;
}

template <typename T>
void ModuleIO::write_dmk(const std::vector<std::vector<T>>& dmk,
        const K_Vectors &kv,
        const int precision,
        const std::vector<double>& efs,
        const UnitCell* ucell,
        const Parallel_Blkcyclic& pv,
        const bool column_major,
        const std::string& dmk_dir,
        const int istep,
        const int out_app_flag,
        std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_dmk");
    ModuleBase::timer::start("ModuleIO", "write_dmk");

    // Modes 0 and 1 reuse a fixed filename (mode 0 overwrites the previous
    // snapshot, mode 1 appends a new record); only mode 2 embeds the ionic
    // step in the filename. The file body always records the real ionic step.
    const int istep_file = (out_app_flag == 2) ? istep : -1;
    // Each multi-step run (MD or relax) appends a new record to the fixed
    // file, mirroring the open_sparse_file convention; single-point SCF
    // passes istep = -1 and starts a fresh file. The first written step
    // (istep == 0) truncates so a fresh run drops records left by a previous
    // run; only later steps (istep > 0) append.
    const bool append_on_restart = (out_app_flag == 1) && (istep > 0);

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &my_rank);
#endif

    bool gamma_only = std::is_same<double, T>::value;
    const int nlocal = pv.nrglb();
    const int nspin = efs.size();
    assert(nspin > 0);
    const int nk = dmk.size() / nspin;

    if (nk * nspin != dmk.size())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dmk", "The size of dmk is not consistent with nspin and nk.");
    }

    Parallel_Blkcyclic pv_glb;
    std::vector<T> dmk_global(my_rank == 0 ? nlocal * nlocal : 0);

    for (int ispin = 0; ispin < nspin; ispin++)
    {
        for (int ik = 0; ik < nk; ik++)
        {
            // gather dmk[ik] to the reused global buffer
#ifdef __MPI
            pv_glb.set(nlocal, nlocal, nlocal, pv.ictxt());
            Cpxgemr2d(
            nlocal,
            nlocal,
            const_cast<T*>(dmk[ik + nk * ispin].data()),
            1,
            1,
            const_cast<int*>(pv.desc()),
            dmk_global.data(),
            1,
            1,
            const_cast<int*>(pv_glb.desc()),
            pv_glb.ictxt());
#else
            dmk_global = dmk[ik + nk * ispin];
#endif

            if (my_rank == 0)
            {
                std::string fn = dmk_dir
                    + dmk_gen_fname(gamma_only, ispin, nspin, ik, istep_file);

                std::ofstream ofs;
                std::ios_base::openmode mode = std::ios::out;
                if (append_on_restart)
                {
                    mode |= std::ios::app;
                }
                ofs.open(fn.c_str(), mode);

                ofs << "Ionic_Step " << istep+1 << std::endl;

                if (!ofs)
                {
                    ModuleBase::WARNING("ModuleIO::write_dmk", "Can't create DENSITY MATRIX File < " + fn + " >.");
                    continue;
                }
                else
                {
                    ModuleIO::write_matrix_notice(ofs_running, "DM(k)", ispin, nspin, fn);
                }


                // information about density matrix at this k-point
                ofs << " " << nspin << " # number of spin directions" << std::endl;
                ofs << " " << ispin+1 << " # spin index" << std::endl;
                ofs << " " << kv.get_nkstot_nospin() << " # total k points " << std::endl;
                ofs << " " << kv.get_nkstot() << " # total k points after symmetrized (if open) " << std::endl;
                ofs << " " << ik+1 << " # k-point index " << std::endl;
                ofs << " " << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z
                    << " # k point coordinate (Cartesian) " << std::endl;
                ofs << " " << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y << " " << kv.kvec_d[ik].z
                    << " # k point coordinate (direct) " << std::endl;
                ofs << " " << kv.wk[ik] << " # weight of this k point" << std::endl;
                ofs << " " << efs[ispin] << " # Fermi energy in Ry " << std::endl;
                // State the element convention in the file itself: dm_from_psi
                // builds DM = sum_n w_n conj(c_n) c_n^T (conj on the first
                // index), i.e. the complex conjugate of the physical one-body
                // density matrix. This line lives inside the header block so
                // read_data_line (which skips '#' lines) absorbs it.
                ofs << "# matrix element (i,j) = sum_n w_n * conj(c_n,i) * c_n,j" << std::endl;
                ofs << " " << nlocal << " # number of localized basis " << std::endl;
                ofs << " " << nlocal << " " << nlocal << " # size of this matrix " << std::endl;
                ofs << std::endl;

                // write ucell
                ModuleIO::UcellIO::write_ucell(ofs, ucell);

                // The file body is row-major; column_major selects the
                // traversal order of the gathered (possibly column-major)
                // DM(k) buffer so the transpose is not emitted.
                ModuleIO::write_dense_text_body(ofs, dmk_global, nlocal, precision, column_major);
                ofs.close();
            } // rank0
        }     // ik
    }         // ispin

    ModuleBase::timer::end("ModuleIO", "write_dmk");
}

template bool ModuleIO::read_dmk<double>(const int nspin,
        const int nk,
        const K_Vectors &kv,
        const Parallel_Blkcyclic& pv,
        const bool column_major,
        const std::string& dmk_dir,
        const int istep,
        const int out_app_flag,
        std::vector<std::vector<double>>& dmk,
        std::ofstream &ofs);

template bool ModuleIO::read_dmk<std::complex<double>>(const int nspin,
        const int nk,
        const K_Vectors &kv,
        const Parallel_Blkcyclic& pv,
        const bool column_major,
        const std::string& dmk_dir,
        const int istep,
        const int out_app_flag,
        std::vector<std::vector<std::complex<double>>>& dmk,
        std::ofstream &ofs);

template void ModuleIO::write_dmk<double>(const std::vector<std::vector<double>>& dmk,
        const K_Vectors &kv,
        const int precision,
        const std::vector<double>& efs,
        const UnitCell* ucell,
        const Parallel_Blkcyclic& pv,
        const bool column_major,
        const std::string& dmk_dir,
        const int istep,
        const int out_app_flag,
        std::ofstream& ofs_running);

template void ModuleIO::write_dmk<std::complex<double>>(const std::vector<std::vector<std::complex<double>>>& dmk,
        const K_Vectors &kv,
        const int precision,
        const std::vector<double>& efs,
        const UnitCell* ucell,
        const Parallel_Blkcyclic& pv,
        const bool column_major,
        const std::string& dmk_dir,
        const int istep,
        const int out_app_flag,
        std::ofstream& ofs_running);


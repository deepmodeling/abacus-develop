#include "write_hs_r.h"

#include "source_base/module_out/sparse_matrix.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/lcao_hs_arrays.hpp"
#include "source_lcao/spar_dh.h"
#include "source_lcao/spar_hsr.h"
#include "source_lcao/spar_st.h"
#include "write_hs_sparse.h"

#include <algorithm>
#include <complex>
#include <fstream>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix.
// If the absolute value of the matrix element is less than or equal to the
// 'sparse_thr', it will be ignored.

void ModuleIO::output_dSR(const int& istep,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr,
                          const int precision)
{
    ModuleBase::TITLE("ModuleIO", "output_dSR");
    ModuleBase::timer::start("ModuleIO", "output_dSR");

    sparse_format::cal_dS(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);

    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary, "s", precision);

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_dSR");
    return;
}

void ModuleIO::output_dHR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr,
                          const int precision)
{
    ModuleBase::TITLE("ModuleIO", "output_dHR");
    ModuleBase::timer::start("ModuleIO", "output_dHR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                         #Print out dH/dR#                          |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    const int nspin = PARAM.inp.nspin;

    if (nspin == 1 || nspin == 4)
    {
        // mohan add 2024-04-01
        const int cspin = 0;

        sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff);
    }
    else if (nspin == 2)
    {
        for (int cspin = 0; cspin < 2; cspin++)
        {
            sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff);
        }
    }
    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary, "h", precision);

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_dHR");
    return;
}

template <typename TK>
void ModuleIO::output_SR(Parallel_Orbitals& pv,
                         const Grid_Driver& grid,
                         hamilt::Hamilt<TK>* p_ham,
                         const std::string& SR_filename,
                         const bool& binary,
                         const double& sparse_thr,
                         const int precision)
{
    ModuleBase::TITLE("ModuleIO", "output_SR");
    ModuleBase::timer::start("ModuleIO", "output_SR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                 #Print out overlap matrix S(R)#                    |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::cout << " Overlap matrix file is in " << SR_filename << std::endl;
    GlobalV::ofs_running << " Overlap matrix file is in " << SR_filename << std::endl;

    LCAO_HS_Arrays HS_Arrays;

    sparse_format::cal_SR(pv,
                          HS_Arrays.all_R_coor,
                          HS_Arrays.SR_sparse,
                          HS_Arrays.SR_soc_sparse,
                          grid,
                          sparse_thr,
                          p_ham);

    const int istep = 0;
    ModuleIO::SparseWriteOptions options;
    options.filename = SR_filename;
    options.label = "S";
    options.threshold = sparse_thr;
    options.binary = binary;
    options.precision = precision;
    options.istep = istep;
    options.reduce = true;
    options.temp_dir = PARAM.globalv.global_out_dir;

    if (PARAM.inp.nspin == 4)
    {
        ModuleIO::save_sparse(HS_Arrays.SR_soc_sparse,
                              HS_Arrays.all_R_coor,
                              pv,
                              options);
    }
    else
    {
        ModuleIO::save_sparse(HS_Arrays.SR_sparse,
                              HS_Arrays.all_R_coor,
                              pv,
                              options);
    }

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_SR");
    return;
}

void ModuleIO::output_TR(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         LCAO_HS_Arrays& HS_Arrays,
                         const Grid_Driver& grid,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const std::string& TR_filename,
                         const bool& binary,
                         const double& sparse_thr,
                         const int precision)
{
    ModuleBase::TITLE("ModuleIO", "output_TR");
    ModuleBase::timer::start("ModuleIO", "output_TR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |           #Print out kinetic energy term matrix T(R)#              |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::stringstream sst;
    if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
    {
        sst << PARAM.globalv.global_matrix_dir << TR_filename << "g" << istep;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }
    else
    {
        sst << PARAM.globalv.global_out_dir << TR_filename;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }

    sparse_format::cal_TR(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);
    ModuleIO::SparseWriteOptions options;
    options.filename = sst.str();
    options.label = "T";
    options.threshold = sparse_thr;
    options.binary = binary;
    options.precision = precision;
    options.istep = istep;
    options.reduce = true;
    options.temp_dir = PARAM.globalv.global_out_dir;

    ModuleIO::save_sparse(HS_Arrays.TR_sparse,
                          HS_Arrays.all_R_coor,
                          pv,
                          options);

    sparse_format::destroy_T_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_TR");
    return;
}

template void ModuleIO::output_SR<double>(Parallel_Orbitals& pv,
                                          const Grid_Driver& grid,
                                          hamilt::Hamilt<double>* p_ham,
                                          const std::string& SR_filename,
                                          const bool& binary,
                                          const double& sparse_thr,
                                          const int precision);
template void ModuleIO::output_SR<std::complex<double>>(Parallel_Orbitals& pv,
                                                        const Grid_Driver& grid,
                                                        hamilt::Hamilt<std::complex<double>>* p_ham,
                                                        const std::string& SR_filename,
                                                        const bool& binary,
                                                        const double& sparse_thr,
                                                        const int precision);

#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_hamilt/module_hcontainer/output_hcontainer.h"
#include "source_cell/ucell_io.h"

std::string ModuleIO::hsr_gen_fname(const std::string& prefix,
                                     const int ispin,
                                     const bool append,
                                     const int istep)
{
    return hsr_gen_fname(prefix, ispin, append, istep, 1);
}

std::string ModuleIO::hsr_gen_fname(const std::string& prefix,
                                     const int ispin,
                                     const bool append,
                                     const int istep,
                                     const int out_type)
{
    const std::string extension = out_type == 2 ? ".dat" : ".csr";
    if (!append && istep >= 0)
    {
        return prefix + std::to_string(ispin + 1) + "g" + std::to_string(istep + 1) + "_nao" + extension;
    }
    else
    {
        return prefix + std::to_string(ispin + 1) + "_nao" + extension;
    }
}

std::string ModuleIO::sr_gen_fname(const bool append, const int istep)
{
    return sr_gen_fname(append, istep, 1);
}

std::string ModuleIO::sr_gen_fname(const bool append, const int istep, const int out_type)
{
    const std::string extension = out_type == 2 ? ".dat" : ".csr";
    if (!append && istep >= 0)
    {
        return "srg" + std::to_string(istep + 1) + "_nao" + extension;
    }
    return "sr_nao" + extension;
}

std::string ModuleIO::dhr_gen_fname(const std::string& prefix,
                                     const int ispin,
                                     const bool append,
                                     const int istep)
{
    std::string fname = prefix + "rs" + std::to_string(ispin + 1);
    if (!append && istep >= 0)
    {
        fname += "g" + std::to_string(istep + 1);
    }
    fname += "_nao.csr";
    return fname;
}

template <typename TR>
void ModuleIO::write_hcontainer_csr(const std::string& fname,
                                     const UnitCell* ucell,
                                     const int precision,
                                     hamilt::HContainer<TR>* mat_serial,
                                     const int istep,
                                     const int ispin,
                                     const int nspin,
                                     const std::string& label,
                                     const std::string& representation_note)
{
    std::ofstream ofs;
    if (istep <= 0)
    {
        ofs.open(fname);
    }
    else
    {
        ofs.open(fname, std::ios::app);
    }
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_hcontainer_csr",
                                 "Cannot open HContainer CSR file: " + fname);
    }

    ofs << " --- Ionic Step " << istep + 1 << " ---" << std::endl;
    ofs << " # print " << label << " matrix in real space " << label << "(R)" << std::endl;
    ofs << " " << nspin << " # number of spin directions" << std::endl;
    ofs << " " << ispin + 1 << " # spin index" << std::endl;
    ofs << " " << mat_serial->get_nbasis() << " # number of localized basis" << std::endl;
    ofs << " " << mat_serial->size_R_loop() << " # number of Bravais lattice vector R" << std::endl;
    ofs << std::endl;

    ModuleIO::UcellIO::write_ucell(ofs, ucell);
    if (!representation_note.empty())
    {
        ofs << "# representation: " << representation_note << std::endl;
    }
    ofs << std::endl;

    const double sparse_threshold = 1e-10;
    hamilt::Output_HContainer<TR> out(mat_serial, ofs, sparse_threshold, precision);
    out.write();
    ofs.close();
}

namespace
{
template <typename T>
void write_native_value(std::ofstream& ofs, const T& value)
{
    ofs.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

void write_native_value(std::ofstream& ofs, const std::complex<double>& value)
{
    const double real = value.real();
    const double imag = value.imag();
    write_native_value(ofs, real);
    write_native_value(ofs, imag);
}

template <typename TR>
std::vector<std::vector<std::pair<int, int>>> index_atom_pair_R_blocks(
    const hamilt::HContainer<TR>* mat_serial,
    const std::vector<std::tuple<int, int, int>>& R_coordinates)
{
    std::vector<std::vector<std::pair<int, int>>> block_indices(R_coordinates.size());
    std::map<std::tuple<int, int, int>, size_t> R_to_index;
    for (size_t index = 0; index < R_coordinates.size(); ++index)
    {
        R_to_index[R_coordinates[index]] = index;
    }

    for (int iap = 0; iap < mat_serial->size_atom_pairs(); ++iap)
    {
        const auto& atom_pair = mat_serial->get_atom_pair(iap);
        for (size_t ir = 0; ir < atom_pair.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R = atom_pair.get_R_index(ir);
            const std::tuple<int, int, int> coordinate = std::make_tuple(R.x, R.y, R.z);
            const std::map<std::tuple<int, int, int>, size_t>::const_iterator target
                = R_to_index.find(coordinate);
            if (target == R_to_index.end())
            {
                continue;
            }
            block_indices[target->second].push_back(std::make_pair(iap, static_cast<int>(ir)));
        }
    }
    return block_indices;
}

template <typename TR>
ModuleIO::SparseMatrix<TR> make_sparse_R_block(
    const hamilt::HContainer<TR>* mat_serial,
    const std::vector<std::pair<int, int>>& block_indices,
    const double sparse_threshold)
{
    const int nbasis = mat_serial->get_nbasis();
    ModuleIO::SparseMatrix<TR> sparse_matrix(nbasis, nbasis);
    sparse_matrix.setSparseThreshold(sparse_threshold);
    for (const std::pair<int, int>& block_index: block_indices)
    {
        const auto& atom_pair = mat_serial->get_atom_pair(block_index.first);
        const auto matrix_info = atom_pair.get_matrix_values(block_index.second);
        const int* index = std::get<0>(matrix_info).data();
        TR* data = std::get<1>(matrix_info);
        for (int irow = index[0]; irow < index[0] + index[1]; ++irow)
        {
            for (int icol = index[2]; icol < index[2] + index[3]; ++icol)
            {
                sparse_matrix.insert(irow, icol, *data);
                ++data;
            }
        }
    }
    return sparse_matrix;
}
} // namespace

template <typename TR>
void ModuleIO::write_hcontainer_csr_binary(const std::string& fname,
                                            hamilt::HContainer<TR>* mat_serial,
                                            const int istep,
                                            const bool append)
{
    std::ios_base::openmode mode = std::ios::out | std::ios::binary;
    if (append && istep > 0)
    {
        mode |= std::ios::app;
    }
    std::ofstream ofs(fname.c_str(), mode);
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_hcontainer_csr_binary",
                                 "Cannot open HContainer binary CSR file: " + fname);
    }

    const size_t nR_size = mat_serial->size_R_loop();
    if (nR_size > static_cast<size_t>(std::numeric_limits<int>::max()))
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_hcontainer_csr_binary",
                                 "Too many R blocks for native binary CSR output: " + fname);
    }

    std::vector<std::tuple<int, int, int>> R_coordinates;
    R_coordinates.reserve(nR_size);
    for (size_t iR = 0; iR < nR_size; ++iR)
    {
        int rx = 0;
        int ry = 0;
        int rz = 0;
        mat_serial->loop_R(iR, rx, ry, rz);
        R_coordinates.push_back(std::make_tuple(rx, ry, rz));
    }
    std::sort(R_coordinates.begin(), R_coordinates.end());

    const double sparse_threshold = 1e-10;
    const std::vector<std::vector<std::pair<int, int>>> block_indices
        = index_atom_pair_R_blocks(mat_serial, R_coordinates);

    const int step = std::max(istep, 0);
    const int nbasis = mat_serial->get_nbasis();
    const int nR = static_cast<int>(R_coordinates.size());
    write_native_value(ofs, step);
    write_native_value(ofs, nbasis);
    write_native_value(ofs, nR);

    for (size_t R_index = 0; R_index < R_coordinates.size(); ++R_index)
    {
        const auto& R_coordinate = R_coordinates[R_index];
        const int rx = std::get<0>(R_coordinate);
        const int ry = std::get<1>(R_coordinate);
        const int rz = std::get<2>(R_coordinate);
        const SparseMatrix<TR> sparse_matrix
            = make_sparse_R_block(mat_serial, block_indices[R_index], sparse_threshold);
        const auto& elements = sparse_matrix.getElements();
        if (elements.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_hcontainer_csr_binary",
                                     "Too many nonzero values for native binary CSR output: " + fname);
        }
        const int nnz = static_cast<int>(elements.size());

        write_native_value(ofs, rx);
        write_native_value(ofs, ry);
        write_native_value(ofs, rz);
        write_native_value(ofs, nnz);

        for (const auto& element: elements)
        {
            write_native_value(ofs, element.second);
        }

        std::vector<long long> row_ptr(nbasis + 1, 0);
        for (const auto& element: elements)
        {
            write_native_value(ofs, element.first.second);
            ++row_ptr[element.first.first + 1];
        }
        for (int irow = 1; irow <= nbasis; ++irow)
        {
            row_ptr[irow] += row_ptr[irow - 1];
        }
        for (const long long pointer: row_ptr)
        {
            write_native_value(ofs, pointer);
        }
    }

    ofs.close();
    if (!ofs)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_hcontainer_csr_binary",
                                 "Failed to write HContainer binary CSR file: " + fname);
    }
}

template <typename TR>
void ModuleIO::write_hsr(const std::vector<hamilt::HContainer<TR>*>& hr_vec,
                          const hamilt::HContainer<TR>* sr,
                          const UnitCell* ucell,
                          const int out_type,
                          const int precision,
                          const Parallel_2D& paraV,
                          const bool append,
                          const bool gamma_only,
                          const int* iat2iwt,
                          const int nat,
                          const int istep)
{
    if (out_type != 1 && out_type != 2)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_hsr", "out_type must be 1 or 2");
    }
    const int nspin = hr_vec.size();
    assert(nspin > 0);
    const std::string representation_note
        = gamma_only
              ? "gamma-only folded matrix; stored R-space contributions are summed into R = (0, 0, 0)"
              : "";

    // Output HR (one file per spin)
    for (int ispin = 0; ispin < nspin; ispin++)
    {
        const int nbasis = hr_vec[ispin]->get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        hamilt::HContainer<TR> hr_serial(&serialV);
        hamilt::gatherParallels(*hr_vec[ispin], &hr_serial, 0);
#else
        hamilt::HContainer<TR> hr_serial(*hr_vec[ispin]);
#endif

        if (GlobalV::MY_RANK == 0)
        {
            std::string fname = PARAM.globalv.global_out_dir
                                + hsr_gen_fname("hrs", ispin, append, istep, out_type);
            if (out_type == 2)
            {
                write_hcontainer_csr_binary(fname, &hr_serial, istep, append);
            }
            else
            {
                write_hcontainer_csr(
                    fname, ucell, precision, &hr_serial, istep, ispin, nspin, "H", representation_note);
            }
        }
    }

    // Output SR (single file)
    {
        const int nbasis = sr->get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        hamilt::HContainer<TR> sr_serial(&serialV);
        hamilt::gatherParallels(*sr, &sr_serial, 0);
#else
        hamilt::HContainer<TR> sr_serial(*sr);
#endif

        if (GlobalV::MY_RANK == 0)
        {
            std::string fname = PARAM.globalv.global_out_dir
                                + sr_gen_fname(append, istep, out_type);
            if (out_type == 2)
            {
                write_hcontainer_csr_binary(fname, &sr_serial, istep, append);
            }
            else
            {
                write_hcontainer_csr(
                    fname, ucell, precision, &sr_serial, istep, 0, 1, "S", representation_note);
            }
        }
    }
}

// Explicit instantiations
template void ModuleIO::write_hcontainer_csr<double>(
    const std::string&, const UnitCell*, const int,
    hamilt::HContainer<double>*, const int, const int, const int, const std::string&, const std::string&);
template void ModuleIO::write_hcontainer_csr<std::complex<double>>(
    const std::string&, const UnitCell*, const int,
    hamilt::HContainer<std::complex<double>>*, const int, const int, const int, const std::string&, const std::string&);

template void ModuleIO::write_hcontainer_csr_binary<double>(
    const std::string&, hamilt::HContainer<double>*, const int, const bool);
template void ModuleIO::write_hcontainer_csr_binary<std::complex<double>>(
    const std::string&, hamilt::HContainer<std::complex<double>>*, const int, const bool);

template void ModuleIO::write_hsr<double>(
    const std::vector<hamilt::HContainer<double>*>&,
    const hamilt::HContainer<double>*,
    const UnitCell*, const int, const int, const Parallel_2D&,
    const bool, const bool, const int*, const int, const int);
template void ModuleIO::write_hsr<std::complex<double>>(
    const std::vector<hamilt::HContainer<std::complex<double>>*>&,
    const hamilt::HContainer<std::complex<double>>*,
    const UnitCell*, const int, const int, const Parallel_2D&,
    const bool, const bool, const int*, const int, const int);


template <typename TR>
void ModuleIO::write_matrix_r(const std::string& matrix_label,
                               const std::string& description,
                               const std::vector<hamilt::HContainer<TR>*>& matrices,
                               const UnitCell* ucell,
                               const int precision,
                               const Parallel_2D& paraV,
                               const bool append,
                               const int* iat2iwt,
                               const int nat,
                               const int istep)
{
    const int nspin = matrices.size();
    assert(nspin > 0);
    
    for (int ispin = 0; ispin < nspin; ispin++)
    {
        const int nbasis = matrices[ispin]->get_nbasis();
        
        // Generate filename
        std::string fname = dhr_gen_fname(matrix_label, ispin, append, istep);
        if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
        {
            fname = PARAM.globalv.global_matrix_dir + fname;
        }
        else
        {
            fname = PARAM.globalv.global_out_dir + fname;
        }
        
        // Gather parallel matrix to serial
#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        
        hamilt::HContainer<TR> matrix_serial(&serialV);
        hamilt::gatherParallels(*matrices[ispin], &matrix_serial, 0);
        
        if (GlobalV::MY_RANK == 0)
        {
            write_hcontainer_csr(
                fname, ucell, precision, &matrix_serial, istep, ispin, nspin, description, "");
        }
#else
        write_hcontainer_csr(
            fname, ucell, precision, matrices[ispin], istep, ispin, nspin, description, "");
#endif
    }
}

// Template instantiations
template void ModuleIO::write_matrix_r<double>(
    const std::string&,
    const std::string&,
    const std::vector<hamilt::HContainer<double>*>&,
    const UnitCell*,
    const int,
    const Parallel_2D&,
    const bool,
    const int*,
    const int,
    const int);

template void ModuleIO::write_matrix_r<std::complex<double>>(
    const std::string&,
    const std::string&,
    const std::vector<hamilt::HContainer<std::complex<double>>*>&,
    const UnitCell*,
    const int,
    const Parallel_2D&,
    const bool,
    const int*,
    const int,
    const int);

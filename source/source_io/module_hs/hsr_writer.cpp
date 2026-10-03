#include "hsr_writer.h"

#include "source_base/module_out/sparse_matrix.h"
#include "source_base/tool_quit.h"
#include "source_cell/ucell_io.h"
#include "source_cell/unitcell.h"
#include "source_base/global_variable.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_hamilt/module_hcontainer/output_hcontainer.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <fstream>
#include <limits>
#include <tuple>
#include <vector>

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
ModuleIO::SparseMatrix<TR> make_sparse_R_block(hamilt::HContainer<TR>* mat_serial,
                                               const int rx,
                                               const int ry,
                                               const int rz,
                                               const double sparse_threshold)
{
    const int nbasis = mat_serial->get_nbasis();
    ModuleIO::SparseMatrix<TR> sparse_matrix(nbasis, nbasis);
    sparse_matrix.setSparseThreshold(sparse_threshold);
    mat_serial->fix_R(rx, ry, rz);

    for (int iap = 0; iap < mat_serial->size_atom_pairs(); ++iap)
    {
        const auto atom_pair = mat_serial->get_atom_pair(iap);
        const int r_index = atom_pair.find_R(rx, ry, rz);
        if (r_index < 0)
        {
            continue;
        }
        const auto matrix_info = atom_pair.get_matrix_values(r_index);
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

    mat_serial->unfix_R();
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

    const int step = std::max(istep, 0);
    const int nbasis = mat_serial->get_nbasis();
    const int nR = static_cast<int>(R_coordinates.size());
    write_native_value(ofs, step);
    write_native_value(ofs, nbasis);
    write_native_value(ofs, nR);

    const double sparse_threshold = 1e-10;
    for (const auto& R_coordinate: R_coordinates)
    {
        const int rx = std::get<0>(R_coordinate);
        const int ry = std::get<1>(R_coordinate);
        const int rz = std::get<2>(R_coordinate);
        const SparseMatrix<TR> sparse_matrix
            = make_sparse_R_block(mat_serial, rx, ry, rz, sparse_threshold);
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
                          const int istep,
                          const std::string& global_out_dir)
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
            std::string fname = global_out_dir
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
            std::string fname = global_out_dir
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
    const bool, const bool, const int*, const int, const int,
    const std::string&);
template void ModuleIO::write_hsr<std::complex<double>>(
    const std::vector<hamilt::HContainer<std::complex<double>>*>&,
    const hamilt::HContainer<std::complex<double>>*,
    const UnitCell*, const int, const int, const Parallel_2D&,
    const bool, const bool, const int*, const int, const int,
    const std::string&);

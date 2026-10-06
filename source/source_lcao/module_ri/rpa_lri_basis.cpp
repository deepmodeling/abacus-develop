#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace RpaLriDetail
{
std::vector<std::vector<int>> collect_abfs_l_nchi(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    std::vector<std::vector<int>> abfs_l_nchi;
    abfs_l_nchi.reserve(abfs.size());
    for (const auto& abfs_type: abfs)
    {
        std::vector<int> shell_counts;
        shell_counts.reserve(abfs_type.size());
        for (const auto& abfs_l: abfs_type)
        {
            shell_counts.push_back(static_cast<int>(abfs_l.size()));
        }
        abfs_l_nchi.push_back(std::move(shell_counts));
    }
    return abfs_l_nchi;
}

std::vector<std::vector<int>> collect_wfc_l_nchi(const UnitCell& ucell)
{
    std::vector<std::vector<int>> wfc_l_nchi;
    wfc_l_nchi.reserve(static_cast<std::size_t>(ucell.ntype));
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const auto& atom = ucell.atoms[itype];
        if (atom.nwl < 0 || atom.l_nchi.size() < static_cast<std::size_t>(atom.nwl + 1))
        {
            throw std::runtime_error("LibRPA v1 basis output found inconsistent AO shell counts.");
        }
        std::vector<int> shell_counts;
        shell_counts.reserve(static_cast<std::size_t>(atom.nwl + 1));
        for (int l = 0; l <= atom.nwl; ++l)
        {
            if (atom.l_nchi[static_cast<std::size_t>(l)] < 0)
            {
                throw std::runtime_error("LibRPA v1 basis output found negative AO shell count.");
            }
            shell_counts.push_back(atom.l_nchi[static_cast<std::size_t>(l)]);
        }
        wfc_l_nchi.push_back(std::move(shell_counts));
    }
    return wfc_l_nchi;
}

int basis_size_from_shell_counts(const std::vector<int>& shell_counts)
{
    int basis_size = 0;
    for (std::size_t l = 0; l < shell_counts.size(); ++l)
    {
        const int count = shell_counts[l];
        if (count < 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found negative shell count.");
        }
        const int shell_size = static_cast<int>(2 * l + 1);
        if (count > 0 && shell_size > std::numeric_limits<int>::max() / count)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 shell sizes.");
        }
        const int increment = count * shell_size;
        if (increment > std::numeric_limits<int>::max() - basis_size)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 shell sizes.");
        }
        basis_size += increment;
    }
    return basis_size;
}

void write_librpa_split_basis_file(const UnitCell& ucell,
                                   const std::vector<int>& type_sizes,
                                   const std::vector<std::vector<int>>& type_shell_counts,
                                   const std::string& filename)
{
    if (type_sizes.size() != static_cast<std::size_t>(ucell.ntype)
        || type_shell_counts.size() != static_cast<std::size_t>(ucell.ntype))
    {
        throw std::runtime_error("LibRPA v1 basis output found inconsistent atom-type count.");
    }

    int total_basis = 0;
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const int type_size = type_sizes[static_cast<std::size_t>(itype)];
        if (type_size <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type basis size.");
        }
        const int shell_size = basis_size_from_shell_counts(type_shell_counts[static_cast<std::size_t>(itype)]);
        if (type_size != shell_size)
        {
            throw std::runtime_error("LibRPA v1 basis output found inconsistent shell layout size.");
        }
        if (ucell.atoms[itype].na < 0 || type_size > std::numeric_limits<int>::max() / ucell.atoms[itype].na)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        const int increment = type_size * ucell.atoms[itype].na;
        if (increment > std::numeric_limits<int>::max() - total_basis)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        total_basis += increment;
    }

    std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }
    ofs << std::setw(10) << ucell.ntype << std::setw(10) << total_basis << "    abacus" << std::endl;
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        ofs << std::setw(10) << itype + 1 << std::setw(10) << type_sizes[static_cast<std::size_t>(itype)] << std::endl;
    }
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const auto& shell_counts = type_shell_counts[static_cast<std::size_t>(itype)];
        int nshell = 0;
        for (const int count: shell_counts)
        {
            nshell += count;
        }
        ofs << std::setw(10) << itype + 1 << std::setw(10) << nshell << std::endl;
        for (std::size_t l = 0; l < shell_counts.size(); ++l)
        {
            for (int iradial = 0; iradial < shell_counts[l]; ++iradial)
            {
                ofs << std::setw(10) << l << std::endl;
            }
        }
    }
}
} // namespace RpaLriDetail

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_librpa_basis_v1(const UnitCell& ucell,
                                            Exx_LRI<double>* exx_lri,
                                            const std::string& aux_filename,
                                            const std::string& legacy_filename)
try
{
    if (this->runtime.rank != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_librpa_basis_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> type_naux(static_cast<std::size_t>(ucell.ntype), 0);
    std::vector<int> type_nw(static_cast<std::size_t>(ucell.ntype), 0);
    int total_wfc = 0;
    int total_aux = 0;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        type_naux[static_cast<std::size_t>(it)] = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(it);
        if (type_naux[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type auxiliary size.");
        }
        type_nw[static_cast<std::size_t>(it)] = ucell.atoms[it].nw;
        if (type_nw[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type wave-function size.");
        }
        total_wfc += type_nw[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
        total_aux += type_naux[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
    }

    const auto wfc_l_nchi = RpaLriDetail::collect_wfc_l_nchi(ucell);
    const auto aux_l_nchi = RpaLriDetail::collect_abfs_l_nchi(exx_lri->abfs);
    const std::string wfc_filename = this->runtime.input.rpa_outdir + "wfc_basis.txt";
    const std::string aux_output_filename = this->runtime.input.rpa_outdir + aux_filename;
    const std::string legacy_output_filename = this->runtime.input.rpa_outdir + legacy_filename;
    RpaLriDetail::write_librpa_split_basis_file(ucell, type_nw, wfc_l_nchi, wfc_filename);
    RpaLriDetail::write_librpa_split_basis_file(ucell, type_naux, aux_l_nchi, aux_output_filename);

    std::ofstream ofs(legacy_output_filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + legacy_output_filename);
    }
    ofs << std::setw(10) << ucell.ntype << std::setw(10) << total_wfc << std::setw(10) << total_aux << "    fallback"
        << std::endl;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        ofs << std::setw(10) << it + 1 << std::setw(10) << type_nw[static_cast<std::size_t>(it)] << std::setw(10)
            << type_naux[static_cast<std::size_t>(it)] << std::endl;
    }
    const auto write_legacy_shells = [&ofs](const std::vector<std::vector<int>>& shell_counts_by_type) {
        for (std::size_t itype = 0; itype < shell_counts_by_type.size(); ++itype)
        {
            const auto& shell_counts = shell_counts_by_type[itype];
            int nshell = 0;
            for (const int count: shell_counts)
            {
                nshell += count;
            }
            ofs << std::setw(10) << itype + 1 << std::setw(10) << nshell << std::endl;
            for (std::size_t l = 0; l < shell_counts.size(); ++l)
            {
                for (int iradial = 0; iradial < shell_counts[l]; ++iradial)
                {
                    ofs << std::setw(10) << l << std::endl;
                }
            }
        }
    };
    write_legacy_shells(wfc_l_nchi);
    write_legacy_shells(aux_l_nchi);
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_librpa_basis_v1");
}

template void RPA_LRI<double, double>::out_librpa_basis_v1(const UnitCell& ucell,
                                                           Exx_LRI<double>* exx_lri,
                                                           const std::string& aux_filename,
                                                           const std::string& legacy_filename);
template void RPA_LRI<std::complex<double>, double>::out_librpa_basis_v1(const UnitCell& ucell,
                                                                         Exx_LRI<double>* exx_lri,
                                                                         const std::string& aux_filename,
                                                                         const std::string& legacy_filename);

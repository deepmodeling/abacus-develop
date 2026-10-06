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

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k_v1(const UnitCell& ucell,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                         std::string filename,
                                         Exx_LRI<double>* exx_lri)
try
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k_v1");
    ModuleBase::timer::start("RPA_LRI", "out_coulomb_k_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    const auto atom_naux = this->collect_atom_naux_(ucell, exx_lri);
    const int all_mu = RpaLriDetail::sum_int_vector(atom_naux);
    const int nks_tot = this->runtime.input.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    const std::size_t natoms = static_cast<std::size_t>(ucell.nat);

    struct V1Block
    {
        int pair_index = 0;
        int I = 0;
        int J = 0;
        std::int64_t offset = 0;
        RI::Tensor<std::complex<double>> tensor;
    };

    for (int ik = 0; ik != nks_tot; ++ik)
    {
        std::vector<V1Block> blocks;
        for (auto& Ip: Vs)
        {
            const int I = static_cast<int>(Ip.first);
            const int mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                const int J = static_cast<int>(JPp.first.first);
                if (J < I)
                {
                    continue;
                }
                if (!RpaLriDetail::has_valid_matrix_shape(JPp.second))
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const auto R = JPp.first.second;
                const double arg
                    = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec)) * ModuleBase::TWO_PI;
                const std::complex<double> kphase = std::complex<double>(std::cos(arg), std::sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                const int J = static_cast<int>(vq_Jp.first);
                auto& vq_J = vq_Jp.second;
                const int nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[J]);
                if (static_cast<int>(vq_J.shape[0]) != mu_num || static_cast<int>(vq_J.shape[1]) != nu_num)
                {
                    throw std::runtime_error("LibRPA v1 Coulomb output encountered an inconsistent tensor shape.");
                }
                V1Block block;
                block.pair_index = static_cast<int>(RpaLriDetail::coulomb_atom_pair_index(static_cast<std::size_t>(I),
                                                                                          static_cast<std::size_t>(J),
                                                                                          natoms));
                block.I = I;
                block.J = J;
                block.tensor = std::move(vq_J);
                blocks.push_back(std::move(block));
            }
        }

        std::sort(blocks.begin(), blocks.end(), [](const V1Block& lhs, const V1Block& rhs) {
            return lhs.pair_index < rhs.pair_index;
        });
        if (blocks.empty())
        {
            continue;
        }
        for (std::size_t ib = 1; ib < blocks.size(); ++ib)
        {
            if (blocks[ib - 1].pair_index == blocks[ib].pair_index)
            {
                throw std::runtime_error("LibRPA v1 Coulomb output found duplicate atom-pair blocks on one MPI rank.");
            }
        }

        const std::int64_t nblocks = static_cast<std::int64_t>(blocks.size());
        const std::int64_t header_size
            = static_cast<std::int64_t>(6 * sizeof(std::int32_t))
              + static_cast<std::int64_t>(ucell.nat * sizeof(std::int32_t))
              + nblocks * static_cast<std::int64_t>(sizeof(std::int32_t) + sizeof(std::int64_t));
        std::int64_t offset = header_size;
        for (auto& block: blocks)
        {
            block.offset = offset;
            const unsigned long long values = RpaLriDetail::checked_mul_u64(
                static_cast<unsigned long long>(atom_naux[static_cast<std::size_t>(block.I)]),
                static_cast<unsigned long long>(atom_naux[static_cast<std::size_t>(block.J)]),
                "LibRPA v1 Coulomb block size");
            const unsigned long long bytes
                = RpaLriDetail::checked_mul_u64(values,
                                                static_cast<unsigned long long>(sizeof(std::complex<double>)),
                                                "LibRPA v1 Coulomb block size");
            offset += RpaLriDetail::checked_i64_from_u64(bytes, "LibRPA v1 Coulomb block size");
        }

        std::stringstream ss;
        ss << filename << ik + 1 << "_r" << this->runtime.rank << ".dat";
        const std::string out_name = this->runtime.input.rpa_outdir + ss.str();
        std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + out_name);
        }

        const std::int32_t marker = RpaLriDetail::LIBRPA_COULOMB_V1_MARKER;
        const std::int32_t iq = ik + 1;
        const std::int32_t naux = all_mu;
        const std::int32_t value_flag = RpaLriDetail::LIBRPA_COULOMB_V1_COMPLEX_FLAG;
        const std::int32_t natom = ucell.nat;
        const std::int32_t nblock_i32 = static_cast<std::int32_t>(nblocks);
        RpaLriDetail::write_scalar(ofs, marker, out_name);
        RpaLriDetail::write_scalar(ofs, iq, out_name);
        RpaLriDetail::write_scalar(ofs, naux, out_name);
        RpaLriDetail::write_scalar(ofs, value_flag, out_name);
        RpaLriDetail::write_scalar(ofs, natom, out_name);
        RpaLriDetail::write_scalar(ofs, nblock_i32, out_name);
        for (const int atom_aux: atom_naux)
        {
            const std::int32_t atom_aux_i32 = atom_aux;
            RpaLriDetail::write_scalar(ofs, atom_aux_i32, out_name);
        }
        for (const auto& block: blocks)
        {
            const std::int32_t pair_index = block.pair_index;
            RpaLriDetail::write_scalar(ofs, pair_index, out_name);
            RpaLriDetail::write_scalar(ofs, block.offset, out_name);
        }
        for (const auto& block: blocks)
        {
            const std::size_t nvalues = block.tensor.get_shape_all();
            if (nvalues == 0)
            {
                throw std::runtime_error("LibRPA v1 Coulomb output encountered an empty tensor payload.");
            }
            RpaLriDetail::checked_write(ofs, block.tensor.ptr(), nvalues * sizeof(std::complex<double>), out_name);
        }
        ofs.close();
    }

    ModuleBase::timer::end("RPA_LRI", "out_coulomb_k_v1");
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_coulomb_k_v1");
}

template void RPA_LRI<double, double>::out_coulomb_k_v1(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Vs,
    std::string filename,
    Exx_LRI<double>* exx_lri);
template void RPA_LRI<std::complex<double>, double>::out_coulomb_k_v1(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Vs,
    std::string filename,
    Exx_LRI<double>* exx_lri);

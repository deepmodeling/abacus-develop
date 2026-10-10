#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_UTILS_LR_UTIL_HCONTAINER_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_UTILS_LR_UTIL_HCONTAINER_H

#include "source_estate/module_dm/density_matrix.h"
#include <numeric>
#include  "source_base/parallel_reduce.h"
#include "source_base/macros.h"
#include <ATen/core/tensor.h>
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_hs/lat_r_csr.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#ifdef __EXX
#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "source_lcao/module_ri/ri_2d_comm.h"
#endif
namespace LR_Util
{
    template <typename T>
    using Real = typename GetTypeReal<T>::type;
    template<typename TR>
    void print_HR(const hamilt::HContainer<TR>& HR, const std::string& label, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        for (int iap = 0; iap < HR.size_atom_pairs(); ++iap)
        {
            auto ap = HR.get_atom_pair(iap);
            const int ia = ap.get_atom_i();
            const int ja = ap.get_atom_j();
            for (int iR = 0;iR < ap.get_R_size();++iR)
            {
                std::cout << "atom pair (" << ia << ", " << ja << "),  "
                    << "R=(" << ap.get_R_index(iR)[0] << ", " << ap.get_R_index(iR)[1] << ", " << ap.get_R_index(iR)[2] << "): \n";
                auto& mat = ap.get_HR_values(iR);
                std::cout << "rowsize=" << ap.get_row_size() << ", colsize=" << ap.get_col_size() << "\n";
                for (int i = 0;i < ap.get_row_size();++i)
                {
                    for (int j = 0;j < ap.get_col_size();++j)
                    {
                        auto& v = mat.get_value(i, j);
                        std::cout << (std::abs(v) > threshold ? v : 0) << " ";
                    }
                    std::cout << "\n";
                }
            }
        }
    }

    template <typename TK, typename TR>
    void print_DMR(const module_dm::DensityMatrix<TK, TR>& DMR, const std::string& label, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        int is = 0;
        for (auto& dr : DMR.get_dmr_vec())
            print_HR(*dr, "DMR[ispin=s" + std::to_string(is++) + "]", threshold);
    }

    /// copy one atom pair's HContainer values from a complex DMR into a real-valued one, by
    /// matching atom-pair/R rather than by nat-based index iteration (so it needs no `nat`).
    template<typename TR_from, typename TR_to>
    inline void copy_hcontainer_real_imag_part(const hamilt::HContainer<TR_from>& HR,
        hamilt::HContainer<TR_to>& HR_real, const bool get_imag)
    {
        for (int iap = 0; iap < HR.size_atom_pairs(); ++iap)
        {
            auto ap = &const_cast<hamilt::HContainer<TR_from>&>(HR).get_atom_pair(iap);
            const int ia = ap->get_atom_i();
            const int ja = ap->get_atom_j();
            auto ap_real = HR_real.find_pair(ia, ja);
            // under MPI-parallel (2D block-cyclic) HContainer, an atom pair not owned by this
            // rank is absent from find_pair() and returns nullptr here; skip it
            if (!ap_real) { continue; }
            for (int iR = 0;iR < ap->get_R_size();++iR)
            {
                // R index may be different between the two HContainers, find by R value instead of R-index
                auto dR = ap->get_R_index(iR);
                auto ptr = ap->get_HR_values(iR).get_pointer();
                auto ptr_real = ap_real->get_HR_values(dR.x, dR.y, dR.z).get_pointer();
                for (int i = 0;i < ap->get_size();++i) { ptr_real[i] = (get_imag ? std::imag(ptr[i]) : std::real(ptr[i])); }
            }
        }
    }

    template<typename T>
    void get_DMR_real_imag_part(const module_dm::DensityMatrix<T, T>& DMR,
        module_dm::DensityMatrix<T, Real<T>>& DMR_real,
        const char& type = 'R')
    {
        assert(DMR.get_dmr_vec().size() == DMR_real.get_dmr_vec().size());
        const bool get_imag = (type == 'I' || type == 'i');
        for (size_t is = 0;is < DMR.get_dmr_vec().size();++is)
        {
            auto dr = DMR.get_dmr_vec()[is]; //get_dmr_ptr() has bug when is=0
            auto dr_real = DMR_real.get_dmr_vec()[is];
            assert(dr != nullptr);
            assert(dr_real != nullptr);
            copy_hcontainer_real_imag_part(*dr, *dr_real, get_imag);
        }
    }

    /// overload: only copy the `is`-th spin channel of DMR (source) into the (single-channel)
    /// DMR_real, to avoid mixing/overlapping spin channels when DMR has more than one spin channel
    template<typename T>
    void get_DMR_real_imag_part(const module_dm::DensityMatrix<T, T>& DMR,
        module_dm::DensityMatrix<T, Real<T>>& DMR_real,
        const int& is,
        const char& type = 'R')
    {
        assert(static_cast<size_t>(is) < DMR.get_dmr_vec().size());
        assert(DMR_real.get_dmr_vec().size() == 1);
        const bool get_imag = (type == 'I' || type == 'i');
        auto dr = DMR.get_dmr_vec()[is]; //get_dmr_ptr() has bug when is=0
        auto dr_real = DMR_real.get_dmr_vec()[0];
        assert(dr != nullptr);
        assert(dr_real != nullptr);
        copy_hcontainer_real_imag_part(*dr, *dr_real, get_imag);
    }

    inline void set_HR_real_imag_part(const hamilt::HContainer<double>& HR_real,
        hamilt::HContainer<std::complex<double>>& HR,
        const char& type)
    {
        bool get_imag = (type == 'I' || type == 'i');
        for (int iap = 0; iap < HR.size_atom_pairs(); ++iap)
        {
            auto ap = &HR.get_atom_pair(iap);
            const int ia = ap->get_atom_i();
            const int ja = ap->get_atom_j();
            auto ap_real = HR_real.find_pair(ia, ja);
            assert(ap_real != nullptr);
            for (int iR = 0;iR < ap->get_R_size();++iR)
            {
                // R index may be different between the two HContainers, find by R value instead of R-index
                auto dR = ap->get_R_index(iR);
                auto ptr = ap->get_HR_values(iR).get_pointer();
                auto ptr_real = ap_real->get_HR_values(dR.x, dR.y, dR.z).get_pointer();
                for (int i = 0;i < ap->get_size();++i) { get_imag ? ptr[i].imag(ptr_real[i]) : ptr[i].real(ptr_real[i]); }
            }
        }
    }

    template <typename T, typename TR>
    void initialize_HR(hamilt::HContainer<TR>& hR,
                       const UnitCell& ucell,
                       const Grid_Driver& gd,
                       const std::vector<double>& orb_cutoff)
    {
        const auto& pmat = *hR.get_paraV();
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            auto tau1 = ucell.get_tau(iat1);
            int T1, I1;
            ucell.iat2iait(iat1, &I1, &T1);
            AdjacentAtomInfo adjs;
            gd.Find_atom(ucell, tau1, T1, I1, &adjs);
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                int iat2 = ucell.itia2iat(T2, I2);
                if (pmat.is_invalid_atom_pair(iat1, iat2)) { continue; }
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
                if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0 >= orb_cutoff[T1] + orb_cutoff[T2]) { continue; }
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, &pmat);
                hR.insert_pair(tmp);
            }
        }
        hR.allocate(nullptr, true);
        // hR.set_paraV(&pmat);
        if (std::is_same<T, double>::value) { hR.fix_gamma(); }
    }
    template <typename T, typename TR>
    void initialize_DMR(module_dm::DensityMatrix<T, TR>& dm,
                        const Parallel_Orbitals& pmat,
                        const UnitCell& ucell,
                        const Grid_Driver& gd,
                        const std::vector<double>& orb_cutoff)
    {
        hamilt::HContainer<TR> hR_tmp(&pmat);
        initialize_HR<T, TR>(hR_tmp, ucell, gd, orb_cutoff);
        dm.init_dmr(hR_tmp);
    }

    /// $\sum_{uvR} H1_{uv}(R) H2_{uv}(R)$
    template<typename TR1, typename TR2>
    TR1 dot_R_matrix(const hamilt::HContainer<TR1>& h1, const hamilt::HContainer<TR2>& h2)
    {
        TR1 sum = 0;
        // Match by atom/R values: containers may store pairs and images in different orders.
        for (int iap = 0; iap < h1.size_atom_pairs(); ++iap)
        {
            const auto& ap1 = h1.get_atom_pair(iap);
            const int iat1 = ap1.get_atom_i();
            const int iat2 = ap1.get_atom_j();
            const auto* ap2 = h2.find_pair(iat1, iat2);
            assert(ap2);
            for (int iR = 0; iR < ap1.get_R_size(); ++iR)
            {
                const auto& R = ap1.get_R_index(iR);
                const auto& mat1 = ap1.get_HR_values(R.x, R.y, R.z);
                const auto& mat2 = ap2->get_HR_values(R.x, R.y, R.z);
                const int size = mat1.get_col_size() * mat1.get_row_size();
                const TR1* begin = mat1.get_pointer();
                const TR1* end = begin + size;
                const TR2* values = mat2.get_pointer();
                const TR1 zero = 0;
                sum += std::inner_product(begin, end, values, zero);
            }
        }
        // Parallel_Reduce::reduce_all(sum);  // not needed, since it will be reduced outside
        return sum;
    }


    // cal_dmr converts column-major DMK into row-major DMR without changing AO indices.
    // It already uses the ground-state DMR convention, including exp(+ik.R) for multi-k.
    // A physical transpose must therefore be applied to DMK, never by swapping atom pairs.
    template<typename TK>
    void transpose_DMR(module_dm::DensityMatrix<TK, double>& dm, const Parallel_Orbitals& pv)
    {
        // 1. transpose dm(k)
        for (auto& dk : dm.get_dmk_vec())
        {
#ifdef __MPI
            // dm(k) is 2D-block-cyclic distributed, so the transpose needs the PBLAS routine
            // `mattrans` (pdtran_/pztranc_) rather than a plain serial swap.
            LR_Util::mattrans(dk.data(), pv.get_global_row_size(), pv);
#else
            // In a serial build the entire square AO matrix is local, in column-major order.
            const int n = pv.get_global_row_size();
            std::vector<TK> transposed(dk.size());
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < n; ++i)
                {
                    transposed[j * n + i] = LR_Util::get_conj(dk[i * n + j]);
                }
            }
            dk.swap(transposed);
#endif
        }

        // 2. FT
        dm.cal_dmr(-1);
    }
    template<typename TK>
    void transpose_DMR(module_dm::DensityMatrix<TK, std::complex<double>>& dm, const Parallel_Orbitals& pv)
    {
        throw std::runtime_error("transpose_DMR is not implemented for complex DMR, due to the lack of minus-sign FT.");
        // 1. dm(k) dagger
        for (auto& dk : dm.get_dmk_vec())
        {
#ifdef __MPI
            LR_Util::mattrans(dk.data(), pv.get_global_row_size(), pv);
#else
            throw std::runtime_error("transpose_DMR requires MPI (PBLAS mattrans) for the 2D-block-cyclic dm(k) transpose.");
#endif
        }

        // 2. FT with the minus sign in the exponent (TO DO)
        dm.cal_dmr(-1);
    }

    template<typename TK, typename TR>
    module_dm::DensityMatrix<TK, TR> build_dm_from_dmk(const std::vector<ct::Tensor>& dmk,
        const Parallel_Orbitals& pmat,
        const int& nk,
        const std::vector<ModuleBase::Vector3<double>>& kvec_d,
        const UnitCell& ucell,
        const Grid_Driver& gd,
        const std::vector<double>& orb_cutoff,
        const bool symmetrize = false,
        const bool cal_dmr = true,
        const bool transpose = false)
    {
        module_dm::DensityMatrix<TK, TR> dm(&pmat, 1, kvec_d, nk);
        initialize_DMR(dm, pmat, ucell, gd, orb_cutoff);

        if (symmetrize)
            for (int ik = 0; ik < nk; ++ik)
#ifdef __MPI
                LR_Util::matsym(dmk[ik].data<TK>(), pmat.get_global_row_size(), pmat);
#else
                LR_Util::matsym(dmk[ik].data<TK>(), pmat.get_global_row_size());
#endif

        for (int ik = 0; ik < nk; ++ik)
            dm.set_dmk_ptr(ik, dmk[ik].data<TK>());

        if (cal_dmr)
        {
            dm.cal_dmr(-1);
        }
        return dm;
    }

    /// @brief Spin-resolved counterpart of `build_dm_from_dmk`: one `dmk` list per spin channel.
    ///
    /// Needed by the open-shell gradient, where $D^X$, $T$, $D^Z$ and the energy-weighted
    /// density matrix all have two independent channels. `DensityMatrix` stores DMK as a flat
    /// `[nspin][nk]` array, so channel `is` starts at `is * nk`.
    template <typename TK, typename TR>
    module_dm::DensityMatrix<TK, TR> build_dm_from_dmk_spin(const std::vector<std::vector<ct::Tensor>>& dmk,
        const Parallel_Orbitals& pmat,
        const int& nk,
        const std::vector<ModuleBase::Vector3<double>>& kvec_d,
        const UnitCell& ucell,
        const Grid_Driver& gd,
        const std::vector<double>& orb_cutoff,
        const bool symmetrize = false,
        const bool cal_dmr = true)
    {
        const int nspin_dm = static_cast<int>(dmk.size());
        module_dm::DensityMatrix<TK, TR> dm(&pmat, nspin_dm, kvec_d, nk);
        initialize_DMR(dm, pmat, ucell, gd, orb_cutoff);
        for (int is = 0; is < nspin_dm; ++is)
        {
            assert(static_cast<int>(dmk[is].size()) >= nk);
            if (symmetrize)
            {
                for (int ik = 0; ik < nk; ++ik)
                {
#ifdef __MPI
                    LR_Util::matsym(dmk[is][ik].data<TK>(), pmat.get_global_row_size(), pmat);
#else
                    LR_Util::matsym(dmk[is][ik].data<TK>(), pmat.get_global_row_size());
#endif
                }
            }
            for (int ik = 0; ik < nk; ++ik) { dm.set_dmk_ptr(is * nk + ik, dmk[is][ik].data<TK>()); }
        }
        if (cal_dmr)
        {
            dm.cal_dmr(-1);
        }
        return dm;
    }

    namespace sparse_format
    {
        // ref: sparse_format::cal_HContainer_d/cd and sparse_format::cal_HSR
        // but more general(not depend on LCAO_HS_Arrays)
        template<typename TR>
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, TR>>>
            get_sparse_format(
                const hamilt::HContainer<TR>& hR,
                const Parallel_Orbitals& pv,
                const double& sparse_thr = 1e-10)
        {
            std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, TR>>> target;
            auto row_indexes = pv.get_indexes_row();
            auto col_indexes = pv.get_indexes_col();
            for (int iap = 0; iap < hR.size_atom_pairs(); ++iap) {
                int atom_i = hR.get_atom_pair(iap).get_atom_i();
                int atom_j = hR.get_atom_pair(iap).get_atom_j();
                int start_i = pv.atom_begin_row[atom_i];
                int start_j = pv.atom_begin_col[atom_j];
                int row_size = pv.get_nrow_atom(atom_i);
                int col_size = pv.get_ncol_atom(atom_j);
                for (int iR = 0; iR < hR.get_atom_pair(iap).get_R_size(); ++iR) {
                    auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
                    const ModuleBase::Vector3<int> r_index
                        = hR.get_atom_pair(iap).get_R_index(iR);
                    Abfs::Vector3_Order<int> dR(r_index.x, r_index.y, r_index.z);
                    for (int i = 0; i < row_size; ++i) {
                        int mu = row_indexes[start_i + i];
                        for (int j = 0; j < col_size; ++j) {
                            int nu = col_indexes[start_j + j];
                            const auto& value_tmp = matrix.get_value(i, j);
                            if (std::abs(value_tmp) > sparse_thr) {
                                target[dR][mu][nu] = value_tmp;
                            }
                        }
                    }
                }
            }
            return target;
        }

        //a more general version of save_HSR_sparse (not depend on LCAO_HS_Arrays)
        template<typename TR>
        void save_sparse(
            const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, TR>>>& smat,
            // const std::set<Abfs::Vector3_Order<int>>& all_R_coor,
            // const bool& binary,
            const std::string& filename,
            const Parallel_Orbitals& pv,
            const std::string& out_dir,
            const int& nlocal,
            const int& my_rank,
            const double& sparse_thr = 1e-10)
        {
            // calculate the total number of non-zero elements of the (nbasis, nbasis) matrix for each R
            std::vector<int> non_zero_counts(smat.size(), 0);   // number of Rs
            int i = 0;
            for (const auto& Rij : smat)
                non_zero_counts[i++] = std::accumulate(Rij.second.begin(), Rij.second.end(), 0,
                    [](int sum, const std::pair<const size_t, std::map<size_t, TR>>& line) { return sum + line.second.size(); });
            Parallel_Reduce::reduce_all(non_zero_counts.data(), non_zero_counts.size());


            std::string out_file = out_dir + filename;
            std::ofstream ofs;
            if (my_rank == 0)
            {
                ofs.open(out_file);
                // if (binary) ofs.open(out_file, std::ios::binary);
                ofs << "STEP: 0" << std::endl;
                ofs << "Matrix Dimension: " << nlocal << std::endl;
                ofs << "Matrix number: " << non_zero_counts.size() << std::endl;
            }
            i = 0;
            for (const auto& Rij : smat)
            {
                const auto& R = Rij.first;
                ofs << R.x << " " << R.y << " " << R.z << " " << non_zero_counts[i++] << std::endl;
                ModuleIO::SparseWriteOptions single_R_options;
                single_R_options.threshold = sparse_thr;
                single_R_options.binary = false;
                ModuleIO::save_lat_r(ofs, Rij.second, pv, single_R_options);
            }
            if (my_rank == 0) { ofs.close(); }
        }
    }

    template<typename TR>
    void save_HR(
        const hamilt::HContainer<TR>& hR,
        // const std::set<Abfs::Vector3_Order<int>>& all_R_coor,
        // const bool& binary,
        const std::string& filename,
        const Parallel_Orbitals& pv,
        const std::string& out_dir,
        const int& nlocal,
        const int& my_rank,
        const double& sparse_thr = 1e-10)
    {
        sparse_format::save_sparse(sparse_format::get_sparse_format(hR, pv, sparse_thr),
            filename, pv, out_dir, nlocal, my_rank, sparse_thr);
    }

    template <typename TK, typename TR>
    void save_DMR(const module_dm::DensityMatrix<TK, TR>& DMR,
        const std::string& filename,
        const Parallel_Orbitals& pv,
        const std::string& out_dir,
        const int& nlocal,
        const int& my_rank,
        const double& sparse_thr = 1e-10)
    {
        int is = 0;
        for (auto& dr : DMR.get_dmr_vec())
            save_HR(*dr, filename + "_s" + std::to_string(is), pv, out_dir, nlocal, my_rank, sparse_thr);
    }

#ifdef __EXX
    // convert DensityMatrix to maps of RI::Tensors
    // return 0.5*D[0]
    template <typename TK, typename TR>
    auto get_exx_Ds_spin1(const module_dm::DensityMatrix<TK, TR>& dm,
        const UnitCell& ucell, const K_Vectors& kv, const Parallel_Orbitals& pmat)
        -> std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<TK>>>
    {
        const int& nk = dm.get_DMK_nks();   // nks/nspin
        std::vector<const std::vector<TK>*> DMk_trans_pointer(nk);
        for (int ik = 0;ik < nk;++ik) { DMk_trans_pointer[ik] = &dm.get_dmk_vec()[ik]; }
        return RI_2D_Comm::split_m2D_ktoR<TK>(ucell, kv, DMk_trans_pointer, pmat, /*nspin=*/1)[0];
    }
    // return SPIN_multiple*D[0] as implemented in split_m2D_ktoR
    // SPIN_multiple = map({ {1,0.5}, {2,1}, {4,1} }).at(nspin)
    template <typename TK, typename TR>
    auto get_exx_Ds_gs(const module_dm::DensityMatrix<TK, TR>& dm,
        const UnitCell& ucell, const K_Vectors& kv, const Parallel_Orbitals& pmat)
        -> std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<TK>>>>
    {
        const int& nspin = dm.get_dmr_vec().size();
        const int nks = dm.get_DMK_nks();
        std::vector<const std::vector<TK>*> DMk_trans_pointer(nks);
        for (int iks = 0;iks < dm.get_DMK_nks();++iks)
            DMk_trans_pointer[iks] = &dm.get_dmk_vec()[iks];
        return RI_2D_Comm::split_m2D_ktoR<TK>(ucell, kv, DMk_trans_pointer, pmat, nspin);
    }
#endif
}

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_UTILS_LR_UTIL_HCONTAINER_H

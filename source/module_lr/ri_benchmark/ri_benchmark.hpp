#pragma once
#include "ri_benchmark.h"
#include "module_base/module_container/base/third_party/blas.h"
namespace RI_Benchmark
{
    // std::cout << "the size of Cs:" << std::endl;
    // for (auto& it1: Cs)
    // {
    //     for (auto& it2: it1.second)
    //     {
    //         std::cout << "iat1=" << it1.first << ", iat2= " << it2.first.first << std::endl;
    //         auto& ts = it2.second.shape;
    //         std::cout << "Tensor shape:" << ts[0] << " " << ts[1] << " " << ts[2] << std::endl; // abf, nw1, nw2
    //     }
    // }

    template <typename TR>
    inline int count_nao_from_Cs(const TLRI<TR>& Cs_ao)
    {
        int naos = 0;
        for (const auto& it2 : Cs_ao.at(0))
        {
            assert(it2.second.shape.size() == 3);
            naos += it2.second.shape[2];
        }
        return naos;
    }

    /// @brief slice the psi from wfc_ks of all the k-points, some bands, some basis functions
    /// @return
    template <typename TK>
    inline std::vector<TK> slice_psi(const psi::Psi<TK>& wfc_ks,
        const int& ibstart,
        const int& nb,
        const int& iwstart,
        const int& nw)
    {
        std::vector<TK> psi_slice(wfc_ks.get_nk() * nb * nw);
        int i = 0;
        for (int ik = 0; ik < wfc_ks.get_nk(); ++ik)
            for (int ib = ibstart; ib < ibstart + nb; ++ib)
                for (int iw = iwstart; iw < iwstart + nw; ++iw)
                    psi_slice[i++] = wfc_ks(ik, ib, iw);
        assert(i == psi_slice.size());
        return psi_slice;
    }

    template <typename TK, typename TR>
    TLRI<TK> cal_Cs_mo(const UnitCell& ucell,
        const TLRI<TR>& Cs_ao,
        const psi::Psi<TK>& wfc_ks,
        const int& nocc,
        const int& nvirt,
        const int& occ_first,
        const bool& read_from_aims,
        const std::vector<int>& aims_nbasis)
    {
        // assert(wfc_ks.get_nk() == 1);   // currently only gamma-only is supported
        assert(nocc + nvirt <= wfc_ks.get_nbands());
        const bool use_aims_nbasis = (read_from_aims && !aims_nbasis.empty());
        TLRI<TK> Cs_mo;
        int iw1 = 0;
        for (auto& c1 : Cs_ao)
        {
            const int& iat1 = c1.first;
            const int& it1 = ucell.iat2it[iat1];
            const int& nw1 = (use_aims_nbasis ? aims_nbasis[it1] : ucell.atoms[it1].nw);
            if (!use_aims_nbasis) { assert(iw1 == ucell.get_iat2iwt()[iat1]); }
            int iw2 = 0;
            for (auto& c2 : c1.second)
            {
                const int& iat2 = c2.first.first;
                const int& it2 = ucell.iat2it[iat2];
                const int& nw2 = (use_aims_nbasis ? aims_nbasis[it2] : ucell.atoms[it2].nw);
                if (!use_aims_nbasis) { assert(iw2 == ucell.get_iat2iwt()[iat2]); }

                const auto& tensor_ao = c2.second;
                const size_t& nabf = tensor_ao.shape[0];
                assert(tensor_ao.shape.size() == 3); // abf, nw1, nw2

                std::vector<TK> psi_a1 = occ_first ? slice_psi(wfc_ks, 0, nocc, iw1, nw1)   //occ
                    : slice_psi(wfc_ks, nocc, nvirt, iw1, nw1); // virt
                std::vector<TK> psi_a2 = occ_first ? slice_psi(wfc_ks, nocc, nvirt, iw2, nw2)   //virt
                    : slice_psi(wfc_ks, 0, nocc, iw2, nw2); // occ

                Cs_mo[c1.first][c2.first] = RI::Tensor<TK>({ nabf, (std::size_t)nocc, (std::size_t)nvirt });
                for (int iabf = 0; iabf < nabf; ++iabf)
                {
                    const auto ptr = &tensor_ao(iabf, 0, 0);
                    std::vector<TK> tmp(nw1 * (occ_first ? nvirt : nocc));
                    // caution: Cs are row-major  (ia2 contiguous)
                    if (occ_first)
                    {
                        container::BlasConnector::gemm('T', 'N', nvirt, nw1, nw2, 1.0, psi_a2.data(), nw2, ptr, nw2, 0.0, tmp.data(), nvirt);
                        container::BlasConnector::gemm('N', 'N', nvirt, nocc, nw1, 1.0, tmp.data(), nvirt, psi_a1.data(), nw1, 0.0, &Cs_mo[c1.first][c2.first](iabf, 0, 0), nvirt);
                    }
                    else
                    {
                        container::BlasConnector::gemm('T', 'N', nw1, nocc, nw2, 1.0, ptr, nw2, psi_a2.data(), nw2, 0.0, tmp.data(), nw1);
                        container::BlasConnector::gemm('T', 'N', nvirt, nocc, nw1, 1.0, psi_a1.data(), nw1, tmp.data(), nw1, 0.0, &Cs_mo[c1.first][c2.first](iabf, 0, 0), nvirt);
                    }
                }
                iw2 += nw2;
            }
            iw1 += nw1;
        }
        return Cs_mo;
    }

    template <typename TK, typename TR>
    std::vector<TK> cal_Amat_full(const TLRI<TK>& Cs_a,
        const TLRI<TK>& Cs_b,
        const TLRI<TR>& Vs)
    {
        assert(Cs_a.size() > 0);
        assert(Cs_b.size() > 0);
        assert(Vs.size() > 0);
        auto& Cs_shape = Cs_a.at(0).begin()->second.shape;
        auto& Vs_shape = Vs.at(0).begin()->second.shape;
        assert(Cs_shape.size() == 3); // abf, nocc, nvirt
        assert(Cs_shape.size() == 3); // abf, nocc, nvirt
        assert(Vs_shape.size() == 2); // abf, abf

        const int& npairs = Cs_shape[1] * Cs_shape[2];
        std::vector<TK> Amat_full(npairs * npairs, 0.0);

        for (auto& itv1 : Vs)
        {
            const int& iat1 = itv1.first;
            for (auto& itv2 : itv1.second)
            {
                const int& iat2 = itv2.first.first;
                const auto& tensor_v = itv2.second; // (nabf2, nabf1), T
                for (auto& itca2 : Cs_a.at(iat1))
                {
                    const int& iat3 = itca2.first.first;
                    const auto& tensor_ca = itca2.second; // (nvirt*nocc, nabf1), N
                    for (auto& itcb2 : Cs_b.at(iat2))
                    {
                        const int& iat4 = itcb2.first.first;
                        const auto& tensor_cb = itcb2.second; // (nvirt*nocc, nabf2), T
                        const int& nabf1 = tensor_v.shape[0];
                        const int& nabf2 = tensor_v.shape[1];
                        assert(tensor_ca.shape[0] == nabf1);   //abf1
                        assert(tensor_cb.shape[0] == nabf2); //abf2
                        std::vector<TK> tmp(npairs * nabf1);
                        container::BlasConnector::gemm('T', 'T', nabf1, npairs, nabf2, 1.0, tensor_v.ptr(), nabf2, tensor_cb.ptr(), npairs, 0.0, tmp.data(), nabf1);
                        container::BlasConnector::gemm('N', 'N', npairs, npairs, nabf1, 2.0/*Hartree to Ry*/, tensor_ca.ptr(), npairs, tmp.data(), nabf1, 1.0, Amat_full.data(), npairs);
                    }
                }
            }
        }
        return Amat_full;
    }
    template <typename TK>
    TLRIX<TK> cal_CsX(const TLRI<TK>& Cs_mo, const TK* X)
    {
        TLRIX<TK> CsX;
        for (auto& it1 : Cs_mo)
        {
            const int& iat1 = it1.first;
            for (auto& it2 : it1.second)
            {
                const int& iat2 = it2.first.first;
                auto& tensor_c = it2.second;
                const int& nabf = tensor_c.shape[0];
                const int& npairs = tensor_c.shape[1] * tensor_c.shape[2];
                std::vector<TK> CX(nabf);
                for (int iabf = 0;iabf < nabf;++iabf)
                    CX[iabf] = container::BlasConnector::dot(npairs, &tensor_c(iabf, 0, 0), 1, X, 1);
                CsX[iat1][it2.first] = CX;
            }
        }
        return CsX;
    }

    template <typename TK, typename TR>
    TLRI<TK> cal_CV(const TLRI<TK>& Cs_a_mo,
        const TLRI<TR>& Vs)
    {
        TLRI<TK> CV;
        for (auto& it1 : Cs_a_mo)   //the atom on which ABFs locate
        {
            const int& iat1 = it1.first;
            for (auto& it2 : it1.second)
            {
                const int& iat2 = it2.first.first;
                // const auto& Rc=it2.first.second;
                auto& tensor_c = it2.second;
                const int& nabf1 = tensor_c.shape[0];
                const int& npairs = tensor_c.shape[1] * tensor_c.shape[2];
                for (auto& it3 : Vs.at(iat1))
                {
                    const int& iat3 = it3.first.first;
                    // const auto& Rv=it3.first.second;
                    const auto& tensor_v = it3.second;
                    assert(nabf1 == tensor_v.shape[0]);
                    const size_t& nabf2 = tensor_v.shape[1];
                    std::vector<TK> tmp(nabf2 * npairs);
                    // const auto& Rcv = (Rv - Rc) % period;
                    if (CV.count(iat2) && CV.at(iat2).count({ iat3, {0, 0, 0} }))   // add-up, sum over iat1
                    {
                        auto& tensor_cv = CV.at(iat2).at({ iat3, {0, 0, 0} });
                        container::BlasConnector::gemm('N', 'T', npairs, nabf2, nabf1, 1.0, tensor_c.ptr(), npairs, tensor_v.ptr(), nabf2, 1.0, tensor_cv.ptr(), npairs);
                    }
                    else
                    {
                        RI::Tensor<TK> tmp({ nabf2, tensor_c.shape[1], tensor_c.shape[2] });  // (nabf2, nocc, nvirt)
                        container::BlasConnector::gemm('N', 'T', npairs, nabf2, nabf1, 1.0, tensor_c.ptr(), npairs, tensor_v.ptr(), nabf2, 0.0, tmp.ptr(), npairs);
                        CV[iat2][{iat3, { 0, 0, 0 }}] = tmp;
                    }
                }
            }
        }
        return CV;
    }
    template <typename TK, typename TR>
    void cal_AX(const TLRI<TK>& Cs_a,
        const TLRIX<TK>& Cs_bX,
        const TLRI<TR>& Vs,
        TK* AX,
        const double& scale)
    {
        const int& npairs = Cs_a.at(0).begin()->second.shape[1] * Cs_a.at(0).begin()->second.shape[2];
        for (auto& itv1 : Vs)
        {
            const int& iat1 = itv1.first;
            for (auto& itv2 : itv1.second)
            {
                const int& iat2 = itv2.first.first;
                const auto& tensor_v = itv2.second; // (nabf2, nabf1), T
                for (auto& itca2 : Cs_a.at(iat1))
                {
                    const int& iat3 = itca2.first.first;
                    const auto& tensor_ca = itca2.second; // (nvirt*nocc, nabf1), N
                    for (auto& itcb2 : Cs_bX.at(iat2))
                    {
                        const int& iat4 = itcb2.first.first;
                        const auto& vector_cb = itcb2.second; // (nvirt*nocc, nabf2), T
                        const int& nabf1 = tensor_v.shape[0];
                        const int& nabf2 = tensor_v.shape[1];
                        assert(tensor_ca.shape[0] == nabf1);   //abf1
                        assert(vector_cb.size() == nabf2); //abf2
                        std::vector<TK> tmp(nabf1);
                        container::BlasConnector::gemv('T', nabf1, nabf2, 1.0, tensor_v.ptr(), nabf2, vector_cb.data(), 1, 0.0, tmp.data(), 1);
                        container::BlasConnector::gemv('N', npairs, nabf1, scale/*Hartree to Ry; singlet*/, tensor_ca.ptr(), npairs, tmp.data(), 1, 1.0, AX, 1);
                    }
                }
            }
        }
    }

    template <typename TK>
    void cal_AX(const TLRI<TK>& CV,
        const TLRIX<TK>& Cs_bX,
        TK* AX,
        const double& scale)
    {
        for (auto& it1 : CV)
        {
            const int& iat1 = it1.first;
            for (auto& it2 : it1.second)
            {
                const int& iat2 = it2.first.first;
                const auto& tensor_cv = it2.second; // (nabf, nocc, nvirt)
                const int& npairs = tensor_cv.shape[1] * tensor_cv.shape[2];
                for (auto& it3 : Cs_bX.at(iat2))
                {
                    const int& iat3 = it3.first.first;
                    const auto& vector_cx = it3.second; // (nabf)
                    const int& nabf = tensor_cv.shape[0];
                    assert(vector_cx.size() == nabf); //abf on at2
                    container::BlasConnector::gemv('N', npairs, nabf, scale/*Hartree to Ry; singlet*/, tensor_cv.ptr(), npairs, vector_cx.data(), 1, 1.0, AX, 1);
                }
            }
        }
    }

    template <typename FPTYPE>
    std::vector<FPTYPE> read_aims_ebands(const std::string& file, const int nocc, const int nvirt, int& ncore)
    {
        std::vector<FPTYPE> bands;
        std::vector<FPTYPE> bands_final;
        std::ifstream ifs;
        ifs.open(file);
        std::string tmp;
        FPTYPE ene, occ;
        for (int i = 0;i < 6;++i) { std::getline(ifs, tmp); } // skip the first 6 lines
        int ivirt = 0;
        while (ifs.peek() != EOF) {
            ifs >> tmp >> occ >> ene >> tmp;
            std::cout << "occ=" << occ << ", ene=" << ene << std::endl;
            bands.push_back(ene * 2);//Hartree to Ry
            if (occ < 0.1) { ++ivirt; }
            if (ivirt == nvirt) { break; }
        }
        ncore = bands.size() - nocc - nvirt;
        std::cout << "bands_final:" << std::endl;
        for (int i = ncore;i < bands.size();++i)
        {
            bands_final.push_back(bands[i]);
            std::cout << bands[i] << "  ";
        }
        std::cout << std::endl;
        return bands_final;
    }
    template <typename TK>
    void read_aims_eigenvectors(psi::Psi<TK>& wfc_ks, const std::string& file, const int ncore, const int nbands, const int nbasis)
    {
        std::ifstream ifs;
        ifs.open(file);
        std::string tmp;
        int nbands_last = 0;
        while (ifs.peek() != EOF)
        {
            std::getline(ifs, tmp); //the first line
            std::stringstream ss(tmp);
            while (std::getline(ss, tmp, ' ')) {};
            int nbands_file = std::stoi(tmp);
            for (int iw = 0;iw < nbasis;++iw)
            {
                ifs >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;  //useless cols
                for (int ib = nbands_last; ib < nbands_file;++ib)
                {
                    ifs >> tmp;
                    if (ib >= ncore && ib < ncore + nbands)
                    {
                        for (int is = 0;is < wfc_ks.get_nk();++is)
                        {   //only for gamma_only and spin degenerate
                            wfc_ks(is, ib - ncore, iw) = std::stod(tmp);
                        }
                    }
                }
            }
            std::getline(ifs, tmp); // the interval line between two blocks
            std::getline(ifs, tmp); // the interval line between two blocks
            nbands_last = nbands_file;
        }
        // output wfc
        std::cout << "wfc_gs_read_from_aims:" << std::endl;
        for (int ib = 0;ib < nbands;++ib)
        {
            for (int iw = 0;iw < nbasis;++iw)
            {
                std::cout << wfc_ks(0, ib, iw) << "  ";
            }
            std::cout << std::endl;
        }
    }
    template < typename TR> // only for blocking by atom pairs
    TLRI<TR> read_coulomb_mat(const std::string& file, const TLRI<TR>& Cs)
    {   //for gamma_only, V(q)=V(R=0)
        std::ifstream ifs;
        ifs.open(file);
        size_t nks = 0, nabf = 0, istart = 0, jstart = 0, iend = 0, jend = 0;
        std::string tmp;
        ifs >> nks;//   nkstot=1
        if (nks > 1) { std::cout << "Warning: nks>1 is not supported yet!" << std::endl; }
        TLRI<TR> Vs;
        const int nat = Cs.size();
        for (int iat1 = 0;iat1 < nat;++iat1)
        {
            const size_t nabf1 = Cs.at(iat1).at({ 0, {0,0,0} }).shape[0];
            for (int iat2 = 0;iat2 < nat;++iat2)
            {
                if (iat1 > iat2)
                {   // coulomb_mat has only the upper triangle part
                    Vs[iat1][{iat2, { 0,0,0 }}] = Vs[iat2][{iat1, { 0,0,0 }}].transpose();
                    continue;
                }
                const size_t nabf2 = Cs.at(iat2).at({ 0, {0,0,0} }).shape[0];
                ifs >> nabf >> istart >> iend >> jstart >> jend >> tmp /*ik*/ >> tmp/*wk*/;
                assert(nabf1 == iend - istart + 1);
                assert(nabf2 == jend - jstart + 1);
                RI::Tensor<TR> t({ nabf1, nabf2 });
                for (int i = 0;i < nabf1;++i)
                {
                    for (int j = 0;j < nabf2;++j)
                    {
                        // t(i, j) = Vq[(istart + i) * nabf + jstart + j];
                        ifs >> t(i, j) >> tmp;
                    }
                }
                Vs[iat1][{iat2, { 0,0,0 }}] = t;
            }
        }
        return Vs;
    }

    template < typename TR> // any blocking
    TLRI<TR> read_coulomb_mat_general(const std::string& file, const TLRI<TR>& Cs)
    {   //for gamma_only, V(q)=V(R=0)
        std::ifstream ifs;
        ifs.open(file);
        size_t nks = 0, nabf = 0, istart = 0, jstart = 0, iend = 0, jend = 0;
        std::string tmp;
        ifs >> nks;//   nkstot=1
        if (nks > 1) { std::cout << "Warning: nks>1 is not supported yet!" << std::endl; }
        TLRI<TR> Vs;
        std::vector<TR> Vq;
        while (ifs.peek() != EOF)
        {
            ifs >> nabf >> istart >> iend >> jstart >> jend >> tmp /*ik*/ >> tmp/*wk*/;
            if (ifs.peek() == EOF) { break; }
            if (Vq.empty()) { Vq.resize(nabf * nabf, 0.0); }
            for (int i = istart - 1;i < iend;++i)
            {
                for (int j = jstart - 1;j < jend;++j)
                {
                    ifs >> Vq[i * nabf + j] >> tmp;
                }
            }
        }
        const int nat = Cs.size();
        istart = 0;    // 
        for (int iat1 = 0;iat1 < nat;++iat1)
        {
            const size_t nabf1 = Cs.at(iat1).at({ 0, {0,0,0} }).shape[0];
            jstart = 0;
            for (int iat2 = 0;iat2 < nat;++iat2)
            {
                const size_t nabf2 = Cs.at(iat2).at({ 0, {0,0,0} }).shape[0];
                if (iat1 > iat2)
                {   // coulomb_mat has only the upper triangle part
                    Vs[iat1][{iat2, { 0,0,0 }}] = Vs[iat2][{iat1, { 0,0,0 }}].transpose();
                }
                else
                {
                    RI::Tensor<TR> t({ nabf1, nabf2 });
                    for (int i = 0;i < nabf1;++i)
                    {
                        for (int j = 0;j < nabf2;++j)
                        {
                            t(i, j) = Vq[(istart + i) * nabf + jstart + j];
                        }
                    }
                    Vs[iat1][{iat2, { 0,0,0 }}] = t;
                }
                jstart += nabf2;
            }
            assert(jstart == nabf);
            istart += nabf1;
        }
        assert(istart == nabf);
        return Vs;
    }

    template < typename TR>
    bool compare_Vs(const TLRI<TR>& Vs1, const TLRI<TR>& Vs2, const double thr)
    {
        for (auto& tmp1 : Vs1)
        {
            const int& iat1 = tmp1.first;
            for (auto& tmp2 : tmp1.second)
            {
                const int& iat2 = tmp2.first.first;
                const RI::Tensor<TR>& t1 = tmp2.second;
                const RI::Tensor<TR>& t2 = Vs2.at(iat1).at({ iat2, {0,0,0} });
                if (t1.shape[0] != t2.shape[0]) { return false; }
                if (t1.shape[1] != t2.shape[1]) { return false; }
                for (int i = 0;i < t1.shape[0];++i)
                {
                    for (int j = 0;j < t1.shape[1];++j)
                    {
                        if (std::abs(t1(i, j) - t2(i, j)) > thr) { std::cout << "element (" << i << ", " << j << ") are differernt: " << t1(i, j) << ", " << t2(i, j) << std::endl;return false; }
                    }
                }
            }
        }
        return true;
    }
    template <typename TR>
    std::vector<TLRI<TR>> split_Ds(const std::vector<std::vector<TR>>& Ds, const std::vector<int>& aims_nbasis, const UnitCell& ucell)
    {
        // Due to the hard-coded constructor of elecstate::DensityMatrix, singlet-triplet with nspin=2 cannot use DM_trans with size 1
        // if(Ds.size()>1) { throw std::runtime_error("split_Ds only supports gamma-only spin-1 Ds now."); }
        std::vector<TLRI<TR>> Ds_split;
        for (const auto& D : Ds)
        {
            TLRI<TR> D_split;
            const int nbasis = std::sqrt(D.size());
            int iw1_start = 0;
            for (int iat1 = 0;iat1 < ucell.nat;++iat1)
            {
                const int& it1 = ucell.iat2it[iat1];
                const size_t& nw1 = aims_nbasis[it1];
                int iw2_start = 0;
                for (int iat2 = 0;iat2 < ucell.nat;++iat2)
                {
                    const int& it2 = ucell.iat2it[iat2];
                    const size_t& nw2 = aims_nbasis[it2];
                    D_split[iat1][{iat2, { 0,0,0 }}] = RI::Tensor<TR>({ nw1, nw2 });
                    for (int i = 0;i < nw1;++i)
                    {
                        for (int j = 0;j < nw2;++j)
                        {
                            D_split[iat1][{iat2, { 0,0,0 }}](i, j) = D[(iw1_start + i)*nbasis+(iw2_start + j)] * 0.5; // consistent with split_m2D_ktoR
                        }
                    }
                    iw2_start += nw2;
                }
                assert(iw2_start == nbasis);
                iw1_start += nw1;
            }
            assert(iw1_start == nbasis);
            Ds_split.push_back(D_split);
        }
        return Ds_split;
    }
}

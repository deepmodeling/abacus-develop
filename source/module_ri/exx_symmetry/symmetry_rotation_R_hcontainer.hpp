#include "symmetry_rotation.h"
#include "module_base/blas_connector.h"
#include "module_base/parallel_reduce.h"
namespace ModuleSymmetry
{
    template<typename TR>
    inline void print_global(const TR* tlocal, const int nr, const int nc, const std::string name)
    {
        GlobalV::ofs_running << name << std::endl;
        for (int i = 0;i < nr;++i)
        {
            for (int j = 0;j < nc;++j) GlobalV::ofs_running << tlocal[j + i * nc] << " ";
            GlobalV::ofs_running << "\n";
        }
    }
    template<typename TR>
    inline void print_local(const Parallel_Orbitals& pv, const TR* tlocal, const std::string name)
    {
        GlobalV::ofs_running << name << std::endl;
        for (int i = 0;i < pv.get_row_size();++i)
        {
            for (int j = 0;j < pv.get_col_size();++j) GlobalV::ofs_running << tlocal[j + i * pv.get_col_size()] << " ";
            GlobalV::ofs_running << "\n";
        }
    }
    template<typename TR>
    inline void print_atompair_local(const Parallel_Orbitals& pv, const int iat1, const int iat2, const TR* tlocal, const std::string name)
    {
        GlobalV::ofs_running << name << std::endl;
        for (int i = 0;i < pv.get_row_size(iat1);++i)
        {
            for (int j = 0;j < pv.get_col_size(iat2);++j) GlobalV::ofs_running << tlocal[j + i * pv.get_col_size(iat2)] << " ";
            GlobalV::ofs_running << "\n";
        }
    }

    template<typename TR>   // HContainer type
    void Symmetry_rotation::restore_HR(
        const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
        const hamilt::HContainer<TR>& HR_irreduceble,
        hamilt::HContainer<TR>& HR_rotated)const
    {
        ModuleBase::TITLE("Symmetry_rotation", "restore_HR");
        ModuleBase::timer::tick("Symmetry_rotation", "restore_HR");
        for (auto& apR_isym_irapR : this->irs_.full_map_to_irreducible_sector_)
        {
            const Tap& ap = apR_isym_irapR.first.first;
            const TC& R = apR_isym_irapR.first.second;
            const int& isym = apR_isym_irapR.second.first;
            const Tap& irap = apR_isym_irapR.second.second.first;
            const TC& irR = apR_isym_irapR.second.second.second;
            assert(irR == this->irs_.rotate_R(symm, isym, ap.first, ap.second, R));
            // get in and out pointer from HContainer
            const hamilt::AtomPair<TR>& irap_hc = HR_irreduceble.get_atom_pair(irap.first, irap.second);
            const int irR_hc = irap_hc.find_R(irR[0], irR[1], irR[2]);
            if (irR_hc < 0) continue;
            const TR* irijR_ptr = irap_hc.get_pointer(irR_hc);
            const hamilt::AtomPair<TR>& ap_hc = HR_rotated.get_atom_pair(ap.first, ap.second);
            const int R_hc = ap_hc.find_R(R[0], R[1], R[2]);
            if (R_hc < 0) continue;
            TR* ijR_ptr = ap_hc.get_pointer(R_hc);
            rotate_atompair_parallel(irijR_ptr, isym, atoms, st, irap, ap, mode, *irap_hc.get_paraV(), ijR_ptr);
        }
        ModuleBase::timer::tick("Symmetry_rotation", "restore_HR");
    }

    inline void set_block(const int starti, const int startj, const int nr, const ModuleBase::ComplexMatrix& block, double* obj)
    {   // row-major (col-contiguous)
        for (int i = 0;i < block.nr;++i)
            for (int j = 0;j < block.nc;++j)
                obj[starti + i + (startj + j) * nr] = block(j, i).real();
    };
    inline void set_block(const int starti, const int startj, const int nr, const ModuleBase::ComplexMatrix& block, std::complex<double>* obj)
    {   // row-major (col-contiguous)
        for (int i = 0;i < block.nr;++i)
            for (int j = 0;j < block.nc;++j)
                obj[starti + i + (startj + j) * nr] = block(j, i);
    };
    template<typename TR>
    void Symmetry_rotation::rotate_atompair_parallel(const TR* Alocal_in, const int isym, const Atom* atoms, const Statistics& st,
        const Tap& ap_in, const Tap& ap_out, const char mode, const Parallel_Orbitals& pv, TR* Alocal_out, const bool output)const
    {
        // all the matrices are row-major (col-contiguous)
        int iat1 = ap_in.first, iat2 = ap_in.second;
        int it1 = st.iat2it[iat1], it2 = st.iat2it[iat2];
        int ia1 = st.iat2ia[iat1], ia2 = st.iat2ia[iat2];
        // contruct T matrix
        std::vector<TR> T1, T2;
        auto set_rotation_matrix = [&](const int& it, const int& ia, std::vector<TR>& T)->int
            {
                T.resize(atoms[it].nw * atoms[it].nw, 0.0);
                const int iwstart = atoms[it].stapos_wf + ia * atoms[it].nw;
                int iw = 0;
                while (iw < atoms[it].nw)
                {
                    int l = atoms[it].iw2l[iw];
                    int nm = 2 * l + 1;
                    // this->set_block_to_mat2d(iwstart, iwstart, this->rotmat_Slm_[isym][l], T_full_2d, pv);
                    set_block(iw, iw, atoms[it].nw, this->rotmat_Slm_[isym][l], T.data());
                    iw += nm;
                }
                return iwstart;
            };
        int iw1start = set_rotation_matrix(it1, ia1, T1);
        int iw2start = set_rotation_matrix(it2, ia2, T2);

        // copy Aocal_in to a global atom-pair matrix
        std::vector<TR> A(atoms[it1].nw * atoms[it2].nw, 0.0);

        int abr1 = pv.atom_begin_row[iat1], abc2 = pv.atom_begin_col[iat2];
        if (abr1 >= 0 && abc2 >= 0)
        { // "pv.local2global_row(i) - iw1start": global index in current atom pair
            for (int j = 0;j < pv.get_col_size(iat2);++j)
                for (int i = 0;i < pv.get_row_size(iat1);++i)
                    A[(pv.local2global_row(i + abr1) - iw1start) * atoms[it2].nw + (pv.local2global_col(j + abc2) - iw2start)]
                    = Alocal_in[j + i * pv.get_col_size(iat2)];
        }
        if (output) print_global(A.data(), atoms[it1].nw, atoms[it2].nw, "A before allreduce");
        Parallel_Reduce::reduce_all(A.data(), A.size());

        // rotate
        const char notrans = 'N', transpose = 'T', dagger = 'C';
        // const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        std::vector<TR> AT2(atoms[it1].nw * atoms[it2].nw, 0.0);
        std::vector<TR> TAT(atoms[it1].nw * atoms[it2].nw, 0.0);
        if (mode == 'H')
        {   // H'=T1^\dagger * H * T2
            BlasConnector::gemm(notrans, notrans, atoms[it1].nw, atoms[it2].nw, atoms[it2].nw,
                1.0, A.data(), atoms[it2].nw, T2.data(), atoms[it2].nw, 0.0, AT2.data(), atoms[it2].nw);
            BlasConnector::gemm(dagger, notrans, atoms[it1].nw, atoms[it2].nw, atoms[it1].nw,
                1.0, T1.data(), atoms[it1].nw, AT2.data(), atoms[it2].nw, 0.0, TAT.data(), atoms[it2].nw);
        }
        else if (mode == 'D')
        {   // D' = T1^T * D * T2^* = T1^T * [T2^\dagger * D^T]^T
            BlasConnector::gemm(dagger, transpose, atoms[it2].nw, atoms[it1].nw, atoms[it2].nw,
                1.0, T2.data(), atoms[it2].nw, A.data(), atoms[it2].nw, 0.0, AT2.data(), atoms[it1].nw);
            BlasConnector::gemm(transpose, transpose, atoms[it1].nw, atoms[it2].nw, atoms[it1].nw,
                1.0, T1.data(), atoms[it1].nw, AT2.data(), atoms[it1].nw, 0.0, TAT.data(), atoms[it2].nw);
        }
        else throw std::invalid_argument("Symmetry_rotation::rotate_atompair_tensor: invalid mode.");

        if (output)
        {
            print_global(A.data(), atoms[it1].nw, atoms[it2].nw, "A");
            print_global(T1.data(), atoms[it1].nw, atoms[it1].nw, "T1");
            print_global(TAT.data(), atoms[it1].nw, atoms[it2].nw, "TAT");
            print_atompair_local(pv, iat1, iat2, Alocal_out, "Alocal_out");
        }

        // copy back to Alocal_out
        iat1 = ap_out.first, iat2 = ap_out.second;
        it1 = st.iat2it[iat1], it2 = st.iat2it[iat2];
        ia1 = st.iat2ia[iat1], ia2 = st.iat2ia[iat2];
        abr1 = pv.atom_begin_row[iat1], abc2 = pv.atom_begin_col[iat2];
        iw1start = atoms[it1].stapos_wf + ia1 * atoms[it1].nw;
        iw2start = atoms[it2].stapos_wf + ia2 * atoms[it2].nw;
        if (abr1 >= 0 && abc2 >= 0)
        {// ap_in index for TAT but ap_out index for Alocal_out
            for (int j = 0;j < pv.get_col_size(iat2);++j)
                for (int i = 0;i < pv.get_row_size(iat1);++i)
                    Alocal_out[j + i * pv.get_col_size(iat2)]
                    = TAT[(pv.local2global_row(i + abr1) - iw1start) * atoms[it2].nw + (pv.local2global_col(j + abc2) - iw2start)];
        }
    }

    template<typename TR>
    void Symmetry_rotation::test_HR_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st,
        const char mode, const hamilt::HContainer<TR>& HR_full)
    {
        ModuleBase::TITLE("Symmetry_rotation", "test_HR_rotation");
        auto get_irreducible_ijR_info = [&HR_full, this]() -> std::vector<int>
            {
                std::vector<int> irreducible_ijR_info;
                irreducible_ijR_info.push_back(this->irs_.irreducible_sector_.size());
                for (auto& irap_irR : this->irs_.irreducible_sector_)
                {
                    const int iat1 = irap_irR.first.first, iat2 = irap_irR.first.second;
                    irreducible_ijR_info.insert(irreducible_ijR_info.end(), { iat1, iat2, 0 });
                    int nR = 0;
                    for (auto& R : irap_irR.second)
                        if (HR_full.get_atom_pair(iat1, iat2).find_R(R[0], R[1], R[2]) >= 0)
                        {
                            irreducible_ijR_info.insert(irreducible_ijR_info.end(), R.begin(), R.end());
                            ++nR;
                        }
                    irreducible_ijR_info[irreducible_ijR_info.size() - 3 * nR - 1] = nR;
                }
                return irreducible_ijR_info;
            };

        const Parallel_Orbitals* pv = HR_full.get_atom_pair(0, 0).get_paraV();
        // 1. pick out H(R) in the irreducible sector from full H(R)
        const std::vector<int>& irreducible_ijR_info = get_irreducible_ijR_info();
        hamilt::HContainer<TR> HR_irreducible(pv, nullptr, &irreducible_ijR_info);
        HR_irreducible.set_zero();
        for (auto& irap_irR : this->irs_.irreducible_sector_)
        {
            const int iat1 = irap_irR.first.first, iat2 = irap_irR.first.second;
            const hamilt::AtomPair<TR>& irap = HR_irreducible.get_atom_pair(iat1, iat2);
            const hamilt::AtomPair<TR>& ap_full = HR_full.get_atom_pair(iat1, iat2);
            for (auto& R : irap_irR.second)
                if (irap.find_R(R[0], R[1], R[2]) >= 0)
                { // out-of-range R can be added to irreducible sector to be calculated, but not in this test for lack of reference
                    TR* irptr = irap.get_HR_values(R[0], R[1], R[2]).get_pointer();
                    TR* ptr = ap_full.get_HR_values(R[0], R[1], R[2]).get_pointer();
                    std::copy(ptr, ptr + pv->get_row_size(iat1) * pv->get_col_size(iat2), irptr);
                }
        }

        //2. rotate
        hamilt::HContainer<TR> HR_rotated(HR_full);
        HR_rotated.set_zero();
        this->restore_HR(symm, atoms, st, mode, HR_irreducible, HR_rotated);
        //3. compare
        for (int iat1 = 0;iat1 < st.nat;++iat1)
        {
            for (int iat2 = 0;iat2 < st.nat;++iat2)
            {
                const hamilt::AtomPair<TR>& ap_full = HR_full.get_atom_pair(iat1, iat2);
                const hamilt::AtomPair<TR>& ap_rotated = HR_rotated.get_atom_pair(iat1, iat2);
                assert(ap_full.get_R_size() == ap_rotated.get_R_size());
                for (int irR = 0;irR < ap_full.get_R_size();++irR)
                {
                    const TR* full_ptr = ap_full.get_pointer(irR);
                    int* R = ap_full.get_R_index(irR);
                    const TR* rotated_ptr = ap_rotated.get_HR_values(R[0], R[1], R[2]).get_pointer();
                    GlobalV::ofs_running << "atom pair: (" << iat1 << ", " << iat2 << "),  R: (" << R[0] << " " << R[1] << " " << R[2] << ")\n";
                    print_atompair_local(*pv, iat1, iat2, rotated_ptr, std::string("R_rot").insert(0, 1, mode));
                    print_atompair_local(*pv, iat1, iat2, full_ptr, std::string("R_ref").insert(0, 1, mode));
                }
            }
        }
    }
}
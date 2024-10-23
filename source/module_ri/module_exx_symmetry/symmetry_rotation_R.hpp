#include "symmetry_rotation.h"
#include "module_ri/RI_Util.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include <array>
#include <RI/global/Global_Func-2.h>
#include <RI/physics/symmetry/Symmetry_Rotation.h>
namespace ModuleSymmetry
{
    template<typename Tdata>
    inline void print_tensor(const RI::Tensor<Tdata>& t, const std::string& name, const double& threshold = 0.0)
    {
        GlobalV::ofs_running << name << ":\n";
        for (int i = 0;i < t.shape[0];++i)
        {
            for (int j = 0;j < t.shape[1];++j) {
                GlobalV::ofs_running << ((std::abs(t(i, j)) > threshold) ? t(i, j) : static_cast<Tdata>(0)) << " ";
}
            GlobalV::ofs_running << std::endl;
        }
    }
    template<typename Tdata>
    inline void print_tensor3(const RI::Tensor<Tdata>& t, const std::string& name, const double& threshold = 0.0)
    {
        GlobalV::ofs_running << name << ":\n";
        for (int a = 0;a < t.shape[0];++a)
        {
            GlobalV::ofs_running << "abf: " << a << '\n';
            for (int i = 0;i < t.shape[1];++i)
            {
                for (int j = 0;j < t.shape[2];++j) {
                    GlobalV::ofs_running << ((std::abs(t(a, i, j)) > threshold) ? t(a, i, j) : static_cast<Tdata>(0)) << " ";
}
                GlobalV::ofs_running << std::endl;
            }
            GlobalV::ofs_running << std::endl;
        }
    }

    template<typename Tdata>
    std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> Symmetry_rotation::restore_HR(
        const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
        const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_irreducible) const
    {
        ModuleBase::TITLE("Symmetry_rotation", "restore_HR");
        ModuleBase::timer::tick("Symmetry_rotation", "restore_HR");
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_full;
        // irreducile-to-full map ver.
        for (auto& tmp1 : HR_irreducible)
        {
            const int& irap1 = tmp1.first;
            for (auto& tmp2 : tmp1.second)
            {
                const int& irap2 = tmp2.first.first;
                const Tap& irap = { irap1, irap2 };
                const TC& irR = tmp2.first.second;
                const TapR& irapR = { irap, irR };
                if (this->irs_.sector_stars_.find(irapR) != this->irs_.sector_stars_.end())
                {
                    for (auto& isym_apR : this->irs_.sector_stars_.at(irapR))
                    {
                        const int& isym = isym_apR.first;
                        const TapR& apR = isym_apR.second;
                        const int& ap1 = apR.first.first;
                        const int& ap2 = apR.first.second;
                        const TC& R = apR.second;
                        HR_full[ap1][{ap2, R}] = rotate_atompair_serial(tmp2.second, isym, atoms[st.iat2it[irap1]], atoms[st.iat2it[irap2]], mode);
                    }
                }
                else { std::cout << "Warning: not found: irreducible atom pair =(" << irap1 << "," << irap2 << "), irR=(" << irR[0] << "," << irR[1] << "," << irR[2] << ")\n";}
            }
        }
        // full-to-irreducible map ver. (problematic in parallel)
        // openmp slows down this for loop, why?
        // for (auto& apR_isym_irapR : this->irs_.full_map_to_irreducible_sector_)
        // {
        //     const Tap& ap = apR_isym_irapR.first.first;
        //     const TC& R = apR_isym_irapR.first.second;
        //     const int& isym = apR_isym_irapR.second.first;
        //     const Tap& irap = apR_isym_irapR.second.second.first;
        //     const TC& irR = apR_isym_irapR.second.second.second;
        //     // rotate the matrix and pack data
        //     // H_12(R)=T^\dagger(V)H_1'2'(VR+O_1-O_2)T(V)
        //     if (HR_irreducible.find(irap.first) != HR_irreducible.end() && HR_irreducible.at(irap.first).find({ irap.second, irR }) != HR_irreducible.at(irap.first).end())
        //         HR_full[ap.first][{ap.second, R}] = rotate_atompair_serial(HR_irreducible.at(irap.first).at({ irap.second, irR }),
        //             isym, atoms[st.iat2it[irap.first]], atoms[st.iat2it[irap.second]], mode);
        //     else
        //         std::cout << "not found: current atom pair =(" << ap.first << "," << ap.second << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "), irreducible atom pair =(" << irap.first << "," << irap.second << "), irR=(" << irR[0] << "," << irR[1] << "," << irR[2] << ")\n";
        // }
        // // test: output HR_irreducile 
        // for (auto& tmp1 : HR_irreducible)
        // {
        //     const int& a1 = tmp1.first;
        //     for (auto& tmp2 : tmp1.second)
        //     {
        //         const int& a2 = tmp2.first.first;
        //         const TC& R = tmp2.first.second;
        //         print_tensor(tmp2.second, "HR_irreducible (" + std::to_string(a1) + ", " + std::to_string(a2) + "), R=(" + std::to_string(R[0]) + " " + std::to_string(R[1]) + " " + std::to_string(R[2]) + ")");
        //     }
        // }
        // // test: output HR 
        // for (auto& tmp1 : HR_full)
        // {
        //     const int& a1 = tmp1.first;
        //     for (auto& tmp2 : tmp1.second)
        //     {
        //         const int& a2 = tmp2.first.first;
        //         const TC& R = tmp2.first.second;
        //         print_tensor(tmp2.second, "HR_full (" + std::to_string(a1) + ", " + std::to_string(a2) + "), R=(" + std::to_string(R[0]) + " " + std::to_string(R[1]) + " " + std::to_string(R[2]) + ")");
        //     }
        // }
        ModuleBase::timer::tick("Symmetry_rotation", "restore_HR");
        return HR_full;
    }

    template<typename Tdata>
    inline void set_block(const int starti, const int startj, const RI::Tensor<std::complex<double>>& block,
        RI::Tensor<Tdata>& obj_tensor)
    {   // no changing row/col order
        for (int i = 0;i < block.shape[0];++i) {
            for (int j = 0;j < block.shape[1];++j) {
                obj_tensor(starti + i, startj + j) = RI::Global_Func::convert<Tdata>(block(i, j));
}
}
    }

    template<typename Tdata>
    RI::Tensor<Tdata> Symmetry_rotation::set_rotation_matrix(const Atom& a, const int& isym)const
    {
        RI::Tensor<Tdata> T({ static_cast<size_t>(a.nw), static_cast<size_t>(a.nw) }); // check if zero
        int iw = 0;
        while (iw < a.nw)
        {
            int l = a.iw2l[iw];
            int nm = 2 * l + 1;
            set_block(iw, iw, this->rotmat_Slm_[isym][l], T);
            iw += nm;
        }
        return T;
    }
    template<typename Tdata>
    RI::Tensor<Tdata> Symmetry_rotation::rotate_atompair_serial(const RI::Tensor<Tdata>& A, const int isym,
        const Atom& a1, const Atom& a2, const char mode, const bool output)const
    {   // due to col-contiguous, actually what we know is T^T and H^T (or D^T), 
        // and what we calculate is(H'^T = T ^ T * H ^ T * T^*) or (D'^T = T ^ \dagger * D ^ T * T)
        assert(mode == 'H' || mode == 'D');
        bool sametype = (a1.label == a2.label);
        assert(A.shape[0] == a1.nw);//col
        assert(A.shape[1] == a2.nw);//row
        // contrut T matrix 
        const RI::Tensor<Tdata>& T1 = this->set_rotation_matrix<Tdata>(a1, isym);
        const RI::Tensor<Tdata>& T2 = sametype ? T1 : this->set_rotation_matrix<Tdata>(a2, isym);
        // rotate
        RI::Tensor<Tdata>TAT(A.shape);
        (mode == 'H') ? RI::Sym::T1_HR_T2(TAT.ptr(), A.ptr(), T1, T2) : RI::Sym::T1_DR_T2(TAT.ptr(), A.ptr(), T1, T2);
        if (output)
        {
            print_tensor(A, "A");
            print_tensor(T1, "T1");
            print_tensor(T2, "T2");
            print_tensor(TAT, "TAT");
        }
        return TAT;
    }
    template<typename Tdata>
    void Symmetry_rotation::rotate_atompair_serial(Tdata* TAT, const Tdata* A,
        const int& nw1, const int& nw2, const int isym,
        const Atom& a1, const Atom& a2, const char mode)const
    {   // due to col-contiguous, actually what we know is T^T and H^T (or D^T), 
        // and what we calculate is(H'^T = T ^ T * H ^ T * T^*) or (D'^T = T ^ \dagger * D ^ T * T)
        assert(mode == 'H' || mode == 'D');
        bool sametype = (a1.label == a2.label);
        assert(nw1 == a1.nw);//col
        assert(nw2 == a2.nw);//row
        // contrut T matrix 
        const RI::Tensor<Tdata>& T1 = this->set_rotation_matrix<Tdata>(a1, isym);
        const RI::Tensor<Tdata>& T2 = sametype ? T1 : this->set_rotation_matrix<Tdata>(a2, isym);
        // rotate
        (mode == 'H') ? RI::Sym::T1_HR_T2(TAT, A, T1, T2) : RI::Sym::T1_DR_T2(TAT, A, T1, T2);
    }


    template<typename Tdata>
    RI::Tensor<Tdata> Symmetry_rotation::set_rotation_matrix_abf(const int& type, const int& isym)const
    {
        int  nabfs = 0;
        for (int l = 0;l < this->abfs_l_nchi_[type].size();++l) {nabfs += this->abfs_l_nchi_[type][l] * (2 * l + 1);
}
        RI::Tensor<Tdata> T({ static_cast<size_t>(nabfs), static_cast<size_t>(nabfs) }); // check if zero
        int iw = 0;
        for (int L = 0;L < this->abfs_l_nchi_[type].size();++L)
        {
            int nm = 2 * L + 1;
            for (int N = 0;N < this->abfs_l_nchi_[type][L];++N)
            {
                set_block(iw, iw, this->rotmat_Slm_[isym][L], T);
                iw += nm;
                // std::cout << "L=" << L << ", N=" << N << ", iw=" << iw << "\n";
            }
        }
        assert(iw == nabfs);
        return T;
    }

    template<typename Tdata>
    RI::Tensor<Tdata> Symmetry_rotation::rotate_singleC_serial(const RI::Tensor<Tdata>& C,
        const int isym, const Atom& a1, const Atom& a2, const int& type1, bool output)const
    {
        assert(this->reduce_Cs_);
        RI::Tensor<Tdata> Cout(C.shape);
        assert(C.shape.size() == 3);
        const int& slice_size = C.shape[1] * C.shape[2];
        // step 1: multiply 2 AOs' rotation matrices
        for (int iabf = 0;iabf < C.shape[0];++iabf) {
            this->rotate_atompair_serial(Cout.ptr() + iabf * slice_size, C.ptr() + iabf * slice_size,
                a1.nw, a2.nw, isym, a1, a2, 'H');
}
        // step 2: multiply the ABFs' rotation matrix from the left
        const RI::Tensor<Tdata>& Tabfs = this->set_rotation_matrix_abf<Tdata>(type1, isym);
        RI::Sym::T1_HR(Cout.ptr(), Cout.ptr(), Tabfs, slice_size);
        return Cout;
    }

    template<typename Tdata>
    void Symmetry_rotation::print_HR(const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR, const std::string name, const double& threshold)
    {
        for (auto& HR_ia1 : HR)
        {
            int iat1 = HR_ia1.first;
            for (auto& HR_ia12R : HR_ia1.second)
            {
                int iat2 = HR_ia12R.first.first;
                TC R = HR_ia12R.first.second;
                const RI::Tensor<Tdata>& HR_tensor = HR_ia12R.second;
                std::cout << "atom pair (" << iat1 << ", " << iat2 << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "), ";
                print_tensor(HR_tensor, name, threshold);
            }
        }
    }

    template<typename Tdata>
    void Symmetry_rotation::test_HR_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
        const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_full)
    {
        ModuleBase::TITLE("Symmetry_rotation", "test_HR_rotation");

        // 1. pick out H(R) in the irreducible sector from full H(R)
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_irreducible;
        for (auto& irap_Rs : this->irs_.irreducible_sector_)
        {
            const Tap& irap = irap_Rs.first;
            for (auto& irR : irap_Rs.second)
            {
                const std::pair<int, TC> a2_irR = { irap.second, irR };
                HR_irreducible[irap.first][a2_irR] = (HR_full.at(irap.first).count(a2_irR) != 0) ?
                    HR_full.at(irap.first).at(a2_irR)
                    : RI::Tensor<Tdata>(HR_full.at(irap.first).begin()->second.shape);
            }
        }
        // 2. rotate
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_rotated = restore_HR(symm, atoms, st, mode, HR_irreducible);
        // 3. compare
        for (auto& HR_ia1 : HR_rotated)
        {
            int iat1 = HR_ia1.first;
            for (auto& HR_ia12R : HR_ia1.second)
            {
                int iat2 = HR_ia12R.first.first;
                TC R = HR_ia12R.first.second;
                const RI::Tensor<Tdata>& HR_rot = HR_ia12R.second;
                if (HR_full.at(iat1).count({ iat2, R }) == 0)// rot back but not found
                {
                    std::cout << "R_rot not found in atom pair (" << iat1 << ", " << iat2 << "):  R=(" << R[0] << "," << R[1] << "," << R[2] << "):\n";
                    continue;
                }
                const RI::Tensor<Tdata>& HR_ref = HR_full.at(iat1).at({ iat2, R });
                assert(HR_rot.shape[0] == HR_ref.shape[0]);
                assert(HR_rot.shape[1] == HR_ref.shape[1]);
                // output 
                std::cout << "atom pair (" << iat1 << ", " << iat2 << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "):\n";
                print_tensor(HR_rot, std::string("R_rot").insert(0, 1, mode));
                print_tensor(HR_ref, std::string("R_ref").insert(0, 1, mode));
            }
        }
    }

    template<typename Tdata>
    void Symmetry_rotation::test_Cs_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st,
        const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& Cs_full)const
    {
        for (auto& sector_pair : this->irs_.full_map_to_irreducible_sector_)
        {
            const TapR& apR = sector_pair.first;
            const int& isym = sector_pair.second.first;
            const TapR& irapR = sector_pair.second.second;
            // if (apR.first != irapR.first || apR.second != irapR.second)
            if (apR.first != irapR.first)
            {
                std::cout << "irapR=(" << irapR.first.first << "," << irapR.first.second << "), (" << irapR.second[0] << "," << irapR.second[1] << "," << irapR.second[2] << "):\n";
                std::cout << "apR=(" << apR.first.first << "," << apR.first.second << "), (" << apR.second[0] << "," << apR.second[1] << "," << apR.second[2] << "):\n";
                const RI::Tensor<Tdata>& Cs_ir = Cs_full.at(irapR.first.first).at({ irapR.first.second,irapR.second });
                const RI::Tensor<Tdata>& Cs_ref = Cs_full.at(apR.first.first).at({ apR.first.second,apR.second });
                const RI::Tensor<Tdata>& Cs_rot = this->rotate_singleC_serial(Cs_ir, isym,
                    atoms[st.iat2it[irapR.first.first]], atoms[st.iat2it[irapR.first.second]], irapR.first.first);
                print_tensor3(Cs_rot, "Cs_rot");
                print_tensor3(Cs_ref, "Cs_ref");
                print_tensor3(Cs_ir, "Cs_irreducible");
                exit(0);
            }
        }

    }

}
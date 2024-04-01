#include "symmetry_rotation.h"
#include "module_ri/RI_Util.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include <array>
#include <RI/global/Global_Func-2.h>

namespace ModuleSymmetry
{
    template<typename Tdata>
    std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> Symmetry_rotation::restore_HR(
        const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
        const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_irreduceble) const
    {
        ModuleBase::TITLE("Symmetry_rotation", "restore_HR");
        ModuleBase::timer::tick("Symmetry_rotation", "restore_HR");
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_full;
        // openmp slows down this for loop, why?
        for (auto& apR_isym_irapR : this->irs_.full_map_to_irreducible_sector_)
        {
            const Tap& ap = apR_isym_irapR.first.first;
            const TC& R = apR_isym_irapR.first.second;
            const int& isym = apR_isym_irapR.second.first;
            const Tap& irap = apR_isym_irapR.second.second.first;
            const TC& irR = apR_isym_irapR.second.second.second;
            // rotate the matrix and pack data
            // H_12(R)=T^\dagger(V)H_1'2'(VR+O_1-O_2)T(V)
            if (HR_irreduceble.find(irap.first) != HR_irreduceble.end() && HR_irreduceble.at(irap.first).find({ irap.second, irR }) != HR_irreduceble.at(irap.first).end())
                HR_full[ap.first][{ap.second, R}] = rotate_atompair_serial(HR_irreduceble.at(irap.first).at({ irap.second, irR }),
                    isym, atoms[st.iat2it[irap.first]], atoms[st.iat2it[irap.second]], mode);
            else
                std::cout << "not found: current atom pair =(" << ap.first << "," << ap.second << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "), irreducible atom pair =(" << irap.first << "," << irap.second << "), irR=(" << irR[0] << "," << irR[1] << "," << irR[2] << ")\n";
        }
        ModuleBase::timer::tick("Symmetry_rotation", "restore_HR");
        return HR_full;
    }

    template<typename Tdata>
    inline void print_tensor(const RI::Tensor<Tdata>& t, const std::string& name, const double& threshold = 0.0)
    {
        std::cout << name << ":\n";
        for (int i = 0;i < t.shape[0];++i)
        {
            for (int j = 0;j < t.shape[1];++j)
                std::cout << ((std::abs(t(i, j)) > threshold) ? t(i, j) : static_cast<Tdata>(0)) << " ";
            std::cout << std::endl;
        }
    }

    template<typename Tdata>
    inline void set_block(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
        RI::Tensor<Tdata>& obj_tensor)
    {   // no changing row/col order
        for (int i = 0;i < block.nr;++i)
            for (int j = 0;j < block.nc;++j)
                obj_tensor(starti + i, startj + j) = RI::Global_Func::convert<Tdata>(block(i, j));
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

    inline RI::Tensor<std::complex<double>> TAT_HR(const RI::Tensor<std::complex<double>>& A,
        const RI::Tensor<std::complex<double>>& T1, const RI::Tensor<std::complex<double>>& T2)
    {
        const char notrans = 'N', transpose = 'T', dagger = 'C';
        const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        // H'^T = T2^T * H^T * T1^*
        const int& nw2 = T2.shape[0], & nw1 = T1.shape[0];
        RI::Tensor<std::complex<double>> AT2(A.shape);
        zgemm_(&notrans, &notrans, &nw2, &nw1, &nw2, &alpha, T2.ptr(), &nw2, A.ptr(), &nw2, &beta, AT2.ptr(), &nw2);
        RI::Tensor<std::complex<double>> TAT(A.shape);
        zgemm_(&notrans, &dagger, &nw2, &nw1, &nw1, &alpha, AT2.ptr(), &nw2, T1.ptr(), &nw1, &beta, TAT.ptr(), &nw2);
        // row-maj version
        // BlasConnector::gemm(notrans, notrans, nw1, nw2, nw2,
        //     alpha, A_complex.ptr(), nw2, T2.ptr(), nw2, beta, AT2.ptr(), nw2);
        // BlasConnector::gemm(dagger, notrans, nw1, nw2, nw1,
        //     alpha, T1.ptr(), nw1, AT2.ptr(), nw2, beta, TAT.ptr(), nw2);
        return TAT;
    }
    inline RI::Tensor<double> TAT_HR(const RI::Tensor<double>& A,
        const RI::Tensor<double>& T1, const RI::Tensor<double>& T2)
    {
        const char notrans = 'N', transpose = 'T', dagger = 'C';
        const double alpha(1.0), beta(0.0);
        // H'^T = T2^T * H^T * T1^*
        const int& nw2 = T2.shape[0], & nw1 = T1.shape[0];
        RI::Tensor<double> AT2(A.shape);
        dgemm_(&notrans, &notrans, &nw2, &nw1, &nw2, &alpha, T2.ptr(), &nw2, A.ptr(), &nw2, &beta, AT2.ptr(), &nw2);
        RI::Tensor<double> TAT(A.shape);
        dgemm_(&notrans, &dagger, &nw2, &nw1, &nw1, &alpha, AT2.ptr(), &nw2, T1.ptr(), &nw1, &beta, TAT.ptr(), &nw2);
        return TAT;
    }
    inline RI::Tensor<std::complex<double>> TAT_DR(const RI::Tensor<std::complex<double>>& A,
        const RI::Tensor<std::complex<double>>& T1, const RI::Tensor<std::complex<double>>& T2)
    {
        const char notrans = 'N', transpose = 'T', dagger = 'C';
        const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        //T2^\dagger * D^T * T1 = [(D^T)^T * (T2^T)^\dagger]^T * (T1^T)^T
        const int& nw2 = T2.shape[0], & nw1 = T1.shape[0];
        RI::Tensor<std::complex<double>> AT2(A.shape);
        zgemm_(&transpose, &dagger, &nw1, &nw2, &nw2, &alpha, A.ptr(), &nw2, T2.ptr(), &nw2, &beta, AT2.ptr(), &nw1);
        RI::Tensor<std::complex<double>> TAT(A.shape);
        zgemm_(&transpose, &transpose, &nw2, &nw1, &nw1, &alpha, AT2.ptr(), &nw1, T1.ptr(), &nw1, &beta, TAT.ptr(), &nw2);
        return TAT;
    }
    inline RI::Tensor<double> TAT_DR(const RI::Tensor<double>& A,
        const RI::Tensor<double>& T1, const RI::Tensor<double>& T2)
    {
        const char notrans = 'N', transpose = 'T', dagger = 'C';
        const double alpha(1.0), beta(0.0);
        //T2^\dagger * D^T * T1 = [(D^T)^T * (T2^T)^\dagger]^T * (T1^T)^T
        const int& nw2 = T2.shape[0], & nw1 = T1.shape[0];
        RI::Tensor<double> AT2(A.shape);
        dgemm_(&transpose, &dagger, &nw1, &nw2, &nw2, &alpha, A.ptr(), &nw2, T2.ptr(), &nw2, &beta, AT2.ptr(), &nw1);
        RI::Tensor<double> TAT(A.shape);
        dgemm_(&transpose, &transpose, &nw2, &nw1, &nw1, &alpha, AT2.ptr(), &nw1, T1.ptr(), &nw1, &beta, TAT.ptr(), &nw2);
        return TAT;
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
        const RI::Tensor<Tdata>& TAT = (mode == 'H') ? TAT_HR(A, T1, T2) : TAT_DR(A, T1, T2);
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
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_irreduceble;
        for (auto& irap_Rs : this->irs_.irreducible_sector_)
        {
            const Tap& irap = irap_Rs.first;
            for (auto& irR : irap_Rs.second)
            {
                const std::pair<int, TC> a2_irR = { irap.second, irR };
                HR_irreduceble[irap.first][a2_irR] = (HR_full.at(irap.first).count(a2_irR) != 0) ?
                    HR_full.at(irap.first).at(a2_irR)
                    : RI::Tensor<Tdata>(HR_full.at(irap.first).begin()->second.shape);
            }
        }
        // 2. rotate
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_rotated = restore_HR(symm, atoms, st, mode, HR_irreduceble);
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

}
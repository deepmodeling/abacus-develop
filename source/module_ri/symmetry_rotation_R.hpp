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
        for (auto& apR_isym_irapR : this->full_map_to_irreducible_sector_)
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
    RI::Tensor<Tdata> Symmetry_rotation::rotate_atompair_serial(const RI::Tensor<Tdata>& A, const int isym,
        const Atom& a1, const Atom& a2, const char mode, const bool output)const
    {   // due to col-contiguous, actually what we know is T^T and H^T (or D^T), 
        // and what we calculate is(H'^T = T ^ T * H ^ T * T^*) or (D'^T = T ^ \dagger * D ^ T * T)
        auto set_block = [](const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
            RI::Tensor<std::complex<double>>& obj_tensor)->void
            {   // both ComplexMatrix and RI::Tensor are row-major (col-contiguous)
                for (int i = 0;i < block.nr;++i)
                    for (int j = 0;j < block.nc;++j)
                        obj_tensor(starti + i, startj + j) = block(i, j);
            };
        auto set_rotation_matrix = [&, this](const Atom& a) -> RI::Tensor<std::complex<double>>
            {
                RI::Tensor<std::complex<double>> T({ static_cast<size_t>(a.nw), static_cast<size_t>(a.nw) }); // check if zero
                int iw = 0;
                while (iw < a.nw)
                {
                    int l = a.iw2l[iw];
                    int nm = 2 * l + 1;
                    set_block(iw, iw, this->rotmat_Slm_[isym][l], T);
                    iw += nm;
                }
                return T;
            };

        bool sametype = (a1.label == a2.label);
        assert(A.shape[0] == a1.nw);//col
        assert(A.shape[1] == a2.nw);//row
        // contrut T matrix 
        const RI::Tensor<std::complex<double>>& T1 = set_rotation_matrix(a1);
        const RI::Tensor<std::complex<double>>& T2 = sametype ? T1 : set_rotation_matrix(a2);

        // A*T_2 (atom 2 is contiguous)
        const char notrans = 'N', transpose = 'T', dagger = 'C';
        const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
        const RI::Tensor<std::complex<double>>& A_complex = RI::Global_Func::convert<std::complex<double>>(A);

        RI::Tensor<std::complex<double>> TAT(A.shape);
        RI::Tensor<std::complex<double>> AT2(A.shape);
        if (mode == 'H')
        {   // H'^T = T2^T * H^T * T1^*
            zgemm_(&notrans, &notrans, &a2.nw, &a1.nw, &a2.nw, &alpha, T2.ptr(), &a2.nw, A_complex.ptr(), &a2.nw, &beta, AT2.ptr(), &a2.nw);
            zgemm_(&notrans, &dagger, &a2.nw, &a1.nw, &a1.nw, &alpha, AT2.ptr(), &a2.nw, T1.ptr(), &a1.nw, &beta, TAT.ptr(), &a2.nw);
            // row-maj version
            // BlasConnector::gemm(notrans, notrans, a1.nw, a2.nw, a2.nw,
            //     alpha, A_complex.ptr(), a2.nw, T2.ptr(), a2.nw, beta, AT2.ptr(), a2.nw);
            // BlasConnector::gemm(dagger, notrans, a1.nw, a2.nw, a1.nw,
            //     alpha, T1.ptr(), a1.nw, AT2.ptr(), a2.nw, beta, TAT.ptr(), a2.nw);

        }
        else if (mode == 'D')
        {   //T2^\dagger * D^T * T1 = [(D^T)^T * (T2^T)^\dagger]^T * (T1^T)^T
            zgemm_(&transpose, &dagger, &a1.nw, &a2.nw, &a2.nw, &alpha, A_complex.ptr(), &a2.nw, T2.ptr(), &a2.nw, &beta, AT2.ptr(), &a1.nw);
            zgemm_(&transpose, &transpose, &a2.nw, &a1.nw, &a1.nw, &alpha, AT2.ptr(), &a1.nw, T1.ptr(), &a1.nw, &beta, TAT.ptr(), &a2.nw);
        }
        else throw std::invalid_argument("Symmetry_rotation::rotate_atompair_tensor: invalid mode.");
        if (output)
        {
            print_tensor(A, "A");
            print_tensor(T1, "T1");
            print_tensor(T2, "T2");
            print_tensor(TAT, "TAT");
        }
        return RI::Global_Func::convert<Tdata>(TAT);
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
        for (auto& irap_Rs : this->irreducible_sector_)
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
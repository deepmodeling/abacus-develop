#include "symmetry_rotation.h"
#include "module_ri/RI_Util.h"
#include "module_base/blas_connector.h"
#include <array>
#include <RI/global/Global_Func-2.h>

namespace ModuleSymmetry
{
    template<typename Tdata>
    std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> Symmetry_rotation::restore_HR(
        const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
        const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_irreduceble)
    {
        ModuleBase::TITLE("Symmetry_rotation", "restore_HR");
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> HR_full;

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
        return HR_full;
    }

    template<typename Tdata>
    inline void print_tensor(const RI::Tensor<Tdata>& t, const std::string& name)
    {
        std::cout << name << ":\n";
        for (int i = 0;i < t.shape[0];++i)
        {
            for (int j = 0;j < t.shape[1];++j)
                std::cout << t(i, j) << " ";
            std::cout << std::endl;
        }
    }

    template<typename Tdata>
    RI::Tensor<Tdata> Symmetry_rotation::rotate_atompair_serial(const RI::Tensor<Tdata>& A, const int isym,
        const Atom& a1, const Atom& a2, const char mode, const bool output)
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
    void Symmetry_rotation::print_HR(const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR, const std::string name)
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
                print_tensor(HR_tensor, name);
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

    template<typename Tdata>
    void Symmetry_rotation::test_sector_equivalence_in_cal_Hs(const TapR& apR_test,
        const std::vector<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>>& Ds_full,
        Exx_LRI<Tdata>& exx_lri, const Parallel_Orbitals& pv)
    {
        // 1. pick out D(R) on one of the irreducible {abR}s from full D(R)
        std::vector<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>> Ds_one_in_sector = Ds_full;
        for (auto& a1_a2R_tensor : Ds_one_in_sector[0])
        {
            const int& iat1 = a1_a2R_tensor.first;
            for (auto& a2R_tensor : a1_a2R_tensor.second)
            {
                const int& iat2 = a2R_tensor.first.first;
                const TC& R = a2R_tensor.first.second;
                const TapR& apR = { {iat1, iat2}, R };
                if (apR.first == apR_test.first && apR.second == apR_test.second)
                    a2R_tensor.second = Ds_full[0].at(iat1).at({ iat2, R });
                else
                    a2R_tensor.second = RI::Tensor<Tdata>(a2R_tensor.second.shape);
            }
        }
        exx_lri.cal_exx_elec(Ds_one_in_sector, pv);
        const std::vector<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>> Hexxs_one_in_sector = exx_lri.get_Hexxs();


        //2. pick out D(R) on the whole sector star of the tested irreducible {abR}
        TapR irapR_test = this->full_map_to_irreducible_sector_.at(apR_test).second;
        std::vector<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>> Ds_sector_star = Ds_full;
        for (auto& a1_a2R_tensor : Ds_sector_star[0])
        {
            const int& iat1 = a1_a2R_tensor.first;
            for (auto& a2R_tensor : a1_a2R_tensor.second)
            {
                const int& iat2 = a2R_tensor.first.first;
                const TC& R = a2R_tensor.first.second;
                const TapR& apR = { {iat1, iat2}, R };
                const TapR& irapR = this->full_map_to_irreducible_sector_.at(apR).second;
                if (irapR.first == irapR_test.first && irapR.second == irapR_test.second)
                    a2R_tensor.second = Ds_full[0].at(iat1).at({ iat2, R });
                else
                    a2R_tensor.second = RI::Tensor<Tdata>(a2R_tensor.second.shape);
            }
        }
        exx_lri.cal_exx_elec(Ds_sector_star, pv);
        const std::vector<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>> Hexxs_sector_star = exx_lri.get_Hexxs();

        // test: output Ds
        std::cout << "Ds_one_in_sector[0].size(): " << Ds_one_in_sector[0].size() << std::endl;
        for (auto& a1_a2R_tensor : Ds_one_in_sector[0])
        {
            const int& iat1 = a1_a2R_tensor.first;
            for (auto& a2R_tensor : a1_a2R_tensor.second)
            {
                const int& iat2 = a2R_tensor.first.first;
                const TC& R = a2R_tensor.first.second;
                std::cout << "atom pair (" << iat1 << ", " << iat2 << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "):\n";
                const RI::Tensor<Tdata>& Ds_one = a2R_tensor.second;
                const RI::Tensor<Tdata>& Ds_star = Ds_sector_star[0].at(iat1).at({ iat2, R });
                print_tensor(Ds_one, "Ds_one");
                print_tensor(Ds_star, "Ds_star");
            }
        }

        // 3. compare
        // get sector star size
        auto get_star_size = [this](const TapR& apR_test)->int
            {
                for (auto& ss : this->sector_stars_)
                {
                    const TapR& star_apR = ss.begin()->second;
                    const TapR& star_irapR = this->full_map_to_irreducible_sector_.at(star_apR).second;
                    if (star_apR.first == apR_test.first && star_apR.second == apR_test.second)
                        return ss.size();
                }
                return 0;
            };

        const int starsize = get_star_size(apR_test);
        std::cout << "star size of tested apR: " << starsize << std::endl;
        std::cout << "Hexxs_one_in_sector[0].size(): " << Hexxs_one_in_sector[0].size() << std::endl;
        for (auto& a1_a2R_tensor : Hexxs_one_in_sector[0])
        {
            const int& iat1 = a1_a2R_tensor.first;
            for (auto& a2R_tensor : a1_a2R_tensor.second)
            {
                const int& iat2 = a2R_tensor.first.first;
                const TC& R = a2R_tensor.first.second;
                std::cout << "atom pair (" << iat1 << ", " << iat2 << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "):\n";
                const RI::Tensor<Tdata>& Hexxs_one = a2R_tensor.second;
                const RI::Tensor<Tdata>& Hexxs_star = Hexxs_sector_star[0].at(iat1).at({ iat2, R });
                print_tensor(Hexxs_one * RI::Global_Func::convert<Tdata>(static_cast<double>(starsize)), "Hexxs_one* starsize");
                print_tensor(Hexxs_star, "Hexxs_star");
            }
        }
    }
}
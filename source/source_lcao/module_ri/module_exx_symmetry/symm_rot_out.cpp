#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include "./symmetry_rotation.h"
namespace ModuleSymmetry
{
    namespace
    {
        int round_near_integer(const double value, const std::string& context)
        {
            const double rounded = std::round(value);
            if (std::abs(value - rounded) > 1e-8)
            {
                throw std::runtime_error("Failed to round " + context + " to an integer in symmetry sidecar output.");
            }
            return static_cast<int>(rounded);
        }

        std::array<int, 3> round_near_integer_vec3(const TCdouble& vec, const std::string& context)
        {
            return {round_near_integer(vec.x, context + " x"),
                    round_near_integer(vec.y, context + " y"),
                    round_near_integer(vec.z, context + " z")};
        }

        std::string int_vec3_fmt(const std::array<int, 3>& vec)
        {
            return "(" + std::to_string(vec[0]) + " " + std::to_string(vec[1]) + " "
                   + std::to_string(vec[2]) + ")";
        }

        std::array<int, 3> build_kspace_fold_G(const TCdouble& k_ibz,
                                               const TCdouble& k_full,
                                               const int isym,
                                               const UnitCell& ucell)
        {
            const bool trs_conj = isym >= ucell.symm.nrotk;
            const int isym_space = trs_conj ? isym - ucell.symm.nrotk : isym;
            TCdouble rotated = k_full * ucell.symm.kgmatrix[isym_space];
            if (trs_conj)
            {
                rotated.x = -rotated.x;
                rotated.y = -rotated.y;
                rotated.z = -rotated.z;
            }

            TCdouble fold_G;
            fold_G.x = rotated.x - k_ibz.x;
            fold_G.y = rotated.y - k_ibz.y;
            fold_G.z = rotated.z - k_ibz.z;
            return round_near_integer_vec3(fold_G, "k-space fold G");
        }

        int count_basis_size(const std::vector<int>& shell_counts)
        {
            int nao = 0;
            for (int l = 0; l < static_cast<int>(shell_counts.size()); ++l)
            {
                nao += shell_counts[static_cast<std::size_t>(l)] * (2 * l + 1);
            }
            return nao;
        }

        std::vector<int> build_type_lmaxs(
            const std::vector<std::vector<std::vector<int>>>& abf_layout_candidates)
        {
            std::vector<int> type_lmaxs(abf_layout_candidates.size(), -1);
            for (std::size_t itype = 0; itype < abf_layout_candidates.size(); ++itype)
            {
                for (const auto& shell_counts : abf_layout_candidates[itype])
                {
                    type_lmaxs[itype] = std::max(type_lmaxs[itype],
                                                 static_cast<int>(shell_counts.size()) - 1);
                }
            }
            return type_lmaxs;
        }

        void write_ao_layout_header(std::ofstream& ofs, const UnitCell& ucell)
        {
            ofs << "AO shell layouts:\n";
            ofs << "# Each line stores one AO-basis layout for one atom type.\n";
            ofs << "# shell_counts[l] is the number of AO radial shells at angular momentum l.\n";
            for (int it = 0; it < ucell.ntype; ++it)
            {
                const auto& atom = ucell.atoms[it];
                ofs << "type " << it + 1
                    << " label " << atom.label
                    << " nao " << atom.nw
                    << " lmax " << atom.nwl
                    << " shell_counts";
                for (int l = 0; l <= atom.nwl; ++l)
                {
                    ofs << " " << atom.l_nchi[l];
                }
                ofs << "\n";
            }
            ofs << "End AO shell layouts\n";
        }

        void write_abf_layout_header(
            std::ofstream& ofs,
            const std::vector<std::string>& type_labels,
            const std::vector<std::vector<std::vector<int>>>& abf_layout_candidates)
        {
            if (type_labels.size() != abf_layout_candidates.size())
            {
                throw std::runtime_error("ABF shell-layout header is inconsistent with the atom-type labels.");
            }

            ofs << "ABF shell layouts:\n";
            ofs << "# Each line stores one auxiliary-basis candidate for one atom type.\n";
            ofs << "# shell_counts[l] is the number of auxiliary radial shells at angular momentum l.\n";
            for (std::size_t itype = 0; itype < abf_layout_candidates.size(); ++itype)
            {
                for (const auto& shell_counts : abf_layout_candidates[itype])
                {
                    ofs << "type " << itype + 1
                        << " label " << type_labels[itype]
                        << " nao " << count_basis_size(shell_counts)
                        << " lmax " << static_cast<int>(shell_counts.size()) - 1
                        << " shell_counts";
                    for (const int shell_count : shell_counts)
                    {
                        ofs << " " << shell_count;
                    }
                    ofs << "\n";
                }
            }
            ofs << "End ABF shell layouts\n";
        }

        void print_symrot_info_k_impl(const std::string& file_name,
                                      const Symmetry_rotation& symrot,
                                      const K_Vectors& kv,
                                      const UnitCell& ucell,
                                      const std::vector<int>* basis_lmaxs,
                                      const std::vector<std::string>* abf_type_labels,
                                      const std::vector<std::vector<std::vector<int>>>* abf_layout_candidates)
        {
            std::ofstream ofs(file_name);
            ofs << std::scientific << std::setprecision(15);
            if (abf_layout_candidates != nullptr && abf_type_labels != nullptr)
            {
                write_abf_layout_header(ofs, *abf_type_labels, *abf_layout_candidates);
            }
            else
            {
                write_ao_layout_header(ofs, ucell);
            }
            ofs << "Number of IBZ k-points (k stars): " << kv.kstars.size() << std::endl;
            ofs << "Format:\n" << "The symmetry operation index to the irreducible k-point. For the irreducible k-points, isym=0.\n\n"
                << "(The direct coordinate of the original k-point)\n"
                << "(The reciprocal-lattice fold G satisfying k_full * S = k_ibz + G; TRS members use -k_full * S = k_ibz + G)\n"
                << "For each atom: \n"
                << "- Original index->transformed index, type and the Lmax\n"
                << "- Exact return lattice O used by ABACUS in the Bloch phase\n"
                << "- Bloch orbital rotation matrix (M) of the given operation and atom, for each angular momentum\n\n";
            for (int istar = 0;istar < kv.kstars.size();++istar)
            {
                ofs << "Star " << istar + 1 << " of IBZ k-point " << vec3_fmt(kv.kstars[istar].at(0)) << ":\n";
                for (const auto& isym_kvd : kv.kstars[istar])
                {
                    const int& isym = isym_kvd.first;
                    const bool trs_conj = isym >= ucell.symm.nrotk;
                    const int isym_space = trs_conj ? isym - ucell.symm.nrotk : isym;
                    ofs << isym;
                    if (trs_conj)
                    {
                        ofs << " (TRS * " << isym_space << ")";
                    }
                    ofs << "\n" << vec3_fmt(isym_kvd.second) << "\n";
                    ofs << "fold_G = " << int_vec3_fmt(build_kspace_fold_G(kv.kstars[istar].at(0),
                                                                            isym_kvd.second,
                                                                            isym,
                                                                            ucell))
                        << "\n";
                    for (int iat1 =0;iat1 < ucell.nat;++iat1)
                    {
                        const int it = ucell.iat2it[iat1];
                        const int iat2 = ucell.symm.get_rotated_atom(isym_space, iat1);
                        const TCdouble return_lattice = symrot.get_return_lattice(iat1, isym_space);
                        const double arg = 2 * ModuleBase::PI * isym_kvd.second * return_lattice;
                        const std::complex<double> phase_factor = std::complex<double>(std::cos(arg), std::sin(arg));

                        int lmax = ucell.atoms[it].nwl;
                        if (basis_lmaxs != nullptr)
                        {
                            if (it >= basis_lmaxs->size())
                            {
                                throw std::runtime_error("ABF shell-count metadata is inconsistent with atom types");
                            }
                            lmax = basis_lmaxs->at(it);
                        }

                        ofs << "atom " << iat1 + 1 << " -> " << iat2 + 1 << " of type " << it + 1 << " with Lmax= " << lmax << "\n";
                        ofs << "return_lattice = "
                            << int_vec3_fmt(round_near_integer_vec3(return_lattice,
                                                                    "atom return lattice"))
                            << "\n";
                        for (int l = 0;l < lmax + 1;++l)
                        {
                            const int nm = 2 * l + 1;
                            const auto& m_block = symrot.rotmat_Slm[isym_space][l];
                            for (int m1 = 0;m1 < nm;++m1)
                            {
                                for (int m2 = 0;m2 < nm;++m2)
                                {
                                    ofs << phase_factor * m_block(m1, m2);
                                }
                                ofs << "\n";
                            }
                        }
                    }
                }
                ofs << "\n";
            }
            ofs.close();
        }
    }

    std::string mat3_fmt(const ModuleBase::Matrix3& m)
    {
        using ModuleSymmetry::scalar_fmt;
        return scalar_fmt(m.e11) + " " + scalar_fmt(m.e12) + " " + scalar_fmt(m.e13) + "\n"
               + scalar_fmt(m.e21) + " " + scalar_fmt(m.e22) + " " + scalar_fmt(m.e23) + "\n"
               + scalar_fmt(m.e31) + " " + scalar_fmt(m.e32) + " " + scalar_fmt(m.e33);
    }

    // needs to calculate Ts from l=0 to l=max(l_ao,l_abf) before

    void print_symrot_info_R(const Symmetry_rotation& symrot, const Symmetry& symm,
        const int lmax_ao, const std::vector<TC>& Rs)
    {
        ModuleBase::TITLE("ModuleSymmetry", "print_symrot_info_R");
        std::ofstream ofs(PARAM.globalv.global_out_dir + "symrot_R.txt");
        ofs << std::scientific << std::setprecision(15);
        // Print the irreducible sector (to be optimized)
        ofs << "Number of irreducible sector: " << symrot.get_irreducible_sector().size() << std::endl;
        ofs << "Lmax of AOs: " << lmax_ao << "\n";
        ofs << "Lmax of ABFs: " << symrot.abfs_Lmax << "\n";
        // print AO rotation matrix T
        ofs << "Format:\n"
            << "The index of the symmetry operation\n"
            << "The rotation matrix of this symmetry operation (3*3)\n"
            << "(The translation vector of this symmetry operation)\n"
            << "Orbital rotation matrix (T) of each angular momentum with size ((2l + 1) * (2l + 1)) \n\n";
        const int lmax = std::max(lmax_ao, symrot.abfs_Lmax);
        for (int isym = 0;isym < symm.nrotk;++isym)
        {
            ofs << isym << "\n" << mat3_fmt(symm.gmatrix[isym]) << "\n"
                << vec3_fmt(symm.gtrans[isym]) << "\n";
            for (int l=0;l <= lmax;++l)
            {
                const int nm = 2 * l + 1;
                // ofs << "l = " << l << ", nm = " << nm << "\n";
                const auto& T_block = symrot.rotmat_Slm[isym][l];
                for (int m1 = 0;m1 < nm;++m1)
                {
                    for (int m2 = 0;m2 < nm;++m2)
                    {
                        //note: the order of m in orbitals may be different from increasing
                        //note: is Ts row- or col-major ?
                        ofs << T_block(m1, m2);
                    }
                    ofs << "\n";
                }
            }
            }
        ofs.close();
    }

    void print_symrot_info_k(const Symmetry_rotation& symrot, const K_Vectors& kv, const UnitCell& ucell)
    {
        ModuleBase::TITLE("Symmetry_rotation", "print_symrot_info_k");
        print_symrot_info_k_impl(PARAM.globalv.global_out_dir + "symrot_k.txt",
                                 symrot, kv, ucell, nullptr, nullptr, nullptr);
        ModuleBase::timer::tick("Symmetry_rotation", "print_symrot_info_k");
    }

    void print_symrot_info_abf_k(const Symmetry_rotation& symrot,
                                 const K_Vectors& kv,
                                 const UnitCell& ucell,
                                 const std::vector<std::string>& type_labels,
                                 const std::vector<std::vector<std::vector<int>>>& abf_layout_candidates)
    {
        ModuleBase::TITLE("Symmetry_rotation", "print_symrot_info_abf_k");
        ModuleBase::timer::tick("Symmetry_rotation", "print_symrot_info_abf_k");
        const auto type_lmaxs = build_type_lmaxs(abf_layout_candidates);
        print_symrot_info_k_impl(PARAM.globalv.global_out_dir + "symrot_abf_k.txt",
                                 symrot, kv, ucell,
                                 &type_lmaxs, &type_labels, &abf_layout_candidates);
        ModuleBase::timer::tick("Symmetry_rotation", "print_symrot_info_abf_k");
    }
}

#include "./symmetry_rotation.h"
namespace ModuleSymmetry
{
    std::string mat3_fmt(const ModuleBase::Matrix3& m)
    {
        auto s = [](auto x) { return std::to_string(x); };
        return s(m.e11) + " " + s(m.e12) + " " + s(m.e13) + "\n" +
               s(m.e21) + " " + s(m.e22) + " " + s(m.e23) + "\n" +
               s(m.e31) + " " + s(m.e32) + " " + s(m.e33);
    }

    // needs to calculate Ts from l=0 to l=max(l_ao,l_abf) before

    void print_symrot_info_R(const Symmetry_rotation& symrot, const Symmetry& symm,
        const int lmax_ao, const std::vector<TC>& Rs)
    {
        ModuleBase::TITLE("ModuleSymmetry", "print_symrot_info_R");
        std::ofstream ofs(PARAM.globalv.global_out_dir + "symrot_R.txt");
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
        std::ofstream ofs(PARAM.globalv.global_out_dir + "symrot_k.txt");
        ofs << "Number of IBZ k-points (k stars): " << kv.kstars.size() << std::endl;
        ofs << "Format:\n" << "The symmetry operation index to the irreducible k-point. For the irreducible k-points, isym=0.\n\n"
            << "(The direct coordinate of the original k-point)\n"
            << "For each atom: \n"
            << "- Original index->transformed index, type and the Lmax\n"
            << "- Bloch orbital rotation matrix (M) of the given operation and atom, for each angular momentum\n\n";
        for (int istar = 0;istar < kv.kstars.size();++istar)
        {
            ofs << "Star " << istar + 1 << " of IBZ k-point " << vec3_fmt(kv.kstars[istar].at(0)) << ":\n";
            for (const auto& isym_kvd : kv.kstars[istar])
            {
                const int& isym = isym_kvd.first;
                ofs << isym << "\n" << vec3_fmt(isym_kvd.second) << "\n";
                for (int iat1 =0;iat1 < ucell.nat;++iat1)
                {
                    const int it = ucell.iat2it[iat1];  // it1=it2
                    const int lmax = ucell.atoms[it].nwl;
                    const int iat2 = ucell.symm.get_rotated_atom(isym, iat1);
                    const double arg = 2 * ModuleBase::PI * isym_kvd.second * symrot.get_return_lattice(iat1,isym);
                    std::complex<double>phase_factor = std::complex<double>(std::cos(arg), std::sin(arg));
                    ofs << "atom " << iat1 + 1 << " -> " << iat2 + 1 << " of type " << it + 1 << " with Lmax= " << lmax << "\n";
                    for (int l = 0;l < lmax + 1;++l)
                    {
                        const int nm = 2 * l + 1;
                        const auto& m_block = symrot.rotmat_Slm[isym][l];
                        for (int m1 = 0;m1 < nm;++m1)
                        {
                            // const int m1_start = m2 * nm;
                            for (int m2 = 0;m2 < nm;++m2)
                            {
                                ofs << phase_factor * m_block(m1, m2);    // row-major
                            }
                            ofs << "\n";
                        }
                    }// end l
                }   // end iat
            }   // end (k, op)
            ofs << "\n";
        }   // end star
        ofs.close();
        ModuleBase::timer::tick("Symmetry_rotation", "print_symrot_info_k");
    }
}
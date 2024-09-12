#include "irreducible_sector.h"
namespace ModuleSymmetry
{
    ModuleBase::Matrix3 Irreducible_Sector::direct_to_cartesian(const ModuleBase::Matrix3& d, const ModuleBase::Matrix3& latvec)const
    {
        return latvec.Inverse() * d * latvec;
    }

    std::vector<int> inline get_isymbvk_to_isym_map(const std::vector<ModuleBase::Matrix3>& bvkgmat, const ModuleSymmetry::Symmetry& symm)
    {
        auto matequal = [&symm](ModuleBase::Matrix3 a, ModuleBase::Matrix3 b)
            {
                return (symm.equal(a.e11, b.e11) && symm.equal(a.e12, b.e12) && symm.equal(a.e13, b.e13) &&
                    symm.equal(a.e21, b.e21) && symm.equal(a.e22, b.e22) && symm.equal(a.e23, b.e23) &&
                    symm.equal(a.e31, b.e31) && symm.equal(a.e23, b.e23) && symm.equal(a.e33, b.e33));
            };
        std::vector<int> isymbvk2isym(bvkgmat.size(), -1);
        for (int isymbvk = 0;isymbvk < bvkgmat.size();++isymbvk)
        {
            for (int isym = 0;isym < symm.nrotk;++isym)
            {
                if (matequal(bvkgmat[isymbvk], symm.gmatrix[isym]))
                {
                    isymbvk2isym[isymbvk] = isym;
                    break;
                }
            }
        }
        return isymbvk2isym;
    }

    inline int gcd(const int a, const int b)
    {
        assert(a > 0 && b > 0);
        int c = a % b;
        return (c == 0) ? b : gcd(b, c);
    }
    void Irreducible_Sector::gen_symmetry_BvK(const ModuleSymmetry::Symmetry& symm, const Atom* atoms, const Lattice& lat, const Statistics& st, const TC bvk_period)
    {
        ModuleBase::TITLE("Irreducible_Sector", "gen_symmetry_BvK");
        auto set_matrix3 = [](const ModuleBase::Vector3<double>& a1, const ModuleBase::Vector3<double>& a2, const ModuleBase::Vector3<double>& a3)
            -> ModuleBase::Matrix3 {return ModuleBase::Matrix3(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z, a3.x, a3.y, a3.z);};

        if (bvk_period[0] == bvk_period[1] && bvk_period[0] == bvk_period[2])
        {   //the BvK supercell has the same symmetry as the original cell
            this->bvk_nsym_ = symm.nrotk;
            this->isymbvk_to_isym_.resize(symm.nrotk);
            for (int isym = 0;isym < symm.nrotk;++isym) { this->isymbvk_to_isym_[isym] = isym;
}
            return;
        }

        // extern lattice to minimal BvK lattice, and set direct coordinates in min BvK lattice 
        int bvk_gcd = gcd(bvk_period[0], gcd(bvk_period[1], bvk_period[2]));
        const TC bvk_min_period = TC({ bvk_period[0] / bvk_gcd, bvk_period[1] / bvk_gcd, bvk_period[2] / bvk_gcd });
        const int bvk_nat = st.nat * bvk_min_period[0] * bvk_min_period[1] * bvk_min_period[2];
        std::vector<int> bvk_na(st.ntype);
        std::vector<int> bvk_istart(st.ntype, 0);
        int bvk_itmin_start = 0, bvk_itmin_type = 0;
        for (int it = 0;it < st.ntype;++it)
        {
            bvk_na[it] = atoms[it].na * bvk_min_period[0] * bvk_min_period[1] * bvk_min_period[2];
            if (it > 0) { bvk_istart[it] = bvk_istart[it - 1] + bvk_na[it - 1];
}
            if (bvk_na[it] < bvk_na[bvk_itmin_type])
            {
                bvk_itmin_type = it;
                bvk_itmin_start = bvk_istart[it];
            }
        }

        std::vector<double> bvk_dpos(3 * bvk_nat);
        std::vector<double> bvk_rot_dpos(3 * bvk_nat);
        std::vector<int> order_index(bvk_nat + 2);
        ModuleBase::Vector3<double> a1, a2, a3, s1, s2, s3; // a: to be optimized; s: original
        s1 = a1 = lat.a1 * static_cast<double>(bvk_min_period[0]);
        s2 = a2 = lat.a2 * static_cast<double>(bvk_min_period[1]);
        s3 = a3 = lat.a3 * static_cast<double>(bvk_min_period[2]);
        ModuleBase::Matrix3 bvk_min_lat = set_matrix3(s1, s2, s3);
        int at = 0;
        for (int it = 0; it < st.ntype; ++it) {
            for (int c1 = 0;c1 < bvk_min_period[0];++c1) {
                for (int c2 = 0;c2 < bvk_min_period[1];++c2) {
                    for (int c3 = 0;c3 < bvk_min_period[2];++c3) {
                        for (int ia = 0; ia < atoms[it].na; ++ia)
                        {
                            bvk_dpos[3 * at] = (static_cast<double> (c1) + atoms[it].taud[ia].x) / static_cast<double>(bvk_min_period[0]);
                            bvk_dpos[3 * at + 1] = (static_cast<double> (c2) + atoms[it].taud[ia].y) / static_cast<double>(bvk_min_period[1]);
                            bvk_dpos[3 * at + 2] = (static_cast<double> (c3) + atoms[it].taud[ia].z) / static_cast<double>(bvk_min_period[2]);
                            for (int k = 0; k < 3; ++k)
                            {
                                symm.check_translation(bvk_dpos[3 * at + k], -floor(bvk_dpos[3 * at + k]));
                                symm.check_boundary(bvk_dpos[3 * at + k]);
                            }
                            ++at;
                        }
}
}
}
}

        // analyze bravis and generate optimized lattice for minimal BvK lattice
        double cel_const[6];
        double pre_const[6];
        int bvk_brav = 0;
        std::string bvk_latname="";
        // bvk_brav = symm.standard_lat(s1, s2, s3, cel_const); //not enough, optimal lattice may change after cell-extension
        symm.lattice_type(a1, a2, a3, s1, s2, s3, cel_const, pre_const, bvk_brav, bvk_latname, nullptr, false, nullptr);
        ModuleBase::Matrix3 bvk_min_optlat = set_matrix3(a1, a2, a3);
        // convert the direct coordinates to the optimized lattice
        for (int i = 0;i < bvk_nat;++i)
        {
            ModuleBase::Vector3<double> taud(bvk_dpos[3 * i], bvk_dpos[3 * i + 1], bvk_dpos[3 * i + 2]);
            taud = taud * bvk_min_lat * bvk_min_optlat.Inverse();
            bvk_dpos[3 * i] = taud.x;
            bvk_dpos[3 * i + 1] = taud.y;
            bvk_dpos[3 * i + 2] = taud.z;
            for (int k = 0; k < 3; ++k)
            {
                symm.check_translation(bvk_dpos[3 * i + k], -floor(bvk_dpos[3 * i + k]));
                symm.check_boundary(bvk_dpos[3 * i + k]);
            }
        }

        // generate symmetry operation of the BvK lattice using the original optlat-direct coordinates
        std::vector<ModuleBase::Matrix3> bvk_op(48);
        int bvk_nop;
        symm.setgroup(bvk_op.data(), bvk_nop, bvk_brav);
        bvk_op.resize(bvk_nop);
        int bvk_npg, bvk_nsg, bvk_pgnum, bvk_sgnum;
        std::string bvk_pgname, bvk_sgname;
        this->bvk_gmatrix_.resize(48);
        this->bvk_gtrans_.resize(48);
        symm.getgroup(bvk_npg, bvk_nsg, GlobalV::ofs_running, bvk_nop,
            bvk_op.data(), this->bvk_gmatrix_.data(), this->bvk_gtrans_.data(),
            bvk_dpos.data(), bvk_rot_dpos.data(), order_index.data(),
            bvk_itmin_type, bvk_itmin_start, bvk_istart.data(), bvk_na.data());
        this->bvk_gmatrix_.resize(bvk_nsg);
        this->bvk_gtrans_.resize(bvk_nsg);
        this->bvk_nsym_ = bvk_nsg;
        symm.pointgroup(bvk_npg, bvk_pgnum, bvk_pgname, this->bvk_gmatrix_.data(), GlobalV::ofs_running);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "POINT GROUP OF BvK SCELL", bvk_pgname);
        symm.pointgroup(bvk_nsg, bvk_sgnum, bvk_sgname, this->bvk_gmatrix_.data(), GlobalV::ofs_running);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "POINT GROUP IN SPACE GROUP OF BvK SCELL", bvk_sgname);
        symm.gmatrix_convert_int(this->bvk_gmatrix_.data(), this->bvk_gmatrix_.data(), bvk_nsg, bvk_min_optlat, lat.latvec);
        symm.gtrans_convert(this->bvk_gtrans_.data(), this->bvk_gtrans_.data(), bvk_nsg, bvk_min_optlat, lat.latvec);
        // get map from bvk-op to original op
        this->isymbvk_to_isym_ = get_isymbvk_to_isym_map(this->bvk_gmatrix_, symm);
        return;
    }

    // std::vector<bool> Irreducible_Sector::in_plain(const ModuleSymmetry::Symmetry& symm, const ModuleBase::Matrix3& latvec)const
    // {
    //     // get euler angel of the cartesian gmatrix in optimal lattice
    //     std::vector<ModuleBase::Matrix3> gmatc(symm.nrotk);
    //     symm.gmatrix_convert_int(symm.gmatrix, gmatc.data(), symm.nrotk, latvec, symm.optlat);
    //     for (auto& g : gmatc) g = direct_to_cartesian(g, symm.optlat);
    //     std::vector<bool> in_plain(symm.nrotk, false);
    //     for (int i = 0;i < symm.nrotk;++i)
    //     {
    //         TCdouble euler_angle = get_euler_angle(gmatc[i]);
    //         if (symm.equal(euler_angle.y, 0.0) || symm.equal(euler_angle.y, ModuleBase::PI)) in_plain[i] = true;
    //     }
    //     return in_plain;
    // }

};
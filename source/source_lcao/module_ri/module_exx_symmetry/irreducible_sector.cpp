#include "source_lcao/module_ri/module_exx_symmetry/irreducible_sector.h"
#include "source_io/module_parameter/parameter.h"
namespace ModuleSymmetry
{
    // Raw-index dispatch shared by the real-space sector helpers, matching the convention used
    // everywhere else (symmetry_rotation.h): isym < nrotk -> unitary gmatrix[isym];
    // isym >= nrotk -> spatial part of the antiunitary element Theta*gmatrix_anti[isym-nrotk].
    // Only the SPATIAL part is needed here: Theta acts on H(R) as sigma_y (.)^* sigma_y and
    // leaves R and the atom pair untouched, so the sector bookkeeping is identical for both kinds.
    static inline const ModuleBase::Matrix3& sector_gmatrix(const Symmetry& symm, const int isym)
    {
        return (isym < symm.nrotk) ? symm.gmatrix[isym] : symm.gmatrix_anti[isym - symm.nrotk];
    }
    static inline int sector_rotated_atom(const Symmetry& symm, const int isym, const int iat)
    {
        return (isym < symm.nrotk) ? symm.get_rotated_atom(isym, iat)
                                   : symm.get_rotated_atom_anti(isym - symm.nrotk, iat);
    }

    TC Irreducible_Sector::rotate_R(const Symmetry& symm,
        const int isym, const int iat1, const int iat2, const TC& R, const char gauge) const
    {
        auto round2int = [symm](const double x) -> int { return x > 0 ? static_cast<int>(x + symm.epsilon) : static_cast<int>(x - symm.epsilon); };
        const TCdouble R_double(static_cast<double>(R[0]), static_cast<double>(R[1]), static_cast<double>(R[2]));
        // return_lattice_ is already sized nrotk+nrotk_anti and indexed by the same raw isym.
        const ModuleBase::Matrix3& gmat = sector_gmatrix(symm, isym);
        const TCdouble Rrot_double = (gauge == 'L')
            ? R_double * gmat + this->return_lattice_[iat1][isym] - this->return_lattice_[iat2][isym]
            : R_double * gmat + this->return_lattice_[iat2][isym] - this->return_lattice_[iat1][isym];
        return { round2int(Rrot_double.x), round2int(Rrot_double.y), round2int(Rrot_double.z) };
    }
    TapR Irreducible_Sector::rotate_apR_by_formula(const Symmetry& symm,
        const int isym, const TapR& apR, const char gauge) const
    {
        const Tap& aprot = { sector_rotated_atom(symm, isym, apR.first.first),
                             sector_rotated_atom(symm, isym, apR.first.second) };
        return { aprot, this->rotate_R(symm, isym, apR.first.first, apR.first.second, apR.second, gauge) };
    }

    TCdouble Irreducible_Sector::get_aRb_direct(const Atom* atoms, const Statistics& st,
        const int iat1, const int iat2, const TCdouble& R, const char gauge)const
    {
        return TCdouble(atoms[st.iat2it[iat1]].taud[st.iat2ia[iat1]] - atoms[st.iat2it[iat2]].taud[st.iat2ia[iat2]]) + (gauge == 'L' ? R : TCdouble(-R));
    }

    TCdouble Irreducible_Sector::get_aRb_direct(const Atom* atoms, const Statistics& st,
        const int iat1, const int iat2, const TC& R, const char gauge) const
    {
        const TCdouble R_double(static_cast<double>(R[0]), static_cast<double>(R[1]), static_cast<double>(R[2]));
        return get_aRb_direct(atoms, st, iat1, iat2, R_double);
    }

    inline void output_return_lattice(const std::vector<std::vector<TCdouble>>& return_lattice)
    {
        std::cout << "return lattice:" << std::endl;
        for (int iat = 0;iat < return_lattice.size();++iat)
        {
            std::cout << "atom" << iat << std::endl;
            for (int isym = 0;isym < return_lattice[iat].size();++isym)
                std::cout << "isym=" << isym << ", return lattice=" <<
                return_lattice[iat][isym].x << " " << return_lattice[iat][isym].y << " " << return_lattice[iat][isym].z << std::endl;
        }
    }

    // Perfoming {R|t} to atom position r in the R=0 lattice, we get Rr+t, which may get out of R=0 lattice, 
    // whose image in R=0 lattice is r'=Rr+t-O. This function is to get O for each atom and each symmetry operation.
    // the range of direct position is [-0.5, 0.5).
    TCdouble Irreducible_Sector::get_return_lattice(const Symmetry& symm,
        const ModuleBase::Matrix3& gmatd, const TCdouble gtransd,
        const TCdouble& posd_a1, const TCdouble& posd_a2)const
    {
        // auto restrict_center = [&symm](const TCdouble& v) -> TCdouble {
        //     // in [-0.5, 0.5)
        //     TCdouble vr;
        //     vr.x = fmod(v.x + 100.5 + 0.5 * symm.epsilon, 1) - 0.5 - 0.5 * symm.epsilon;
        //     vr.y = fmod(v.y + 100.5 + 0.5 * symm.epsilon, 1) - 0.5 - 0.5 * symm.epsilon;
        //     vr.z = fmod(v.z + 100.5 + 0.5 * symm.epsilon, 1) - 0.5 - 0.5 * symm.epsilon;
        //     if (std::abs(vr.x) < symm.epsilon) vr.x = 0.0;
        //     if (std::abs(vr.y) < symm.epsilon) vr.y = 0.0;
        //     if (std::abs(vr.z) < symm.epsilon) vr.z = 0.0;
        //     return vr;
        //     };
        auto restrict_center = [&symm](const TCdouble& v) -> TCdouble {
            // in [0,1)
            TCdouble vr;
            vr.x = fmod(v.x + 100 + symm.epsilon, 1) - symm.epsilon;
            vr.y = fmod(v.y + 100 + symm.epsilon, 1) - symm.epsilon;
            vr.z = fmod(v.z + 100 + symm.epsilon, 1) - symm.epsilon;
            if (std::abs(vr.x) < symm.epsilon) vr.x = 0.0;
            if (std::abs(vr.y) < symm.epsilon) vr.y = 0.0;
            if (std::abs(vr.z) < symm.epsilon) vr.z = 0.0;
            return vr;
            };
        auto check_integer = [&symm](const double x) -> void {
            assert(symm.equal(x, std::round(x)));
            };
        TCdouble rotpos1 = restrict_center(posd_a1) * gmatd + restrict_center(gtransd);  // row vector
        TCdouble return_lattice_double = rotpos1 - restrict_center(posd_a2);
#ifdef __DEBUG
        check_integer(return_lattice_double.x);
        check_integer(return_lattice_double.y);
        check_integer(return_lattice_double.z);
#endif
        return TCdouble(std::round(return_lattice_double.x), std::round(return_lattice_double.y), std::round(return_lattice_double.z));
    }

    void Irreducible_Sector::cal_return_lattice_all(const Symmetry& symm, const Atom* atoms, const Statistics& st)
    {
        ModuleBase::TITLE("Symmetry_rotation", "cal_return_lattice_all");
        // Columns [0, nrotk) are the unitary operations; columns [nrotk, nrotk+nrotk_anti) are the
        // spatial parts of the antiunitary elements Theta*g of the Shubnikov group (nspin=4 magnetic),
        // so that Symmetry_rotation can address both with one raw index. 
        this->return_lattice_.resize(st.nat, std::vector<TCdouble>(symm.nrotk + symm.nrotk_anti));
        for (int iat1 = 0;iat1 < st.nat;++iat1)
        {
            int it = st.iat2it[iat1];
            int ia1 = st.iat2ia[iat1];
            for (int isym = 0;isym < symm.nrotk;++isym)
            {
                int iat2 = symm.get_rotated_atom(isym, iat1);
                int ia2 = st.iat2ia[iat2];
                this->return_lattice_[iat1][isym] = get_return_lattice(symm, symm.gmatrix[isym], symm.gtrans[isym], atoms[it].taud[ia1], atoms[it].taud[ia2]);
            }
            for (int j = 0;j < symm.nrotk_anti;++j)
            {
                int iat2 = symm.get_rotated_atom_anti(j, iat1);
                int ia2 = st.iat2ia[iat2];
                this->return_lattice_[iat1][symm.nrotk + j] = get_return_lattice(symm, symm.gmatrix_anti[j], symm.gtrans_anti[j], atoms[it].taud[ia1], atoms[it].taud[ia2]);
            }
        }
        // test: output return_lattice
        // output_return_lattice(this->return_lattice_);
    }

    void Irreducible_Sector::output_full_map_to_irreducible_sector(const int nat)
    {
        std::cout << "Final map to irreducible sector: " << std::endl;
        for (auto& apR_isym_irapR : this->full_map_to_irreducible_sector_)
        {
            const Tap& ap = apR_isym_irapR.first.first;
            const TC& R = apR_isym_irapR.first.second;
            const Tap& irap = apR_isym_irapR.second.second.first;
            const TC& irR = apR_isym_irapR.second.second.second;
            std::cout << "atompair (" << ap.first << ", " << ap.second << "), R=(" << R[0] << ", " << R[1] << ", " << R[2] << ") -> "
                << "isym=" << apR_isym_irapR.second.first << " -> irreducible atompair (" << irap.first << ", " << irap.second << "), irreducible R=("
                << irR[0] << ", " << irR[1] << ", " << irR[2] << ")" << std::endl;
        }
    }

    void Irreducible_Sector::output_sector_star()
    {
        std::cout << "Found " << this->sector_stars_.size() << " irreducible sector stars:" << std::endl;
        // for (auto& irs_star : this->sector_stars_)
        for (auto& irap_star : this->sector_stars_)
        {
            const TapR& irapR = irap_star.first;
            const auto& star = irap_star.second;
            const Tap& irap = irapR.first;
            const TC& irR = irapR.second;
            std::cout << "in star of irreducible atompair=(" << irap.first << ", " << irap.second << "), R=(" << irR[0] << ", " << irR[1] << ", " << irR[2] << ") with size " << star.size() << ":\n";
            for (auto& isym_ap_R : star)
                std::cout << "isym=" << isym_ap_R.first << ", atompair=(" << isym_ap_R.second.first.first << ", " << isym_ap_R.second.first.second << "), R=("
                << isym_ap_R.second.second[0] << ", " << isym_ap_R.second.second[1] << ", " << isym_ap_R.second.second[2] << ")" << std::endl;
        }
        // print irreducible sector
        std::cout << "irreducible sector: " << std::endl;
        for (auto& irap_irR : this->irreducible_sector_)
        {
            for (auto& irR : irap_irR.second) {std::cout << "atompair (" << irap_irR.first.first << ", " << irap_irR.first.second << "), R = (" << irR[0] << ", " << irR[1] << ", " << irR[2] << ") \n";}
            std::cout << std::endl;
        }
    }
    void Irreducible_Sector::write_irreducible_sector()
    {
        if(GlobalV::MY_RANK == 0)
        {
            std::ofstream ofs;
            ofs.open(PARAM.globalv.global_out_dir + "irreducible_sector.txt");
            for (auto& irap_irR : this->irreducible_sector_)
            {
                for (auto& irR : irap_irR.second){ofs << "atompair (" << irap_irR.first.first << ", " << irap_irR.first.second << "), R = (" << irR[0] << ", " << irR[1] << ", " << irR[2] << ") \n";}
            }
            ofs.close();
        }
    }

    void Irreducible_Sector::find_irreducible_sector(const Symmetry& symm, const Atom* atoms, const Statistics& st, const std::vector<TC>& Rs, const TC& period, const Lattice& lat)
    {
        this->full_map_to_irreducible_sector_.clear();
        this->irreducible_sector_.clear();
        this->sector_stars_.clear();

        if (this->return_lattice_.empty()) this->cal_return_lattice_all(symm, atoms, st);
        // if (this->atompair_stars_.empty()) this->find_irreducible_atom_pairs(symm);

        // contruct {atom pair, R} set
        // constider different number of Rs for different atom pairs later.
        std::map<Tap, std::set<TC, len_less_func>> apR_all;
        for (int iat1 = 0;iat1 < st.nat; iat1++)
            for (int iat2 = 0; iat2 < st.nat; iat2++)
                for (auto& R : Rs)
                    apR_all[{iat1, iat2}].insert(R);

        // get invmap over the operation set actually used by the sector search.
        // For nspin=4 magnetic that is the full Shubnikov group H (union) A, laid out as
        // [gmatrix[0..nrotk) | gmatrix_anti[0..nrotk_anti)]. gmatrix_invmap needs no change:
        // it searches the whole array for s[i]*s[j]==I, and A is closed under inversion
        // (if g in A had g^-1 in H then g = (g^-1)^-1 would be in H, contradicting H n A = {}),
        // so the concatenated array is exactly the parent group and every inverse is found,
        // with the inverse of a coset member landing inside the coset block.
        if (this->invmap_.empty())
        {
            const int nop = symm.nrotk + symm.nrotk_anti;
            std::vector<ModuleBase::Matrix3> gmat_all(nop);
            for (int i = 0; i < symm.nrotk; ++i) { gmat_all[i] = symm.gmatrix[i]; }
            for (int j = 0; j < symm.nrotk_anti; ++j) { gmat_all[symm.nrotk + j] = symm.gmatrix_anti[j]; }
            this->invmap_.resize(nop);
            symm.gmatrix_invmap(gmat_all.data(), nop, invmap_.data());
        }

        // get symmetry of BvK supercell
        if (this->isymbvk_to_isym_.empty())
            this->gen_symmetry_BvK(symm, atoms, lat, st, period);
        assert(!this->isymbvk_to_isym_.empty());
        // std::vector<bool> in_2d_plain;
        // const bool judge_2d = (symm.real_brav == 4);
        // if (judge_2d) in_2d_plain = this->in_plain(symm, lat.latvec);

        while (!apR_all.empty())
        {
            const Tap irap = apR_all.begin()->first;
            const TC irR = *apR_all[irap].begin();
            const TapR& irapR = { irap, irR };
            std::map<int, TapR> sector_star;
            for (int isymbvk = 0;isymbvk < this->bvk_nsym_;++isymbvk)
            {
                // if (judge_2d)
                // {
                //     std::cout << "isym=" << isym << ", inv[isym]=" << invmap_[isym] << std::endl;
                //     assert(in_2d_plain[isym] == in_2d_plain[invmap_[isym]]);
                //     if (!in_2d_plain[isym]) continue;
                // }
                const int& isym = this->isymbvk_to_isym_[isymbvk];
                // A BvK operation with no counterpart in the unit-cell operation set is marked -1.
                // For a magnetic nspin=4 system the unit-cell set is the Shubnikov group H (union) A,
                // which is generally a PROPER subset of the crystallographic group (operations that
                // merely tilt the moment belong to neither), so an unmatched BvK operation is a
                // normal outcome, not an error: it is simply not a symmetry of the magnetic system.
                // Skipping it only costs reduction, never correctness.
                if (isym < 0) { continue; }
                const TapR& apRrot = this->rotate_apR_by_formula(symm, this->invmap_[isym], irapR);
                const Tap& aprot = apRrot.first;
                const TC& Rrot = apRrot.second;
                if (apR_all.count(aprot) && apR_all.at(aprot).count(Rrot))
                {
                    this->full_map_to_irreducible_sector_.insert({ apRrot, {isym, irapR} });
                    sector_star.insert({ isym, apRrot });
                    apR_all[aprot].erase(Rrot);
                    if (apR_all.at(aprot).empty()) apR_all.erase(aprot);
                    if (apR_all.empty()) break;
                }
            }// end for isym
            if (!sector_star.empty())
            {
                const TapR& irapR = { irap, irR };
                this->sector_stars_.emplace(irapR, sector_star);
                if (this->irreducible_sector_.count(irap))
                    this->irreducible_sector_.at(irap).insert(irR);
                else
                    this->irreducible_sector_.insert({ irap, {irR} });
            }
        }
        // test
        int total_apR_in_star = 0;
        for (auto& sector : this->sector_stars_)
            total_apR_in_star += sector.second.size();
        assert(total_apR_in_star == this->full_map_to_irreducible_sector_.size());
        // this->output_full_map_to_irreducible_sector(st.nat);
        // this->output_sector_star();
        this->write_irreducible_sector();
    }
}

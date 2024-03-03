#include "symmetry_rotation.h"

namespace ModuleSymmetry
{
    int Symmetry_rotation::round2int(const double x) const
    {
        return x > 0 ? static_cast<int>(x + eps_) : static_cast<int>(x - eps_);
    }

    TC Symmetry_rotation::rotate_R_by_formula(const Symmetry& symm,
        const int isym, const int iat1, const int iat2, const TC& R, const char gauge) const
    {
        const TCdouble R_double(static_cast<double>(R[0]), static_cast<double>(R[1]), static_cast<double>(R[2]));
        const TCdouble Rrot_double = (gauge == 'L')
            ? R_double * symm.gmatrix[isym] + this->return_lattice_[iat1][isym] - this->return_lattice_[iat2][isym]
            : R_double * symm.gmatrix[isym] + this->return_lattice_[iat2][isym] - this->return_lattice_[iat1][isym];
        return { round2int(Rrot_double.x), round2int(Rrot_double.y), round2int(Rrot_double.z) };
    }
    TapR Symmetry_rotation::rotate_R_by_formula(const Symmetry& symm,
        const int isym, const TapR& apR, const char gauge) const
    {
        const Tap& aprot = { symm.get_rotated_atom(isym, apR.first.first), symm.get_rotated_atom(isym, apR.first.second) };
        return { aprot, this->rotate_R_by_formula(symm, isym, apR.first.first, apR.first.second, apR.second, gauge) };
    }

    TCdouble Symmetry_rotation::get_aRb_direct(const Atom* atoms, const Statistics& st,
        const int iat1, const int iat2, const TCdouble& R, const char gauge)const
    {
        return TCdouble(atoms[st.iat2it[iat1]].taud[st.iat2ia[iat1]] - atoms[st.iat2it[iat2]].taud[st.iat2ia[iat2]]) + (gauge == 'L' ? R : TCdouble(-R));
    }

    TCdouble Symmetry_rotation::get_aRb_direct(const Atom* atoms, const Statistics& st,
        const int iat1, const int iat2, const TC& R, const char gauge) const
    {
        const TCdouble R_double(static_cast<double>(R[0]), static_cast<double>(R[1]), static_cast<double>(R[2]));
        return get_aRb_direct(atoms, st, iat1, iat2, R_double);
    }

    // void Symmetry_rotation::find_irreducible_atom_pairs(const Symmetry& symm)
    // {
    //     ModuleBase::TITLE("Symmetry_rotation", "find_irreducible_atom_pairs");
    //     this->eps_ = symm.epsilon;
    //     for (int iat1 = 0;iat1 < symm.nat;++iat1)
    //         for (int iat2 = 0;iat2 < symm.nat;++iat2)
    //         {
    //             Tap pair = { iat1, iat2 };
    //             bool exist = false;
    //             for (int isym = 0;isym < symm.nrotk;++isym)
    //             {
    //                 Tap rotpair = { symm.get_rotated_atom(isym,iat1), symm.get_rotated_atom(isym,iat2) };
    //                 for (int ip = 0;ip < this->atompair_stars_.size();++ip) // current irreduceble pairs
    //                 {
    //                     if (rotpair == this->atompair_stars_[ip].at(0))
    //                     {
    //                         this->atompair_stars_[ip].insert({ isym, pair });
    //                         exist = true;
    //                         break;
    //                     }
    //                 }
    //                 if (exist)break;
    //             }
    //             if (!exist)this->atompair_stars_.push_back({ {0, pair} });
    //         }
    // }

    // void Symmetry_rotation::find_irreducible_atom_pairs_set(const Symmetry& symm)
    // {
    //     ModuleBase::TITLE("Symmetry_rotation", "find_irreducible_atom_pairs_set");
    //     this->eps_ = symm.epsilon;
    //     if (this->invmap_.empty())
    //     {
    //         this->invmap_.resize(symm.nrotk);
    //         symm.gmatrix_invmap(symm.gmatrix, symm.nrotk, invmap_.data());
    //     }
    //     // contruct initial ap-set
    //     std::set<Tap, ap_less_func> ap_set;
    //     for (int iat1 = 0; iat1 < symm.nat; ++iat1)
    //         for (int iat2 = 0; iat2 < symm.nat; ++iat2)
    //             ap_set.insert({ iat1, iat2 });
    //     while (!ap_set.empty())
    //     {
    //         Tap ap = *ap_set.begin();
    //         std::map<int, Tap> ap_star;
    //         for (int isym = 0; isym < symm.nrotk; ++isym)
    //         {
    //             Tap rotpair = { symm.get_rotated_atom(isym,ap.first), symm.get_rotated_atom(isym,ap.second) };
    //             if (ap_set.find(rotpair) != ap_set.end())
    //             {
    //                 ap_star.insert({ this->invmap_[isym], rotpair });
    //                 ap_set.erase(rotpair);
    //             }
    //         }
    //         this->atompair_stars_.push_back(ap_star);
    //     }
    // }

    // inline void output_atompair_stars(const std::vector<std::map<int, Tap>>& ap_stars)
    // {
    //     std::cout << "stars of irreducible atom pairs: " << std::endl;
    //     for (int ip = 0; ip < ap_stars.size(); ++ip)
    //     {
    //         std::cout << "atom pair star " << ip << ": " << std::endl;
    //         for (const auto& ap : ap_stars[ip])
    //             std::cout << "isym=" << ap.first << ", atompair=(" << ap.second.first << ", " << ap.second.second << ") " << std::endl;
    //     }
    // }

    // void Symmetry_rotation::test_irreducible_atom_pairs(const Symmetry& symm)
    // {
    //     std::cout << "Algorithm 1: find irreducible atom pairs by set" << std::endl;
    //     this->find_irreducible_atom_pairs_set(symm);
    //     output_atompair_stars(this->atompair_stars_);

    //     std::cout << std::endl;
    //     this->atompair_stars_.clear();
    //     std::cout << "Algorithm 2: find irreducible atom pairs by judge" << std::endl;
    //     this->find_irreducible_atom_pairs(symm);
    //     output_atompair_stars(this->atompair_stars_);
    // }

    std::vector<TC> Symmetry_rotation::get_Rs_from_BvK(const K_Vectors& kv) const
    {
        const TC& period = RI_Util::get_Born_vonKarmen_period(kv);
        return RI_Util::get_Born_von_Karmen_cells(period);
    }
    std::vector<TC> Symmetry_rotation::get_Rs_from_adjacent_list(const UnitCell& ucell, Grid_Driver& gd, const Parallel_Orbitals& pv) const
    {
        // find the union set of Rs for all the atom pairs
        std::set<TC> Rs_set;
        for (int iat1 = 0;iat1 < ucell.nat;++iat1)
        {
            auto tau1 = ucell.get_tau(iat1);
            int it1 = ucell.iat2it[iat1], ia1 = ucell.iat2ia[iat1];
            AdjacentAtomInfo adjs;
            gd.Find_atom(ucell, tau1, it1, ia1, &adjs);
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int it2 = adjs.ntype[ad];
                const int ia2 = adjs.natom[ad];
                int iat2 = ucell.itia2iat(it2, ia2);
                if (pv.get_row_size(iat1) && pv.get_col_size(iat2))
                {
                    const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
                    if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0
                        < ucell.atoms[it1].Rcut + ucell.atoms[it2].Rcut)
                        Rs_set.insert({ R_index.x, R_index.y, R_index.z });
                }
            }
        }
        // set to vector
        std::vector<TC> Rs(Rs_set.size());
        for (auto& R : Rs_set) Rs.push_back(R);
        return Rs;
    }

    void Symmetry_rotation::output_full_map_to_irreducible_sector(const int nat)
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

    void Symmetry_rotation::output_sector_star()
    {
        std::cout << "Found " << this->sector_stars_.size() << " irreducible sector stars:" << std::endl;
        // for (auto& irs_star : this->sector_stars_)
        for (int istar = 0;istar < this->sector_stars_.size();++istar)
        {
            std::cout << "in star " << istar << " with size " << this->sector_stars_[istar].size() << ":\n";
            for (auto& isym_ap_R : this->sector_stars_[istar])
                std::cout << "isym=" << isym_ap_R.first << ", atompair=(" << isym_ap_R.second.first.first << ", " << isym_ap_R.second.first.second << "), R=("
                << isym_ap_R.second.second[0] << ", " << isym_ap_R.second.second[1] << ", " << isym_ap_R.second.second[2] << ")" << std::endl;
        }
        // print irreducible sector
        std::cout << "irreducible sector: " << std::endl;
        for (auto& irap_irR : this->irreducible_sector_)
        {
            for (auto& irR : irap_irR.second)
                std::cout << "atompair (" << irap_irR.first.first << ", " << irap_irR.first.second << "), R = (" << irR[0] << ", " << irR[1] << ", " << irR[2] << ") \n";
            std::cout << std::endl;
        }
    }

    void Symmetry_rotation::find_irreducible_sector(const Symmetry& symm, const Atom* atoms, const Statistics& st, const std::vector<TC>& Rs)
    {
        this->full_map_to_irreducible_sector_.clear();
        this->irreducible_sector_.clear();
        this->sector_stars_.clear();

        if (this->return_lattice_.empty()) this->get_return_lattice_all(symm, atoms, st);
        // if (this->atompair_stars_.empty()) this->find_irreducible_atom_pairs(symm);

        // contruct {atom pair, R} set
        // constider different number of Rs for different atom pairs later.
        std::map<Tap, std::set<TC, len_less_func>> apR_all;
        for (int iat1 = 0;iat1 < st.nat; iat1++)
            for (int iat2 = 0; iat2 < st.nat; iat2++)
                for (auto& R : Rs)
                    apR_all[{iat1, iat2}].insert(R);

        // get invmap
        if (this->invmap_.empty())
        {
            this->invmap_.resize(symm.nrotk);
            symm.gmatrix_invmap(symm.gmatrix, symm.nrotk, invmap_.data());
        }

        while (!apR_all.empty())
        {
            const Tap irap = apR_all.begin()->first;
            const TC irR = *apR_all[irap].begin();
            const TapR& irapR = { irap, irR };
            std::map<int, TapR> sector_star;
            for (int isym = 0;isym < symm.nrotk;++isym)
            {
                const TapR& apRrot = this->rotate_R_by_formula(symm, this->invmap_[isym], irapR);
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
                this->sector_stars_.push_back(sector_star);
                if (this->irreducible_sector_.count(irap))
                    this->irreducible_sector_.at(irap).insert(irR);
                else
                    this->irreducible_sector_.insert({ irap, {irR} });
            }
        }
        // test
        int total_apR_in_star = 0;
        for (auto& sector : this->sector_stars_)
            total_apR_in_star += sector.size();
        assert(total_apR_in_star == this->full_map_to_irreducible_sector_.size());
        this->output_full_map_to_irreducible_sector(st.nat);
        this->output_sector_star();
    }
}
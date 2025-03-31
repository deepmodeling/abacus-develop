#include "symmetry.h"
using namespace ModuleSymmetry;

#include <set>

void Symmetry::analyze_magnetic_group(const Atom* atoms, const Statistics& st, int& nrot_out, int& nrotk_out)
{
    // 1. classify atoms with different magmom
    //  (use symmetry_prec to judge if two magmoms are the same)
    std::vector<std::set<int>> mag_type_atoms;
    for (int it = 0;it < ntype;++it)
    {
        for (int ia = 0; ia < atoms[it].na; ++ia)
        {
            bool find = false;
            for (auto& mt : mag_type_atoms)
            {
                const int mag_iat = *mt.begin();
                const int mag_it = st.iat2it[mag_iat];
                const int mag_ia = st.iat2ia[mag_iat];
                if (it == mag_it && this->equal(atoms[it].mag[ia], atoms[mag_it].mag[mag_ia]))
                {
                    mt.insert(st.itia2iat(it, ia));
                    find = true;
                    break;
                }
            }
            if (!find)
            {
                mag_type_atoms.push_back(std::set<int>({ st.itia2iat(it,ia) }));
            }
        }
    }

    // 2. get the start index, number of atoms and positions for each mag_type
    std::vector<int> mag_istart(mag_type_atoms.size());
    std::vector<int> mag_na(mag_type_atoms.size());
    std::vector<double> mag_pos;
	int mag_itmin_type = 0;
	int mag_itmin_start = 0;
	for (int mag_it = 0;mag_it < mag_type_atoms.size(); ++mag_it)
	{
		mag_na[mag_it] = mag_type_atoms.at(mag_it).size();
		if (mag_it > 0)
		{
			mag_istart[mag_it] = mag_istart[mag_it - 1] + mag_na[mag_it - 1];
		}
		if (mag_na[mag_it] < mag_na[itmin_type])
		{
			mag_itmin_type = mag_it;
			mag_itmin_start = mag_istart[mag_it];
		}
		for (auto& mag_iat : mag_type_atoms.at(mag_it))
		{
			// this->newpos have been ordered by original structure(ntype, na), it cannot be directly used here.
			// we need to reset the calculate again the coordinate of the new structure.
			const ModuleBase::Vector3<double> direct_tmp = atoms[st.iat2it[mag_iat]].tau[st.iat2ia[mag_iat]] * this->optlat.Inverse();
			std::array<double, 3> direct = { direct_tmp.x, direct_tmp.y, direct_tmp.z };
			for (int i = 0; i < 3; ++i)
			{
				this->check_translation(direct[i], -floor(direct[i]));
				this->check_boundary(direct[i]);
				mag_pos.push_back(direct[i]);
			}
		}
	}

	// 3. analyze the effective structure
	this->getgroup(nrot_out, nrotk_out, GlobalV::ofs_running, 
			this->nop, this->symop, this->gmatrix, 
			this->gtrans, mag_pos.data(), this->rotpos, 
			this->index, mag_type_atoms.size(), mag_itmin_type, 
			mag_itmin_start, mag_istart.data(), mag_na.data());

}

bool Symmetry::magmom_same_check(const Atom* atoms)const
{
    ModuleBase::TITLE("Symmetry", "magmom_same_check");
    bool pricell_loop = true;
    for (int it = 0;it < ntype;++it)
    {
        if (pricell_loop) {
            for (int ia = 1;ia < atoms[it].na;++ia)
            {
                if (!equal(atoms[it].m_loc_[ia].x, atoms[it].m_loc_[0].x) ||
                    !equal(atoms[it].m_loc_[ia].y, atoms[it].m_loc_[0].y) ||
                    !equal(atoms[it].m_loc_[ia].z, atoms[it].m_loc_[0].z))
                {
                    pricell_loop = false;
                    break;
                }
            }
        }
    }
    return pricell_loop;
}


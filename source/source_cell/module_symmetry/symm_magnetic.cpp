#include "symmetry.h"
using namespace ModuleSymmetry;

#include "symmetry_rotation_spin.h"
#include "source_io/module_parameter/parameter.h"

#include <set>
#include <vector>

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

void Symmetry::analyze_magnetic_group_soc(const Atom* atoms, const Statistics& st, const ModuleBase::Matrix3& latvec)
{
    // Restrict the space group to the unitary magnetic subgroup (nspin=4 / SOC):
    // operation g survives if it preserves the magnetic configuration as a pseudovector,
    // i.e. W(g) m_i = m_{g(i)} for every atom, with W(g) = spin_so3(gmatc). 
    // Operations that reverse the moment (only symmetries together with time reversal) are dropped, 
    // so they are no longer applied in k-reduction or density symmetrization. 
    // Non-magnetic (m_i=0) keeps every operation.
    const ModuleBase::Matrix3 ilatvec = latvec.Inverse();
    std::vector<int> keep;
    keep.reserve(this->nrotk);
    int nrot_new = 0;
    for (int isym = 0; isym < this->nrotk; ++isym)
    {
        const ModuleBase::Matrix3 gmatc = ilatvec * this->gmatrix[isym] * latvec;
        const ModuleBase::Matrix3 W = ModuleSymmetry::SpinRotation::spin_so3(gmatc);
        bool ok = true;
        for (int iat = 0; iat < this->nat && ok; ++iat)
        {
            const ModuleBase::Vector3<double>& m = atoms[st.iat2it[iat]].m_loc_[st.iat2ia[iat]];
            // pseudovector-rotated moment W*m (column-vector convention: m'^i = W_ij m^j)
            const double mx = W.e11 * m.x + W.e12 * m.y + W.e13 * m.z;
            const double my = W.e21 * m.x + W.e22 * m.y + W.e23 * m.z;
            const double mz = W.e31 * m.x + W.e32 * m.y + W.e33 * m.z;
            const int jat = this->get_rotated_atom(isym, iat);
            const ModuleBase::Vector3<double>& mj = atoms[st.iat2it[jat]].m_loc_[st.iat2ia[jat]];
            if (!this->equal(mx, mj.x) || !this->equal(my, mj.y) || !this->equal(mz, mj.z)) { ok = false; }
        }
        if (ok)
        {
            keep.push_back(isym);
            if (isym < this->nrot) { ++nrot_new; }   // pure point-group rotations are the first nrot ops
        }
    }

    const int nrotk_new = static_cast<int>(keep.size());
    if (nrotk_new == this->nrotk) { return; }   // nothing removed (non-magnetic or fully-preserving group)

    // compact the operation arrays in ascending order (keeps the rotations-first layout).
    for (int i = 0; i < nrotk_new; ++i)
    {
        const int isym = keep[i];
        if (i != isym)
        {
            this->gmatrix[i] = this->gmatrix[isym];
            this->kgmatrix[i] = this->kgmatrix[isym];
            this->gtrans[i] = this->gtrans[isym];
            this->isym_rotiat_[i] = this->isym_rotiat_[isym];
        }
    }
    this->isym_rotiat_.resize(nrotk_new);
    this->nrot = nrot_new;
    this->nrotk = nrotk_new;

    // refresh the point-/space-group labels for the reduced (unitary magnetic) group
    this->pointgroup(this->nrot, this->pgnumber, this->pgname, this->gmatrix, GlobalV::ofs_running, nullptr);
    this->pointgroup(this->nrotk, this->spgnumber, this->spgname, this->gmatrix, GlobalV::ofs_running, nullptr);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "MAGNETIC POINT GROUP (unitary, nspin=4)", this->pgname);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "MAGNETIC SPACE GROUP OPERATIONS", this->nrotk);
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


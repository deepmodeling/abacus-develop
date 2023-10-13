#include "write_orb_info.h"

#include "module_base/name_angular.h"
#include "module_cell/atom_spec.h"

void ModuleIO::write_orb_info(const UnitCell* ucell)
{
    std::stringstream os;
    os << GlobalV::global_out_dir << "Orbital";
    std::ofstream out(os.str().c_str());
    out << std::setw(5) << "#io" << std::setw(8) << "spec" << std::setw(5) << "l" << std::setw(5) << "m" << std::setw(5)
        << "z" << std::setw(5) << "sym" << std::endl;

    for (int i = 0; i < ucell->nat; i++)
    {
        int t = ucell->iat2it[i];
        Atom* atom1 = &ucell->atoms[t];
        for (int j = 0; j < atom1->nw; ++j)
        {
            const int L1 = atom1->iw2l[j];
            const int N1 = atom1->iw2n[j];
            const int m1 = atom1->iw2m[j];
            out << std::setw(5) << i << std::setw(8) << ucell->atoms[t].label << std::setw(5) << L1
                << std::setw(5) << m1 << std::setw(5) << N1 + 1 << std::setw(15) << ModuleBase::Name_Angular[L1][m1]
                << std::endl;
        }
    }
    out << std::endl << std::endl;
    out << std::setw(5) << "#io" << std::setw(2) << "=" << std::setw(2) << "Orbital index in supercell" << std::endl;
    out << std::setw(5) << "#spec" << std::setw(2) << "=" << std::setw(2) << "Atomic species label" << std::endl;
    out << std::setw(5) << "#l" << std::setw(2) << "=" << std::setw(2) << "Angular mumentum quantum number" << std::endl;
    out << std::setw(5) << "#m" << std::setw(2) << "=" << std::setw(2) << "Magnetic quantum number" << std::endl;
    out << std::setw(5) << "#z" << std::setw(2) << "=" << std::setw(2) << "Zeta index of orbital" << std::endl;
    out << std::setw(5) << "#sym" << std::setw(2) << "=" << std::setw(2) << "Symmetry name of real orbital" << std::endl;
    out.close();

    return;
}

#include "ucell_io.h"
#include "source_base/constants.h"

#include <sstream>

namespace ModuleIO {

void UcellIO::write_ucell(std::ofstream& ofs, const UnitCell* ucell)
{
    // write the UnitCell information
    ofs << " " << ucell->latName << std::endl;
    ofs << " " << ucell->lat0 * ModuleBase::BOHR_TO_A << std::endl;
    ofs << " " << ucell->latvec.e11 << " " << ucell->latvec.e12 << " " << ucell->latvec.e13 << std::endl;
    ofs << " " << ucell->latvec.e21 << " " << ucell->latvec.e22 << " " << ucell->latvec.e23 << std::endl;
    ofs << " " << ucell->latvec.e31 << " " << ucell->latvec.e32 << " " << ucell->latvec.e33 << std::endl;
    for (int it = 0; it < ucell->ntype; it++)
    {
        ofs << " " << ucell->atoms[it].label;
    }
    ofs << std::endl;
    for (int it = 0; it < ucell->ntype; it++)
    {
        ofs << " " << ucell->atoms[it].na;
    }
    ofs << std::endl;
    ofs << " Direct" << std::endl;
    for (int it = 0; it < ucell->ntype; it++)
    {
        Atom* atom = &ucell->atoms[it];
        ofs << std::setprecision(15);
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            ofs << " " << atom->taud[ia].x << " " << atom->taud[ia].y << " " << atom->taud[ia].z << std::endl;
        }
    }
}

bool UcellIO::read_ucell(std::ifstream& ifs)
{
    std::string tmp;

    // --- six metadata lines: latName + lat0 + three lattice vectors + atom labels ---
    for (int i = 0; i < 6; ++i)
    {
        if (!std::getline(ifs, tmp))
        {
            return false;
        }
    }

    // --- per-type atom counts, must yield at least one positive integer ---
    if (!std::getline(ifs, tmp))
    {
        return false;
    }
    std::istringstream iss(tmp);
    int natom = 0;
    int total_natom = 0;
    while (iss >> natom)
    {
        if (natom < 0)
        {
            // Negative counts make the block size indeterminate.
            return false;
        }
        total_natom += natom;
    }
    if (total_natom <= 0)
    {
        // A malformed count line leaves total_natom at zero; the stream would
        // then be sitting on the coordinate rows, and a subsequent small-matrix
        // read could happily consume them as density-matrix entries. Reject
        // this case explicitly so the caller can abort with a clean diagnostic.
        return false;
    }

    // --- "Direct" marker + total_natom coordinate lines ---
    for (int i = 0; i < total_natom + 1; ++i)
    {
        if (!std::getline(ifs, tmp))
        {
            return false;
        }
    }
    return true;
}

} // namespace ModuleIO

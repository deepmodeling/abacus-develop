#include "read_orb.h"
#include "module_base/formatter.h"

namespace elecstate {
    void read_orb_file(int it, std::string &orb_file, std::ofstream &ofs_running, Atom* atom)
    {
    // the maximum L is 9 like cc-pV9Z, according to the 
    // basissetexchange https://www.basissetexchange.org/
    // there is no orbitals with L>9 presently
    const std::string spectrum = "SPDFGHIKLM";
    std::ifstream ifs(orb_file.c_str(), std::ios::in);  // pengfei 2014-10-13
    // mohan add return 2021-04-26
    if (!ifs)
    {
        std::cout << " Element index " << it+1 << std::endl;
        std::cout << " orbital file: " << orb_file << std::endl;
        ModuleBase::WARNING_QUIT("UnitCell::read_orb_file", 
                                "ABACUS Cannot find the ORBITAL file (basis sets)");
    }
    std::string word;
    atom->nw = 0;
    while (ifs.good())
    {
        ifs >> word;
        if (word == "Element")         // pengfei Li 16-2-29
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->label_orb);
        }
        if (word == "Lmax")
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->nwl);
            atom->l_nchi.resize(atom->nwl+1, 0);
        }
        // assert(atom->nwl<10); // cannot understand why restrict the maximum value of atom->nwl
        if (word == "Cutoff(a.u.)")         // pengfei Li 16-2-29
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->Rcut);
        }
        if (FmtCore::endswith(word, "orbital-->"))
        {
            bool valid = false;
            for (int i = 0; i < spectrum.size(); i++)
            {
                if (word == spectrum.substr(i, 1) + "orbital-->")
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->l_nchi[i]);
                    atom->nw += (2*i + 1) * atom->l_nchi[i];
                    std::stringstream ss;
                    ss << "L=" << i << ", number of zeta";
                    ModuleBase::GlobalFunc::OUT(ofs_running,ss.str(),atom->l_nchi[i]);
                    valid = true;
                    break;
                }
            }
            if (!valid)
            {
                ModuleBase::WARNING_QUIT("UnitCell::read_orb_file", 
                                         "ABACUS does not support numerical atomic orbital with L > 9, "
                                         "or an invalid orbital label is found in the ORBITAL file.");
            }
        }
    }
    ifs.close();
    if(!atom->nw)
    {
        ModuleBase::WARNING_QUIT("UnitCell::read_orb_file","get nw = 0");
    }
}

}
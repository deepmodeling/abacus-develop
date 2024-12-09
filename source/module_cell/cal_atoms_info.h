#ifndef CAL_ATOMS_INFO_H
#define CAL_ATOMS_INFO_H
#include "module_parameter/parameter.h"
#include "module_elecstate/cal_nelec_nband.h"
class CalAtomsInfo
{
  public:
    CalAtomsInfo(){};
    ~CalAtomsInfo(){};

    /**
     * @brief Calculate the atom information from pseudopotential to set Parameter
     *
     * @param atoms [in] Atom pointer
     * @param ntype [in] number of atom types
     * @param para [out] Parameter object
     */
    void cal_atoms_info(const Atom* atoms, const int& ntype, Parameter& para)
    {
        // calculate initial total magnetization when NSPIN=2
        if (para.inp.nspin == 2 && !para.globalv.two_fermi)
        {
            for (int it = 0; it < ntype; ++it)
            {
                for (int ia = 0; ia < atoms[it].na; ++ia)
                {
                    para.input.nupdown  += atoms[it].mag[ia];
                }
            }
            GlobalV::ofs_running << " The readin total magnetization is " << para.inp.nupdown  << std::endl;
        }

        if (!para.inp.use_paw)
        {
            // decide whether to be USPP
            for (int it = 0; it < ntype; ++it)
            {
                if (atoms[it].ncpp.tvanp)
                {
                    para.sys.use_uspp = true;
                }
            }
    
            // calculate the total number of local basis
            para.sys.nlocal = 0;
            for (int it = 0; it < ntype; ++it)
            {
                const int nlocal_it = atoms[it].nw * atoms[it].na;
                if (para.inp.nspin != 4)
                {
                    para.sys.nlocal += nlocal_it;
                }
                else
                {
                    para.sys.nlocal += nlocal_it * 2; // zhengdy-soc
                }
            }
        }

        // calculate the total number of electrons
        elecstate::cal_nelec(atoms, ntype, para.input.nelec);

        // autoset and check GlobalV::NBANDS
        std::vector<double> nelec_spin(2, 0.0);
        if (para.inp.nspin == 2)
        {
            nelec_spin[0] = (para.inp.nelec + para.inp.nupdown ) / 2.0;
            nelec_spin[1] = (para.inp.nelec - para.inp.nupdown ) / 2.0;
        }
        elecstate::cal_nbands(para.inp.nelec, para.sys.nlocal, nelec_spin, para.input.nbands);
        return;
    }
};
#endif
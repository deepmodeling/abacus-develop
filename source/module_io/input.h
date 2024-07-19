#ifndef INPUT_H
#define INPUT_H

#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#include "module_base/vector3.h"
#include "module_parameter/md_parameter.h"

class Input
{
  public:
    ~Input()
    {
        delete[] hubbard_u;
        delete[] orbital_corr;
    }
    
    // They will be removed.
    int cond_dtbatch;  
    int nche_sto;
    double md_tfirst;
    double bessel_nao_rcut;      // radial cutoff for spherical bessel functions(a.u.)
    MD_para mdp;
    int* orbital_corr = nullptr; ///< which correlated orbitals need corrected ;
    double* hubbard_u = nullptr; ///< Hubbard Coulomb interaction parameter U(ev)
    std::string stru_file;     // file contains atomic positions -- xiaohui modify

};

extern Input INPUT;
#endif // INPUT
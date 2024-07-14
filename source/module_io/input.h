#ifndef INPUT_H
#define INPUT_H

#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "input_conv.h"
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
    std::string bands_to_print_; // specify the bands to be calculated in the get_pchg

  public:
    // Return the const string pointer of private member bands_to_print_
    // Not recommended to use this function directly, use get_out_band_kb()
    // instead
    const std::string* get_bands_to_print() const
    {
        return &bands_to_print_;
    }
    // Return parsed bands_to_print_ as a vector of integers
    std::vector<int> get_out_band_kb() const
    {
        std::vector<int> out_band_kb;
        Input_Conv::parse_expression(bands_to_print_, out_band_kb);
        return out_band_kb;
    }
};

extern Input INPUT;
#endif // INPUT
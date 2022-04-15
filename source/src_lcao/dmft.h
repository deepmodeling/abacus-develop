//==========================================================
// Author: Xin Qu
// DATE : 2021-08
//==========================================================
#ifndef DFT_DMFT_INTERFACE
#define DFT_DMFT_INTERFACE

#include "../input.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../module_base/complexmatrix.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/local_orbital_charge.h"

#include <string>
#include <vector>

namespace ModuleDMFT
{
  class DFT_DMFT_interface
  {
    public:
    DFT_DMFT_interface(){;}
    ~DFT_DMFT_interface(){;}

    public:
    void init(Input& in, UnitCell_pseudo &cell);
    void out_to_dmft(Local_Orbital_wfc &lowf,
        LCAO_Matrix &lm);
    void out_kvector();
    void out_correlated_atom_info();
    void out_eigen_vector(Local_Orbital_wfc &lowf);

    static int mag_num2m_index(const int m);

    //TESTS
    void out_Sk(LCAO_Matrix &lm);

    private:
    void out_k_weight();
    void out_bands(double** ekb, const double Ef, const double Nelec);

    private:
    std::string out_path;

    std::vector<int> corr_L; //corr_L[it]
    std::vector<double> U;  //U[it]
    std::vector<double> J;  //J[it]

    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> iatlnmipol2iwt;   //iatlnm2iwt[iat][l][n][m][ipol]
    
  };
}

namespace GlobalC
{
  extern ModuleDMFT::DFT_DMFT_interface dmft;
}

#endif
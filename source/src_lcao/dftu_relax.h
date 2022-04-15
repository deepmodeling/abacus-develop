//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_RELAX_H
#define DFTU_RELAX_H

#include <vector>
#include "dftu_yukawa.h"
#include "../module_base/complexmatrix.h"
#include "src_lcao/LCAO_matrix.h"
using namespace std;

//==========================================================
// CLASS :
// NAME : DFTU_RELAX
//==========================================================
namespace ModuleDFTU{
class DFTU_RELAX : public DFTU_Yukawa
{

public:

    DFTU_RELAX();
    ~DFTU_RELAX();

    void force_stress(std::vector<ModuleBase::matrix>& dm_gamma,
        std::vector<ModuleBase::ComplexMatrix>& dm_k,
        LCAO_Matrix &lm);
    void cal_force_k(const int ik, const std::complex<double>* rho_VU);
    void cal_stress_k(const int ik, const std::complex<double>* rho_VU);
    void cal_force_gamma(const double* rho_VU);
    void cal_stress_gamma(const double* rho_VU);

    void fold_dSR_gamma(const int dim1, const int dim2, double* dSR_gamma);
    void fold_dSm_k(const int ik, const int dim, std::complex<double>* dSm_k);
    void fold_dSR_k(const int ik, const int dim1, const int dim2, std::complex<double>* dSR_gamma);

    void cal_VU_pot_mat_complex(const int spin, const bool newlocale, std::complex<double>* VU);
    void cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU);

    double get_onebody_eff_pot(
      const int T, const int iat,
	    const int L, const int N, const int spin, 
	    const int m0,
      const int m1,
      const int type, const bool newlocale);

    //forces and stress
    std::vector<std::vector<double>> force_dftu;      //force_dftu[iat][dim] 
    std::vector<std::vector<double>> stress_dftu;

    //transform between iwt index and it, ia, L, N and m index
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> iatlnmipol2iwt;   //iatlnm2iwt[iat][l][n][m][ipol]
    std::vector<int> iwt2it;                               //iwt2it[iwt]
    std::vector<int> iwt2l;                                //iwt2l[iwt]
    std::vector<int> iwt2n;                                //iwt2n[iwt]
    std::vector<int> iwt2m;                                //iwt2m[iwt]
    std::vector<int> iwt2ipol;                             //iwt2ipol[iwt]
    std::vector<int> iat2it; 

    //local occupancy matrix of the correlated subspace
    //locale: the out put local occupation number matrix of correlated electrons in the current electronic step
    //locale_save: the input local occupation number matrix of correlated electrons in the current electronic step
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale;            // locale[iat][l][n][spin](m1,m2)
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale_save;       // locale_save[iat][l][n][spin](m1,m2)
protected:
    LCAO_Matrix* LM;
};
}
#endif

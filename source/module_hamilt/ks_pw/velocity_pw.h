#ifndef VELOCITY_PW_H
#define VELOCITY_PW_H
#include "operator_pw.h"
#include "module_cell/unitcell_pseudo.h"
#include "src_pw/VNL_in_pw.h"
#include "module_pw/pw_basis_k.h"
namespace hamilt
{

//velocity operator mv = im/\hbar * [H,r] =  p + im/\hbar [V_NL, r] 
class Velocity
{
    public:
    Velocity(
        const ModulePW::PW_Basis_K* wfcpw_in,
        const int* isk_in,
        pseudopot_cell_vnl* ppcell_in,
        const UnitCell_pseudo* ucell_in,
        const bool nonlocal_in = true
    );

    ~Velocity(){};

    void init(const int ik_in);

    //input: n_npwx, output: 3*n_npwx
    void act
    (
        const psi::Psi<std::complex<double>> *psi_in, 
        const int n_npwx, 
        const std::complex<double>* tmpsi_in, 
        std::complex<double>* tmhpsi,
        const bool add = false
    )const;

    private:
    bool nonlocal = true;
    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    const int* isk = nullptr;

    pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell_pseudo* ucell = nullptr;

    int ik;

    double tpiba;
};
}
#endif
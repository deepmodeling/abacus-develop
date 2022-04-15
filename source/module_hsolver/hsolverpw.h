#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "src_pw/pw_basis.h"

namespace ModuleHSolver
{

class HSolverPW : public HSolver
{
    public:

    HSolverPW(const PW_Basis* pbas_in){
        this->pbas = pbas_in;
        /*this->init(pbas_in);*/}

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/
    
    void solve(
        ModuleHamilt::Hamilt* pHamilt, 
        ModulePsi::Psi<std::complex<double>>& psi, 
        ModuleElecS::ElecState* pes) override;
    

    private:
    void hamiltSolvePsiK(ModuleHamilt::Hamilt* hm, ModulePsi::Psi<std::complex<double>>& psi, double* eigenvalue);

    const PW_Basis* pbas;

    //calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<double> h_diag, const int npw, const double* g2kin);
};

}

#endif
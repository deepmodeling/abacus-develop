#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "hsolver.h"
#include "src_lcao/local_orbital_wfc.h"

namespace ModuleHSolver
{

class HSolverLCAO : public HSolver
{
    public:

    HSolverLCAO(
        Local_Orbital_wfc* lowf_in,
        double** ekb_in
    )
    {
        this->lowf = lowf_in;
        this->ekb = ekb_in;
    }

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

    void solve(
        ModuleHamilt::Hamilt* pHamilt, 
        ModulePsi::Psi<double>& psi, 
        ModuleElecS::ElecState* pes) override;

    static int out_wfc_lcao;
    static int out_mat_hs; // mohan add 2010-09-02
    static int out_mat_hsR; // LiuXh add 2019-07-16
    

    private:
    void hamiltSolvePsiK(ModuleHamilt::Hamilt* hm, ModulePsi::Psi<std::complex<double>>& psi, double* eigenvalue);
    void hamiltSolvePsiK(ModuleHamilt::Hamilt* hm, ModulePsi::Psi<double>& psi, double* eigenvalue);

    template<typename T>
    void solveTemplate(
        ModuleHamilt::Hamilt* pHamilt, 
        ModulePsi::Psi<T>& psi, 
        ModuleElecS::ElecState* pes
    );
    /*void solveTemplate(
        ModuleHamilt::Hamilt* pHamilt, 
        ModulePsi::Psi<std::complex<double>>& psi, 
        ModuleElecS::ElecState* pes
    );*/

    Local_Orbital_wfc* lowf;
    double** ekb;
};

}

#endif
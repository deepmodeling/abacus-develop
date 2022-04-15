#ifndef MODULEHAMILT_H
#define MODULEHAMILT_H

#include<complex>
#include "module_psi/psi.h"

#include "matrixblock.h"

namespace ModuleHamilt
{


class Hamilt
{
    public:
    
    //construct Hamiltonian matrix with inputed electonic density 
    virtual void constructHamilt(const int iter, const MatrixBlock<double> rho)=0;
    
    //for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik)=0;
    
    //core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const ModulePsi::Psi<std::complex<double>>& psi, ModulePsi::Psi<std::complex<double>>& hpsi)const=0;
    
    //core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    virtual void matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)=0;
    virtual void matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)=0;
    
    protected:
    //array, save operations from each operators
    //would be implemented later
    //vector<Operator*> p_operators;
};

}


#endif
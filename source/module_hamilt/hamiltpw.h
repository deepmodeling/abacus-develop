#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "hamilt.h"


namespace ModuleHamilt
{

class HamiltPW : public Hamilt
{
    public:
    HamiltPW(){};
    ~HamiltPW(){};

    //construct Hamiltonian matrix with inputed electonic density 
    void constructHamilt(const int iter, const MatrixBlock<double> rho) override
    {
        this->ch_mock();
    };
    
    //for target K point, update consequence of hPsi() and matrix()
    void updateHk(int ik) override
    {
        this->hk_mock();
    };
    
    //core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const ModulePsi::Psi<std::complex<double>>& psi, ModulePsi::Psi<std::complex<double>>& hpsi) const override
    {
        this->hpsi_mock(psi, hpsi);
    };
    
    //core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    //not used in PW base
    void matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)override{return;}
    void matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)override{return;}

    private:
    void ch_mock();
    void hk_mock();
    void hpsi_mock(const ModulePsi::Psi<std::complex<double>>& psi, ModulePsi::Psi<std::complex<double>>& hpsi) const;
};

}


#endif
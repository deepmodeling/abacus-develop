#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "hsolver.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

class HSolverLCAO : public HSolver<std::complex<double>>
{
  public:
    HSolverLCAO(const Parallel_Orbitals* ParaV_in)
    {
      this->classname = "HSolverPW"; 
      this->ParaV = ParaV_in;
      }
    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/

    void solve(hamilt::Hamilt<std::complex<double>>* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

    void solve(hamilt::Hamilt<std::complex<double>>* pHamilt, psi::Psi<double>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

    static int out_mat_hs; // mohan add 2010-09-02
    static int out_mat_hsR; // LiuXh add 2019-07-16
    static int out_mat_t;
    static int out_mat_dh;

  private:
    void hamiltSolvePsiK(hamilt::Hamilt<std::complex<double>>* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue);
    void hamiltSolvePsiK(hamilt::Hamilt<std::complex<double>>* hm, psi::Psi<double>& psi, double* eigenvalue);

    template <typename T> void solveTemplate(hamilt::Hamilt<std::complex<double>>* pHamilt, psi::Psi<T>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge);
    /*void solveTemplate(
        hamilt::Hamilt* pHamilt,
        psi::Psi<std::complex<double>>& psi,
        elecstate::ElecState* pes
    );*/
    const Parallel_Orbitals* ParaV;
};

} // namespace hsolver

#endif
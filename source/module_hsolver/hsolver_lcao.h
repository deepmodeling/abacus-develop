#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "hsolver.h"
#include "module_orbital/parallel_orbitals.h"

namespace hsolver
{

class HSolverLCAO : public HSolver
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

    void solve(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

    void solve(hamilt::Hamilt* pHamilt, psi::Psi<double>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

    static int out_mat_hs_k; // mohan add 2010-09-02
    static int out_mat_hs_r; // LiuXh add 2019-07-16

  private:
    void hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue);
    void hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<double>& psi, double* eigenvalue);

    template <typename T> void solveTemplate(hamilt::Hamilt* pHamilt, psi::Psi<T>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge);
    /*void solveTemplate(
        hamilt::Hamilt* pHamilt,
        psi::Psi<std::complex<double>>& psi,
        elecstate::ElecState* pes
    );*/
    const Parallel_Orbitals* ParaV;
};

} // namespace hsolver

#endif
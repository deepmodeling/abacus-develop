#ifndef OUTPUT_DM1_H
#define OUTPUT_DM1_H

#include "module_cell/klist.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"

namespace ModuleIO
{

/// @brief the output interface to write the sparse density matrix
class Output_DM1
{
  public:
    Output_DM1(int nspin,
               int istep,
               Local_Orbital_Charge& LOC,
               Record_adj& RA,
               K_Vectors& kv,
               const elecstate::DensityMatrix<std::complex<double>, double>* DM);

    void write();

  private:
    int _nspin;
    int _istep;
    Local_Orbital_Charge& loc;
    Record_adj& _RA;
    K_Vectors& _kv;
    const elecstate::DensityMatrix<std::complex<double>, double>* _DM;
};

} // namespace ModuleIO

#endif

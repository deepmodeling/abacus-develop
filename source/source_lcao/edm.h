#ifndef LCAO_EDM_H
#define LCAO_EDM_H

#include "source_base/global_function.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_psi/psi.h"

template <typename T>
class Force_Stress_LCAO;

template <typename T>
class CalEDM
{
  public:
    friend class Force_Stress_LCAO<T>;

    CalEDM(){};
    ~CalEDM(){};

  private:
    const Parallel_Orbitals* ParaV = nullptr;

    elecstate::Potential* pot = nullptr;

    elecstate::DensityMatrix<T, double> cal_edm(const elecstate::ElecState* pelec,
                                                const psi::Psi<T>& psi,
                                                const elecstate::DensityMatrix<T, double>& dm,
                                                const K_Vectors& kv,
                                                const Parallel_Orbitals& pv,
                                                const int& nspin,
                                                const int& nbands,
                                                const UnitCell& ucell,
                                                Record_adj& ra) const;
};

#endif

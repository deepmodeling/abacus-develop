#ifndef CAL_R_OVERLAP_R_H
#define CAL_R_OVERLAP_R_H

#include "source_base/vector3.h"

#include <memory>
#include <string>
#include <vector>

class UnitCell;
class Parallel_Orbitals;
class LCAO_Orbitals;
class Grid_Driver;
class PosOpBasis;
class PosOpCalc;
class PosOpWriter;

// output r_R matrix, added by Jingan
// Facade that delegates to PosOpBasis / PosOpCalc / PosOpWriter.
class Position_op
{

  public:
    Position_op();
    ~Position_op();

    double kmesh_times = 4;
    double sparse_threshold = 1e-10;
    bool binary = false;

    void init(const UnitCell& ucell,
              const Parallel_Orbitals& pv,
              const LCAO_Orbitals& orb,
              const bool cal_force,
              const int nlocal);
    void init_nonlocal(const UnitCell& ucell,
                       const Parallel_Orbitals& pv,
                       const LCAO_Orbitals& orb,
                       const bool cal_force,
                       const int nlocal);
    ModuleBase::Vector3<double> get_psi_r_psi(const ModuleBase::Vector3<double>& R1,
                                              const int& T1,
                                              const int& L1,
                                              const int& m1,
                                              const int& N1,
                                              const ModuleBase::Vector3<double>& R2,
                                              const int& T2,
                                              const int& L2,
                                              const int& m2,
                                              const int& N2);
    ModuleBase::Vector3<double> get_psi_r_gradpsi(const ModuleBase::Vector3<double>& R1,
                                                  const int& T1,
                                                  const int& L1,
                                                  const int& m1,
                                                  const int& N1,
                                                  const ModuleBase::Vector3<double>& R2,
                                                  const int& T2,
                                                  const int& L2,
                                                  const int& m2,
                                                  const int& N2,
                                                  const ModuleBase::Vector3<double>& Efield,
                                                  const ModuleBase::Vector3<double>& dR);
    void get_psi_r_beta(const UnitCell& ucell,
                        std::vector<std::vector<double>>& nlm,
                        const ModuleBase::Vector3<double>& R1,
                        const int& T1,
                        const int& L1,
                        const int& m1,
                        const int& N1,
                        const ModuleBase::Vector3<double>& R2,
                        const int& T2);
    void out_rR(const UnitCell& ucell,
                const Grid_Driver& gd,
                const int& istep,
                const int precision,
                const std::string& global_out_dir,
                const std::string& global_matrix_dir,
                const std::string& calculation,
                const bool out_app_flag,
                const int nlocal,
                const int npol);

  private:
    std::unique_ptr<PosOpBasis> basis_;
    std::unique_ptr<PosOpCalc> calc_;
    std::unique_ptr<PosOpWriter> writer_;
};
#endif

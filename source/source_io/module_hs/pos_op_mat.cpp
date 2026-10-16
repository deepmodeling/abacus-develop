#include "pos_op_mat.h"

#include "pos_op_basis.h"
#include "pos_op_calc.h"
#include "pos_op_writer.h"
#include "source_base/timer.h"

Position_op::Position_op()
{
}

Position_op::~Position_op()
{
}

void Position_op::init(const UnitCell& ucell,
                           const Parallel_Orbitals& pv,
                           const LCAO_Orbitals& orb,
                           const bool cal_force,
                           const int nlocal)
{
    basis_ = std::unique_ptr<PosOpBasis>(new PosOpBasis());
    basis_->build(ucell, orb, cal_force, nlocal);

    calc_ = std::unique_ptr<PosOpCalc>(new PosOpCalc(*basis_));
    writer_ = std::unique_ptr<PosOpWriter>(new PosOpWriter(*basis_, *calc_, pv));
}

void Position_op::init_nonlocal(const UnitCell& ucell,
                                    const Parallel_Orbitals& pv,
                                    const LCAO_Orbitals& orb,
                                    const bool cal_force,
                                    const int nlocal)
{
    basis_ = std::unique_ptr<PosOpBasis>(new PosOpBasis());
    basis_->build_nonlocal(ucell, orb, cal_force, nlocal);

    calc_ = std::unique_ptr<PosOpCalc>(new PosOpCalc(*basis_));
    writer_ = std::unique_ptr<PosOpWriter>(new PosOpWriter(*basis_, *calc_, pv));
}

ModuleBase::Vector3<double> Position_op::get_psi_r_psi(const ModuleBase::Vector3<double>& R1,
                                                           const int& T1,
                                                           const int& L1,
                                                           const int& m1,
                                                           const int& N1,
                                                           const ModuleBase::Vector3<double>& R2,
                                                           const int& T2,
                                                           const int& L2,
                                                           const int& m2,
                                                           const int& N2)
{
    return calc_->pos_matrix(R1, T1, L1, m1, N1, R2, T2, L2, m2, N2);
}

ModuleBase::Vector3<double> Position_op::get_psi_r_gradpsi(const ModuleBase::Vector3<double>& R1,
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
                                                               const ModuleBase::Vector3<double>& dR)
{
    return calc_->pos_grad_matrix(R1, T1, L1, m1, N1, R2, T2, L2, m2, N2, Efield, dR);
}

void Position_op::get_psi_r_beta(const UnitCell& ucell,
                                     std::vector<std::vector<double>>& nlm,
                                     const ModuleBase::Vector3<double>& R1,
                                     const int& T1,
                                     const int& L1,
                                     const int& m1,
                                     const int& N1,
                                     const ModuleBase::Vector3<double>& R2,
                                     const int& T2)
{
    calc_->pos_beta_matrix(ucell, nlm, R1, T1, L1, m1, N1, R2, T2);
}

void Position_op::out_rR(const UnitCell& ucell,
                             const Grid_Driver& gd,
                             const int& istep,
                             const int precision,
                             const std::string& global_out_dir,
                             const std::string& global_matrix_dir,
                             const std::string& calculation,
                             const bool out_app_flag,
                             const int nlocal,
                             const int npol,
                             std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Position_op", "out_rR");
    ModuleBase::timer::start("Position_op", "out_rR");

    writer_->out_lat_r(ucell,
                       gd,
                       istep,
                       precision,
                       global_out_dir,
                       global_matrix_dir,
                       calculation,
                       out_app_flag,
                       nlocal,
                       npol,
                       sparse_threshold,
                       binary,
                       ofs_running);

    ModuleBase::timer::end("Position_op", "out_rR");
}

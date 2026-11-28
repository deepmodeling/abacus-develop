#ifndef MOMENTUM_OP_MAT_H
#define MOMENTUM_OP_MAT_H

#include <string>

struct Input_para;
class UnitCell;
class Grid_Driver;
class Parallel_Orbitals;
class LCAO_Orbitals;
class TwoCenterBundle;

namespace Momentum_op
{

/**
 * @brief Compute and write the momentum operator matrix p = -i nabla + i [V_nl, r]
 * in the local atomic basis as three sparse R-space matrices (px, py, pz).
 *
 * The matrix elements are evaluated with the existing Velocity_op kernel from
 * module_momentum_op, which contains both the kinetic gradient term and the
 * nonlocal-pseudopotential commutator term required for TDSE applications.
 *
 * The three Cartesian components are written as CSR files named
 * prx_nao.csr, pry_nao.csr and prz_nao.csr (or prxg{step+1}_nao.csr etc.
 * when out_app_flag = 2), using the same header and layout as the H(R)/S(R)
 * sparse output.
 *
 * @param istep            current ionic step index (zero-based; the g{step} suffix uses istep + 1;
 *                         pass -1 from the get_s workflow to skip the suffix)
 * @param ucell            unit cell information
 * @param gd               grid driver for atom-pair search
 * @param two_center_bundle two-center integrator bundle
 * @param pv               parallel orbitals for 2D block-cyclic distribution
 * @param orb              local atomic orbital descriptor
 * @param output_dir       directory prefix for the output files (e.g. "MAT/")
 * @param inp              INPUT parameters, including out_mat_p and out_app_flag
 */
void calculate(const int istep,
               const UnitCell& ucell,
               Grid_Driver& gd,
               TwoCenterBundle& two_center_bundle,
               const Parallel_Orbitals& pv,
               const LCAO_Orbitals& orb,
               const std::string& output_dir,
               const Input_para& inp);

} // namespace Momentum_op

#endif // MOMENTUM_OP_MAT_H

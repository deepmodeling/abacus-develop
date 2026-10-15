#ifndef POS_OP_WRITER_H
#define POS_OP_WRITER_H

#include "source_base/vector3.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "pos_op_basis.h"
#include "pos_op_calc.h"

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

class Grid_Driver;
class Parallel_Orbitals;

/**
 * @brief Output the position-operator matrix <phi_mu | r_hat | phi_nu>
 *        as a sparse matrix for each lattice vector R (lat_r).
 *
 * The writer loops over all periodic images, assembles the sparse blocks,
 * and writes the final rr.csr file via ModuleIO helpers.
 */
class PosOpWriter
{
  public:
    PosOpWriter(const PosOpBasis& basis, const PosOpCalc& calc, const Parallel_Orbitals& pv);

    void out_lat_r(const UnitCell& ucell,
                   const Grid_Driver& gd,
                   int istep,
                   int precision,
                   const std::string& global_out_dir,
                   const std::string& global_matrix_dir,
                   const std::string& calculation,
                   bool out_app_flag,
                   int nlocal,
                   int npol,
                   double sparse_threshold,
                   bool binary,
                   std::ofstream& ofs_running);

  private:
    const PosOpBasis& basis_;
    const PosOpCalc& calc_;
    const Parallel_Orbitals& pv_;
};

#endif // POS_OP_WRITER_H

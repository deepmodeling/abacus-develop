#include "rpa_lri_detail.h"
#include "librpa_stru_symmetry.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace
{
bool preserves_mesh(const ModuleBase::Matrix3& kg, const K_Vectors& kv, const UnitCell& cell)
{
    const double rows[3][3] = {{kg.e11, kg.e12, kg.e13}, {kg.e21, kg.e22, kg.e23}, {kg.e31, kg.e32, kg.e33}};
    // ABACUS and LibRPA rotate row vectors: k' = k * kg.  Therefore the
    // coefficient taking source component j on an nmp[j] mesh to target
    // component i on an nmp[i] mesh is kg[j][i] * nmp[i] / nmp[j].
    for (int j = 0; j < 3; ++j)
    {
        if (kv.nmp[j] <= 0) throw std::runtime_error("Spatial symmetry requires a uniform k-grid");
        for (int i = 0; i < 3; ++i)
        {
            const double step = rows[j][i] * kv.nmp[i] / kv.nmp[j];
            if (std::abs(step - std::round(step)) > 1e-8) return false;
        }
    }
    // Check the actual mesh offset as well as its spacing.
    const auto origin = kv.kvec_c_full.at(0) * cell.G.Inverse();
    const auto delta = origin * kg - origin;
    const double shifts[3] = {delta.x * kv.nmp[0], delta.y * kv.nmp[1], delta.z * kv.nmp[2]};
    for (const double shift : shifts)
    {
        if (std::abs(shift - std::round(shift)) > 1e-8) return false;
    }
    return true;
}

ModuleBase::Matrix3 anti_k_action(const ModuleBase::Matrix3& spatial)
{
    ModuleBase::Matrix3 action = spatial;
    action *= -1.0;
    return action;
}
}

void RpaLriDetail::write_stru_sym(std::ostream& output, const UnitCell& cell, const K_Vectors& kv, const int nspin)
{
    const auto& symm = cell.symm;
    std::vector<LibRpaSymmetryOperation> unitary_operations;
    std::vector<LibRpaSymmetryOperation> antiunitary_operations;
    std::vector<ModuleBase::Matrix3> unitary_rotations;
    std::vector<ModuleBase::Matrix3> antiunitary_rotations;
    for (int isym = 0; isym < symm.nrotk; ++isym)
    {
        // Use the same reciprocal-space matrices that ABACUS used to reduce
        // this mesh.  Reconstructing them from gmatrix can change the basis
        // convention for non-orthogonal primitive cells.
        if (preserves_mesh(symm.kgmatrix[isym], kv, cell))
        {
            unitary_operations.push_back(make_librpa_symmetry_operation(symm.gmatrix[isym], symm.gtrans[isym]));
            unitary_rotations.push_back(symm.gmatrix[isym]);
        }
    }
    if (unitary_operations.empty())
    {
        throw std::runtime_error("Missing mesh-preserving spatial symmetry operations");
    }

    if (nspin == 4 && symm.magnetic_nspin4)
    {
        for (int isym = 0; isym < symm.nrotk_anti; ++isym)
        {
            // The exported row stores the spatial part g. The k-point action
            // of Theta*g is -k*g, so the mesh test must include time reversal.
            if (preserves_mesh(anti_k_action(symm.kgmatrix_anti[isym]), kv, cell))
            {
                antiunitary_operations.push_back(
                    make_librpa_symmetry_operation(symm.gmatrix_anti[isym], symm.gtrans_anti[isym]));
                antiunitary_rotations.push_back(symm.gmatrix_anti[isym]);
            }
        }
    }

    write_librpa_symmetry_rows(output, unitary_operations, antiunitary_operations);
    if (nspin == 4)
    {
        const auto spin_operations
            = make_librpa_spin_symmetry_operations(cell.latvec, unitary_rotations, antiunitary_rotations);
        const int grey_group = symm.magnetic_nspin4 ? 0 : 1;
        write_librpa_spin_symmetry(output, grey_group, 1, spin_operations);
    }
}

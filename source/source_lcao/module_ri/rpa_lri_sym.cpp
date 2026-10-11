#include "rpa_lri_detail.h"
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
}

void RpaLriDetail::write_stru_sym(std::ostream& output, const UnitCell& cell, const K_Vectors& kv)
{
    const auto& symm = cell.symm;
    std::vector<int> operations;
    for (int isym = 0; isym < symm.nrotk; ++isym)
    {
        // Use the same reciprocal-space matrices that ABACUS used to reduce
        // this mesh.  Reconstructing them from gmatrix can change the basis
        // convention for non-orthogonal primitive cells.
        if (preserves_mesh(symm.kgmatrix[isym], kv, cell)) operations.push_back(isym);
    }
    if (operations.empty()) throw std::runtime_error("Missing mesh-preserving spatial symmetry operations");
    // Use the subgroup that preserves this mesh. The full crystal group may
    // rotate an anisotropic mesh outside itself and cannot accelerate its RPA.
    output << operations.size() << " row\n";
    for (const int isym : operations)
    {
        const auto& r = symm.gmatrix[isym];
        const double entries[9] = {r.e11, r.e12, r.e13, r.e21, r.e22, r.e23, r.e31, r.e32, r.e33};
        for (const double entry : entries)
        {
            if (!std::isfinite(entry) || std::abs(entry - std::round(entry)) > 1e-8)
                throw std::runtime_error("Nonintegral fractional symmetry rotation");
            output << std::setw(4) << static_cast<int>(std::round(entry));
        }
        const auto& t = symm.gtrans[isym];
        output << std::scientific << std::setprecision(15) << std::setw(24) << t.x
               << std::setw(24) << t.y << std::setw(24) << t.z << '\n';
    }
}

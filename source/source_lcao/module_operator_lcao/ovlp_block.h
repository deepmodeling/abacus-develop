#ifndef ABACUS_LCAO_OVLP_BLOCK_H
#define ABACUS_LCAO_OVLP_BLOCK_H
namespace ModuleBase { template <typename T> class Vector3; }
class UnitCell;
class Parallel_Orbitals;
class TwoCenterIntegrator;
namespace hamilt
{
// Add one ordered atom-pair overlap block in local row-major AO layout.
// displacement is ket-center minus bra-center in Bohr, including any lattice image
// or cross-geometry shift. No neighbor search, coordinate reconstruction or symmetrization.
// For npol=2, add the same spatial integral to the two spin-diagonal entries.
template <typename TR>
void cal_overlap_block(const UnitCell& cell, const TwoCenterIntegrator& integrator,
    int atom_i, int atom_j, const Parallel_Orbitals& distribution,
    const ModuleBase::Vector3<double>& displacement, TR* values);
}
#endif

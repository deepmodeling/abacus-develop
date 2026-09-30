#ifndef SURCHEM_INPUT_H
#define SURCHEM_INPUT_H

struct Input_para;
struct SurchemParameters;
class UnitCell;
class K_Vectors;

#include <string>

namespace ModuleSymmetry
{
class Symmetry;
}

// Adapter from explicit host inputs to solvent/PCC configuration; no global reads.
namespace ModuleSurchem
{
SurchemParameters make_parameters(const Input_para& input,
                                  const UnitCell& cell,
                                  const double electron_count,
                                  const bool use_uspp,
                                  const int pool_process_count);
void validate_kpoints(const SurchemParameters& parameters, const K_Vectors& kpoints);

// The PCC corrections are invariant only under symmetry operations that keep
// their frame: for pcc_2d the open lattice axis must map onto +-itself and no
// fractional translation may act along it; pcc_0d allows no pure translation.
// Returns an empty string when every analyzed operation is compatible, else
// a description of the first violation.
std::string pcc_symmetry_violation(const SurchemParameters& parameters,
                                   const ModuleSymmetry::Symmetry& symmetry);
// WARNING_QUIT on a pcc_symmetry_violation.
void validate_symmetry(const SurchemParameters& parameters, const ModuleSymmetry::Symmetry& symmetry);
} // namespace ModuleSurchem

#endif

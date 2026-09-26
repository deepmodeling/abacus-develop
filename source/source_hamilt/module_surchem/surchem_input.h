#ifndef SURCHEM_INPUT_H
#define SURCHEM_INPUT_H

struct Input_para;
struct SurchemParameters;
class UnitCell;
class K_Vectors;

// Adapter from explicit host inputs to solvent/PCC configuration; no global reads.
namespace ModuleSurchem
{
SurchemParameters make_parameters(const Input_para& input,
                                  const UnitCell& cell,
                                  const double electron_count,
                                  const bool use_uspp,
                                  const int pool_process_count);
void validate_kpoints(const SurchemParameters& parameters, const K_Vectors& kpoints);
} // namespace ModuleSurchem

#endif

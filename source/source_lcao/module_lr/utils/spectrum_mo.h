#ifndef ABACUS_LR_SPECTRUM_MO_H
#define ABACUS_LR_SPECTRUM_MO_H
#include <complex>
#include <string>
#include <vector>
#include "source_psi/psi.h"
class UnitCell;
class Grid_Driver;
class TwoCenterBundle;
class Parallel_Orbitals;
class Parallel_2D;
class K_Vectors;
namespace LR_Util
{
/// Compute C^dagger v C for the selected KS window; return replicated spin/cartesian/k blocks.
template<typename T>
std::vector<std::complex<double>> cal_velocity_mo(const UnitCell& ucell,
                                                const Grid_Driver& gd,
                                                const TwoCenterBundle& two_center_bundle,
                                                const Parallel_Orbitals& pmat,/*nbasis×nbasis*/
                                                const Parallel_2D& pc,/*nbasis×nbands*/
                                                const K_Vectors& kv,
                                                const psi::Psi<T>& psi_ks,
                                                const int nk,
                                                const int nspin_tmp,
                                                const int nbasis,
                                                const std::vector<int> nocc,
                                                const std::vector<int> nvirt);
void output_spectrum_mo(const std::vector<std::complex<double>>& out_spectrum_mo,
                        const std::string& filename,
                        const double* const eig_ks,
                        const int nk,
                        const int nspin_tmp,
                        const int KS_num,
                        const K_Vectors& kv);
void output_spectrum_mo_librpa(const std::vector<std::complex<double>>& out_spectrum_mo,
                        const std::string& filename,
                        const int nk, const int nspin_tmp, const int nbands, const int nbasis,
                        const int KS_num/*= nocc+nvirt, can be different from nbands*/);
}
#endif

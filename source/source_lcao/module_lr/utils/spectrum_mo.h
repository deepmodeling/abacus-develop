#ifndef ABACUS_LR_SPECTRUM_MO_H
#define ABACUS_LR_SPECTRUM_MO_H
#include <complex>
#include <iosfwd>
#include <string>
#include <vector>
#include "source_psi/psi.h"
class UnitCell;
class Grid_Driver;
class TwoCenterBundle;
class Parallel_Orbitals;
class Parallel_Blkcyclic; // Defined in source/source_base/parallel_blkcyc.h
class K_Vectors;
class LCAO_Orbitals;
namespace hamilt
{
template <typename T> class HContainer;
}
namespace LR_Util
{
/// Compute C^dagger v C for the selected KS window; return replicated spin/cartesian/k blocks.
template<typename T>
std::vector<std::complex<double>> cal_velocity_mo(const UnitCell& ucell,
                                                const Grid_Driver& gd,
                                                const TwoCenterBundle& two_center_bundle,
                                                const Parallel_Orbitals& pmat,/*nbasis×nbasis*/
                                                const Parallel_Blkcyclic& pc,/*nbasis×nbands*/
                                                const K_Vectors& kv,
                                                const psi::Psi<T>& psi_ks,
                                                const int nk,
                                                const int nspin_tmp,
                                                const int nbasis,
                                                const std::vector<int> nocc,
                                                const std::vector<int> nvirt);
/// Compute the same KS window with the full generalized velocity i[H,r].
template <typename T>
std::vector<std::complex<double>> cal_commutator_velocity_mo(
    const UnitCell& ucell, const Grid_Driver& gd, const LCAO_Orbitals& orb,
    const Parallel_Orbitals& pmat, const Parallel_Blkcyclic& pc, const K_Vectors& kv,
    const psi::Psi<T>& psi_ks, const int nk, const int nspin,
    const std::vector<int>& nocc, const std::vector<int>& nvirt, const bool cal_force,
    const hamilt::HContainer<double>& sr,
    const std::vector<const hamilt::HContainer<double>*>& hr);
void read_velocity_matrix_R(hamilt::HContainer<double>& matrix, const std::string& path,
                            const UnitCell& ucell, const int nbasis, const int rank, std::ostream& log);
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

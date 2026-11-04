#ifndef OUTPUT_MULLIKEN_H
#define OUTPUT_MULLIKEN_H
#include "source_base/complexmatrix.h"
#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/cell_index.h"

#include <map>
#include <vector>

namespace ModuleIO
{

template <typename TK>
class Output_Sk;

template <typename TK>
class Output_DMK;

/**
 * @brief Accumulate the nspin=4 Mulliken charge and Mx/My/Mz components.
 *
 * @details Decomposes one k-point mud = DM*S^T into the Pauli components using
 * the conj-first DM convention (DM = conj(P), see dm_from_psi.h):
 * Mx = Re(mud_ud + mud_du), My = Im(mud_ud - mud_du), Mz = Re(mud_uu - mud_dd).
 *
 * @param MecMulP accumulators [4 x nw/2]: 0 = charge, 1/2/3 = Mx/My/Mz
 * @param mud     DM*S^T stored in ScaLAPACK column-major layout
 * @param nw      total number of spinor orbitals (2 * number of spatial orbitals)
 * @param pv      orbital distribution mapping
 */
void collect_MW_nspin4(ModuleBase::matrix& MecMulP,
                       const ModuleBase::ComplexMatrix& mud,
                       int nw,
                       const Parallel_Orbitals* pv);

/// @brief the output interface to write the Mulliken population charges
template <typename TK>
class Output_Mulliken
{
  public:
    /// constructor of Output_Mulliken
    Output_Mulliken(Output_Sk<TK>* output_sk,
                    Output_DMK<TK>* output_dmk,
                    Parallel_Orbitals* ParaV,
                    CellIndex* cell_index,
                    const std::vector<int>& isk,
                    int nspin);
    /// the outer interface to write the Mulliken population charges
    void write(int istep, std::string out_dir);
    /// print atom mag to running log file
    void print_atom_mag(const std::vector<std::vector<double>>& atom_chg, std::ostream& os);
    /// get total charge
    std::vector<double> get_tot_chg();
    /// get atom charge
    std::vector<std::vector<double>> get_atom_chg();
    /// get orbital charge
    std::map<std::vector<int>, double> get_orb_chg();
    /// returun atom_mulliken for updateing STRU file
    std::vector<std::vector<double>> get_atom_mulliken(std::vector<std::vector<double>>& atom_chg);

  private:
    /******************************************************************
     * private functions
     *******************************************************************/
    /// write mulliken.txt for the case of nspin=1
    void write_mulliken_nspin1(int istep,
                               const std::vector<double>& tot_chg,
                               const std::vector<std::vector<double>>& atom_chg,
                               std::map<std::vector<int>, double> orb_chg,
                               std::ofstream& os);
    /// write mulliken.txt for the case of nspin=2
    void write_mulliken_nspin2(int istep,
                               const std::vector<double>& tot_chg,
                               const std::vector<std::vector<double>>& atom_chg,
                               std::map<std::vector<int>, double> orb_chg,
                               std::ofstream& os);
    /// write mulliken.txt for the case of nspin=4
    void write_mulliken_nspin4(int istep,
                               const std::vector<double>& tot_chg,
                               const std::vector<std::vector<double>>& atom_chg,
                               std::map<std::vector<int>, double> orb_chg,
                               std::ofstream& os);
    /// set nspin
    void set_nspin(int nspin_in);
    /// set orbital parallel info
    void set_ParaV(Parallel_Orbitals* ParaV_in);
    /// collect_mw from matrix multiplication result
    void collect_MW(ModuleBase::matrix& MecMulP, const ModuleBase::ComplexMatrix& mud, int nw, int isk);
    /// mulliken population = trace(dm*overlap)
    void cal_orbMulP();

  private:
    /******************************************************************
     * private variables
     *******************************************************************/
    Output_Sk<TK>* output_sk_ = nullptr;
    Output_DMK<TK>* output_dmk_ = nullptr;
    Parallel_Orbitals* ParaV_ = nullptr;
    CellIndex* cell_index_ = nullptr;
    const std::vector<int>& isk_;
    int nspin_;
    ModuleBase::matrix orbMulP_;
};

} // namespace ModuleIO

#endif // OUTPUT_MULLIKEN_H

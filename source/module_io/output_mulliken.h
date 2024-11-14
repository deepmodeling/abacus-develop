#ifndef OUTPUT_MULLIKEN_H
#define OUTPUT_MULLIKEN_H
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/cell_index.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_io/output_dmk.h"
#include "module_io/output_sk.h"

#include <map>
#include <vector>

namespace ModuleIO
{

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

template <typename TK>
void cal_mag(Parallel_Orbitals* pv,
             hamilt::Hamilt<TK>* p_ham,
             K_Vectors& kv,
             elecstate::ElecState* pelec,
             UnitCell& ucell,
             const int istep,
             const bool print)
{
    auto cell_index
        = CellIndex(ucell.get_atomLabels(), ucell.get_atomCounts(), ucell.get_lnchiCounts(), PARAM.inp.nspin);
    auto out_sk = ModuleIO::Output_Sk<TK>(p_ham, pv, PARAM.inp.nspin, kv.get_nks());
    auto out_dmk = ModuleIO::Output_DMK<TK>(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(pelec)->get_DM(),
                                            pv,
                                            PARAM.inp.nspin,
                                            kv.get_nks());
    auto mulp = ModuleIO::Output_Mulliken<TK>(&(out_sk), &(out_dmk), pv, &cell_index, kv.isk, PARAM.inp.nspin);
    auto atom_chg = mulp.get_atom_chg();
    /// used in updating mag info in STRU file
    ucell.atom_mulliken = mulp.get_atom_mulliken(atom_chg);
    if (print && GlobalV::MY_RANK == 0)
    {
        /// write the Orbital file
        cell_index.write_orb_info(PARAM.globalv.global_out_dir);
        /// write mulliken.txt
        mulp.write(istep, PARAM.globalv.global_out_dir);
        /// write atomic mag info in running log file
        mulp.print_atom_mag(atom_chg, GlobalV::ofs_running);
    }
}

} // namespace ModuleIO

#endif // OUTPUT_MULLIKEN_H

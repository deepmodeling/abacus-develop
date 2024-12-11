#ifndef OUTPUT_MULLIKEN_H
#define OUTPUT_MULLIKEN_H
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/cell_index.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_io/output_dmk.h"
#include "module_io/output_sk.h"
#include "module_base/formatter.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dspin_lcao.h"

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
             const TwoCenterBundle& two_center_bundle,
             const LCAO_Orbitals& orb,
             UnitCell& ucell,
             Grid_Driver& gd,
             const int istep,
             const bool print)
{
    // 1) calculate and output Mulliken population charges and magnetic moments
    if (PARAM.inp.out_mul)
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
    // 2) calculate and output the magnetizations of each atom with projection method
    if (PARAM.inp.onsite_radius > 0)
    {
        std::vector<std::vector<double>> atom_mag(ucell.nat, std::vector<double>(PARAM.inp.nspin, 0.0));
        std::vector<ModuleBase::Vector3<int>> constrain(ucell.nat, ModuleBase::Vector3<int>(1, 1, 1));
        const hamilt::HContainer<double>* dmr
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(pelec)->get_DM()->get_DMR_pointer(1);
        std::vector<double> moments;
        std::vector<double> mag_x(ucell.nat, 0.0);
        std::vector<double> mag_y(ucell.nat, 0.0);
        std::vector<double> mag_z(ucell.nat, 0.0);
        auto atomLabels = ucell.get_atomLabels();
        if(PARAM.inp.nspin == 2)
        {
            auto sc_lambda
                = new hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, double>>(nullptr,
                                                                          kv.kvec_d,
                                                                          nullptr,
                                                                          ucell,
                                                                          &gd,
                                                                          two_center_bundle.overlap_orb_onsite.get(),
                                                                          orb.cutoffs());
            dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(pelec)->get_DM()->switch_dmr(2);
            moments = sc_lambda->cal_moment(dmr, constrain);
            dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(pelec)->get_DM()->switch_dmr(0);
            delete sc_lambda;
            //const std::vector<std::string> title = {"Total Magnetism (uB)", ""};
            //const std::vector<std::string> fmts = {"%-26s", "%20.10f"};
            //FmtTable table(title, ucell.nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
            for(int iat=0;iat<ucell.nat;iat++)
            {
                atom_mag[iat][0] = 0.0;
                atom_mag[iat][1] = moments[iat];
            //    mag_z[iat] = moments[iat];
            }
            //table << atomLabels << mag_z;
            //GlobalV::ofs_running << table.str() << std::endl;
        }
        else if(PARAM.inp.nspin == 4)
        {
            auto sc_lambda = new hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(
                nullptr,
                kv.kvec_d,
                nullptr,
                ucell,
                &gd,
                two_center_bundle.overlap_orb_onsite.get(),
                orb.cutoffs());
            moments = sc_lambda->cal_moment(dmr, constrain);
            delete sc_lambda;
            //const std::vector<std::string> title = {"Total Magnetism (uB)", "", "", ""};
            //const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
            //FmtTable table(title, ucell.nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
            for(int iat=0;iat<ucell.nat;iat++)
            {
                atom_mag[iat][0] = 0.0;
                atom_mag[iat][1] = moments[iat*3];
                atom_mag[iat][2] = moments[iat*3+1];
                atom_mag[iat][3] = moments[iat*3+2];
                //mag_x[iat] = moments[iat*3];
                //mag_y[iat] = moments[iat*3+1];
                //mag_z[iat] = moments[iat*3+2];
            }
            //table << atomLabels << mag_x << mag_y << mag_z;
            //GlobalV::ofs_running << table.str() << std::endl;
        }
        ucell.atom_mulliken = atom_mag;
    }
}

} // namespace ModuleIO

#endif // OUTPUT_MULLIKEN_H

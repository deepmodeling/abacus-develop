#ifndef CAL_MAG_H
#define CAL_MAG_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/cell_index.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/hamilt.h"
#include "source_io/module_mulliken/output_dmk.h"
#include "source_io/module_mulliken/output_mulliken.h"
#include "source_io/module_mulliken/output_sk.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"

#include <vector>

namespace ModuleIO
{

template <typename TK>
void cal_mag(Parallel_Orbitals* pv,
             hamilt::Hamilt<TK>* p_ham,
             K_Vectors& kv,
             elecstate::DensityMatrix<TK,double>* dm,
             const TwoCenterBundle& two_center_bundle,
             const LCAO_Orbitals& orb,
             UnitCell& ucell,
             const Grid_Driver& gd,
             const int istep,
             const bool print)
{
    if (PARAM.inp.out_mul)
    {
        auto cell_index
            = CellIndex(ucell.get_atomLabels(),
			ucell.get_atomCounts(), ucell.get_lnchiCounts(), PARAM.inp.nspin);
        auto out_s_k = ModuleIO::Output_Sk<TK>(p_ham, pv, PARAM.inp.nspin, kv.get_nks());
        auto out_dm_k = ModuleIO::Output_DMK<TK>(dm, pv, PARAM.inp.nspin, kv.get_nks());

        auto mulp = ModuleIO::Output_Mulliken<TK>(&(out_s_k),
			&(out_dm_k), pv, &cell_index, kv.isk, PARAM.inp.nspin);
        auto atom_chg = mulp.get_atom_chg();
        ucell.atom_mulliken = mulp.get_atom_mulliken(atom_chg);
        if (print && GlobalV::MY_RANK == 0)
        {
            cell_index.write_orb_info(PARAM.globalv.global_out_dir);
            mulp.write(istep, PARAM.globalv.global_out_dir);
            mulp.print_atom_mag(atom_chg, GlobalV::ofs_running);
        }
    }
    if (PARAM.inp.onsite_radius > 0)
    {
        std::vector<std::vector<double>> atom_mag(ucell.nat, std::vector<double>(PARAM.inp.nspin, 0.0));
        std::vector<ModuleBase::Vector3<int>> constrain(ucell.nat, ModuleBase::Vector3<int>(1, 1, 1));
        const hamilt::HContainer<double>* dmr = dm->get_DMR_pointer(1);
        std::vector<double> moments;
        std::vector<double> mag_x(ucell.nat, 0.0);
        std::vector<double> mag_y(ucell.nat, 0.0);
        std::vector<double> mag_z(ucell.nat, 0.0);
        auto atomLabels = ucell.get_atomLabels();

        if(PARAM.inp.nspin == 2)
        {
            auto sc_lambda = new hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, double>>(nullptr,
		kv.kvec_d,
		dynamic_cast<hamilt::HamiltLCAO<TK, double>*>(p_ham)->getHR(),
		ucell,
		&gd,
		two_center_bundle.overlap_orb_onsite.get(),
		orb.cutoffs());

	    dm->switch_dmr(2);
	    moments = sc_lambda->cal_moment(dmr, constrain);
	    dm->switch_dmr(0);

	    delete sc_lambda;

            for(int iat=0;iat<ucell.nat;iat++)
            {
                atom_mag[iat][0] = 0.0;
                atom_mag[iat][1] = moments[iat];
            }
        }
        else if(PARAM.inp.nspin == 4)
        {
            auto sc_lambda = new hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(
                nullptr,
                kv.kvec_d,
                dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->getHR(),
                ucell,
                &gd,
                two_center_bundle.overlap_orb_onsite.get(),
                orb.cutoffs());
            moments = sc_lambda->cal_moment(dmr, constrain);
            delete sc_lambda;

            for(int iat=0;iat<ucell.nat;iat++)
            {
                atom_mag[iat][0] = 0.0;
                atom_mag[iat][1] = moments[iat*3];
                atom_mag[iat][2] = moments[iat*3+1];
                atom_mag[iat][3] = moments[iat*3+2];
            }
        }
        ucell.atom_mulliken = atom_mag;
    }
}

} // namespace ModuleIO

#endif // CAL_MAG_H

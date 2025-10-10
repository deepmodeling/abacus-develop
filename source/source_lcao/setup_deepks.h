#ifndef SETUP_DEEPKS_H
#define SETUP_DEEPKS_H

#include "source_cell/unitcell.h" // use unitcell
#include "source_io/module_parameter/input_parameter.h" // Input_para
#include "source_basis/module_ao/parallel_orbitals.h" // parallel orbitals
#include "source_basis/module_ao/ORB_read.h" // orb
#include "source_basis/module_nao/two_center_integrator.h" // overlap_orb_alpha 
#include "source_cell/module_neighbor/sltk_grid_driver.h" // grid driver

#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h" // deepks
#endif


template <typename TK>
class Setup_DeePKS
{
    public:

    Setup_DeePKS();
    ~Setup_DeePKS();

#ifdef __MLALGO
    LCAO_Deepks<TK> ld;
#endif

    std::string dpks_out_type;

	void before_runner(
			const UnitCell &ucell, // unitcell
			const int nks, // k points
            const LCAO_Orbitals &orb, // orbital info
			Parallel_Orbitals &pv, // parallel orbitals
			const Input_para &inp);

    void build_overlap(
		const UnitCell &ucell,
		const LCAO_Orbitals &orb,
		const Parallel_Orbitals &pv,
		const Grid_Driver &gd,
        TwoCenterIntegrator &overlap_orb_alpha,
		const Input_para &inp);

    void write_forces(
		const ModuleBase::matrix &fcs,
		const ModuleBase::matrix &fvnl_dalpha,
		const Input_para &inp);

    void write_stress(
		const ModuleBase::matrix &scs,
		const ModuleBase::matrix &svnl_dalpha,
		const double &omega,
		const Input_para &inp);

};


#endif

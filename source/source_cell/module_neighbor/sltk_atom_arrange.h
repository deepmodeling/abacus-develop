#ifndef ATOM_ARRANGE_H
#define ATOM_ARRANGE_H

#include "sltk_grid.h"
#include "sltk_grid_driver.h"


class atom_arrange
{
public:

	atom_arrange();
	~atom_arrange();
	
    /**
     * Build a replicated neighbor Grid.
     *
     * search_radius_bohr is the physical cutoff. skin_bohr extends only the
     * candidate-list build radius; returned neighbors are still filtered by
     * the physical cutoff. A zero skin preserves full-rebuild behavior.
     */
	static void search(
		const bool flag,
		std::ofstream &ofs,
		Grid_Driver &grid_d, 
		const UnitCell &ucell, 
		const double& search_radius_bohr, 
		const int &test_atom_in,
		const bool test_only = false,
        const double skin_bohr = 0.0);

    /**
     * Build a local-center MPI Grid using the same Bohr-to-lat0 conversion as
     * search(). search_radius_bohr and skin_bohr have the same meaning as in
     * search(). This grid is intended for consumers that understand spatial
     * ownership; nonlocal center queries are rejected.
     */
    static void search_mpi(const bool pbc_flag,
                           std::ofstream& ofs,
                           Grid_Driver& grid_d,
                           const UnitCell& ucell,
                           const double& search_radius_bohr,
                           ModuleNeighbor::NeighborMpiComm communicator,
                           ModuleNeighbor::MpiGhostExchangeStats* stats = nullptr,
                           const double skin_bohr = 0.0);

	//caoyu modify 2021-05-24
	static double set_sr_NL(
		std::ofstream &ofs_in,
		const std::string &output_level,
		const double& rcutmax_Phi, 
		const double& rcutmax_Beta, 
		const bool gamma_only_local);
};

#endif

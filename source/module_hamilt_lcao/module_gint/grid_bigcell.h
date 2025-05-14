#ifndef GRID_BIGCELL_H
#define GRID_BIGCELL_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"
#include "grid_meshcell.h"

class Grid_BigCell: public Grid_MeshCell
{
	public:
		Grid_BigCell();
		~Grid_BigCell();
		// number atoms and type.
		int nat=0;
		// save the relative cartesian position
		// to bigcell of each atom.
        std::vector<std::vector<double>> tau_in_bigcell;

        /// move operator for the next ESolver to directly use its infomation
        Grid_BigCell& operator=(Grid_BigCell&& rhs) = default;

      protected:
        // get the max radius of all orbitals
		// which will use to generate grid expansion,
		// and the meshball.
		double orbital_rmax;
		
		// the added number of bigcelli each direction.
		int dxe;
		int dye;
		int dze;

		// expansion grid dimension.
		int nxe;
		int nye;
		int nze;
		int nxyze;
		
        std::vector<int> index_atom;

        // save the position of base vector of bigcell.
        std::vector<double> bigcell_vec1;
        std::vector<double> bigcell_vec2;
        std::vector<double> bigcell_vec3;

		ModuleBase::Matrix3 bigcell_latvec0;
		ModuleBase::Matrix3 bigcell_GT;
		
		//---------------------------------
		void grid_expansion_index(bool f2normal, int *target)const;
		//---------------------------------
		void init_big_latvec(const UnitCell &ucell);
		//---------------------------------
		void init_tau_in_bigcell(const UnitCell& ucell);
		//---------------------------------
		void init_grid_expansion(const UnitCell& ucell,double* rcut);
};
#endif

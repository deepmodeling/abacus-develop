#ifndef GRID_MESHBALL_H
#define GRID_MESHBALL_H

#include "grid_bigcell.h"

class Grid_MeshBall : public Grid_BigCell
{
	public:
		Grid_MeshBall();
		~Grid_MeshBall();
		// cartesian coordinates of meshball.
        std::vector<std::vector<double>> meshball_positions;

        /// move operator for the next ESolver to directly use its infomation
        Grid_MeshBall& operator=(Grid_MeshBall&& rhs) = default;

      protected:
		// number of meshcells in meshball.
		int meshball_ncells=0;
		// used in index2normal
		std::vector<int> index_ball;
		// search each meshcell of this meshball.
		void init_meshball(void);

	private:
		// init the meshball radius.
		double meshball_radius=0.0;
		// Handle as a truncation function.
		double deal_with_atom_spillage(const double* pos);
	
};
#endif

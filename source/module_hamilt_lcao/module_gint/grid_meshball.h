#ifndef GRID_MESHBALL_H
#define GRID_MESHBALL_H

#include "grid_bigcell.h"

class Grid_MeshBall : public Grid_BigCell
{
	public:
	// number of meshcells in meshball.
	int meshball_ncells;

	// cartesian coordinates of meshball.
	std::vector<std::vector<double>> meshball_positions;	
	
	protected:

	Grid_MeshBall();
	~Grid_MeshBall();	

	// used in index2normal
	std::vector<int> index_ball;

	// init the meshball radius,
	// search each meshcell of this meshball.
	void init_meshball(void);

	private:

	double deal_with_atom_spillage(const double* pos);

	bool flag_mp;

	double meshball_radius;
	
};
#endif

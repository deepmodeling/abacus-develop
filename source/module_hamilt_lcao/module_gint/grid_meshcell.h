#ifndef GRID_MESHCELL_H
#define GRID_MESHCELL_H
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"
#include "grid_meshk.h"
#include "module_cell/unitcell.h"
class Grid_MeshCell: public Grid_MeshK
{
	public:
	Grid_MeshCell();
	~Grid_MeshCell();
	
	int ncx,ncy,ncz,ncxyz;
	int bx=1,by=1,bz=1,bxyz=1;
	int nbx,nby,nbz,nbxyz;
	int nbxx;
	int nbzp_start,nbzp;
	// save the position of each meshcell.
	std::vector<std::vector<double>> meshcell_pos;
	ModuleBase::Matrix3 meshcell_latvec0;
	ModuleBase::Matrix3 meshcell_GT;
	
	protected:

	std::vector<double> meshcell_vec1;
	std::vector<double> meshcell_vec2;
	std::vector<double> meshcell_vec3;

	void set_grid_dim(
			const int &ncx_in,
			const int &ncy_in,
			const int &ncz_in,
			const int &bx_in,
			const int &by_in,
			const int &bz_in,
			const int &nbx_in,
			const int &nby_in,
			const int &nbz_in,
			const int &nbxx_in,
			const int &nbzp_start_in,
			const int &nbzp_in);

	void init_latvec(const UnitCell &ucell);
    void init_meshcell_pos(void);

};

#endif

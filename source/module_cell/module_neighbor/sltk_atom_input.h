#ifndef ATOM_INPUT_H
#define ATOM_INPUT_H

#include "sltk_atom.h"
#include "module_cell/unitcell.h"

class Atom_input
{
public:
//==========================================================
// Constructors and destructor
//==========================================================
	Atom_input
	(
		std::ofstream &ofs_in,
		const UnitCell &ucell,
		const int amount = 0,	//number of atoms
	    const int ntype = 0,	//number of atom_types
	    const bool boundary = true,	// 1 : periodic ocndition
		const double radius_in = 0, // searching radius
		const int &test_atom_in = 0	//caoyu reconst 2021-05-24
	);
	~Atom_input();
//==========================================================
// Manipulators
//==========================================================
	void set_FAtom(const UnitCell &ucell, FAtom& a)const ;

	double vec1[3];
	double vec2[3];
	double vec3[3];

public:
	bool getExpandFlag()const {return expand_flag;}

	int getAmount() const
	{
		if (!expand_flag) { return d_amount;
		} else { return d_amount_expand;
}
	}

	int getBoundary()const { return periodic_boundary;}

	double getLatNow() const { return lat_now;}

	double getRadius() const {return radius;}

//==========================================================
//
//==========================================================
	double getCellXLength() const
	{
		if (!expand_flag) { return radius;
		} else { return 1;
}
	}

	double getCellYLength() const
	{
		if (!expand_flag) { return radius;
		} else { return 1;
}
	}

	double getCellZLength() const
	{
		if (!expand_flag) { return radius;
		} else { return 1;
}
	}

//==========================================================
//
//==========================================================
	double Clength0() const { return (glayerX + glayerX_minus) * clength0;}
	double Clength1() const { return (glayerY + glayerY_minus) * clength1;}
	double Clength2() const { return (glayerZ + glayerZ_minus) * clength2;}

//==========================================================
//
//==========================================================
	double minX() const
	{
		if (!expand_flag) { return x_min;
		} else { return (double)(-glayerX_minus);
}
	}

	double minY() const
	{
		if (!expand_flag) { return y_min;
		} else { return (double)(-glayerY_minus);
}
	}

	double minZ() const
	{
		if (!expand_flag) { return z_min;
		} else { return (double)(-glayerZ_minus);
}
	}

//==========================================================
//
//==========================================================
	int getCellX() const { return cell_nx; }

	int getCellY() const { return cell_ny; }

	int getCellZ() const { return cell_nz; }

//==========================================================
//
//==========================================================
	int getGrid_layerX() const { return glayerX;}

	int getGrid_layerX_minus() const { return glayerX_minus;}

	int getGrid_layerY() const { return glayerY;}

	int getGrid_layerY_minus() const { return glayerY_minus;}

	int getGrid_layerZ() const { return glayerZ;}

	int getGrid_layerZ_minus() const { return glayerZ_minus;}

	FAtom getFakeAtom(const int index) const { return fake_atoms[index];}

private:
	int test_atom_input;	//caoyu reconst 2021-05-24
	int d_amount;//number of atoms.
	int d_amount_expand;
	bool periodic_boundary;

	double lat_now;
	double radius;

	double clength0;
	double clength1;
	double clength2;

	double x_min;
	double y_min;
	double z_min;
	double x_max;
	double y_max;
	double z_max;
//==========================================================
// MEMBRE FUNCTION :
// NAME : Check_Expand_Condition
//==========================================================
	void Check_Expand_Condition(const UnitCell& ucell);
	bool expand_flag;
	int glayerX;
	int glayerX_minus;
	int glayerY;
	int glayerY_minus;
	int glayerZ;
	int glayerZ_minus;
//==========================================================
// MEMBRE FUNCTION :
// NAME : Expand_Grid
//==========================================================
	void Expand_Grid(const UnitCell& ucell, const int ntype);

	std::vector<FAtom> fake_atoms;

	double x_min_expand;
	double y_min_expand;
	double z_min_expand;
	double x_max_expand;
	double y_max_expand;
	double z_max_expand;
//==========================================================
// MEMBRE FUNCTION :
// NAME : Expand_Grid
//==========================================================
	void calculate_cells();
	int cell_nx;
	int cell_ny;
	int cell_nz;
//==========================================================
// MEMBRE FUNCTION :
// NAME : Load_atom
//==========================================================

	mutable int d_current;
	mutable int type;
	mutable int natom;
};

#endif

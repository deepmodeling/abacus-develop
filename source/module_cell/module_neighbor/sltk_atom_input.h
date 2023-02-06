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
	    const bool boundary = 1,	// 1 : periodic ocndition
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
	bool getExpandFlag(void)const {return expand_flag;}

	int getAmount() const
	{
		if (!expand_flag) return d_amount;
		else return d_amount_expand;
	}

	int getBoundary()const { return periodic_boundary;}

	double getLatNow(void) const { return lat_now;}

	double getRadius() const {return radius;}

//==========================================================
//
//==========================================================
	double getCellXLength(void) const
	{
		if (!expand_flag) return radius;
		else return 1;
	}

	double getCellYLength(void) const
	{
		if (!expand_flag) return radius;
		else return 1;
	}

	double getCellZLength(void) const
	{
		if (!expand_flag) return radius;
		else return 1;
	}

//==========================================================
//
//==========================================================
	double Clength0(void) const { return (glayerX + glayerX_minus) * clength0;}
	double Clength1(void) const { return (glayerY + glayerY_minus) * clength1;}
	double Clength2(void) const { return (glayerZ + glayerZ_minus) * clength2;}

//==========================================================
//
//==========================================================
	double minX(void) const
	{
		if (!expand_flag) return x_min;
		else return (double)(-glayerX_minus);
	}

	double minY(void) const
	{
		if (!expand_flag) return y_min;
		else return (double)(-glayerY_minus);
	}

	double minZ(void) const
	{
		if (!expand_flag) return z_min;
		else return (double)(-glayerZ_minus);
	}

//==========================================================
//
//==========================================================
	int getCellX(void) const { return cell_nx; }

	int getCellY(void) const { return cell_ny; }

	int getCellZ(void) const { return cell_nz; }

//==========================================================
//
//==========================================================
	int getGrid_layerX(void) const { return glayerX;}

	int getGrid_layerX_minus(void) const { return glayerX_minus;}

	int getGrid_layerY(void) const { return glayerY;}

	int getGrid_layerY_minus(void) const { return glayerY_minus;}

	int getGrid_layerZ(void) const { return glayerZ;}

	int getGrid_layerZ_minus(void) const { return glayerZ_minus;}

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
	double* store_x;
	double* store_y;
	double* store_z;
	int* store_cell_x;
	int* store_cell_y;
	int* store_cell_z;
	int* store_type;
	int* store_natom;
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
	void calculate_cells(void);
	int cell_nx;
	int cell_ny;
	int cell_nz;
//==========================================================
// MEMBRE FUNCTION :
// NAME : Load_atom
//==========================================================
	void Load_atom(const UnitCell& ucell)const;
	mutable int d_current;
	mutable double x;
	mutable double y;
	mutable double z;
	mutable int type;
	mutable int natom;
};

#endif

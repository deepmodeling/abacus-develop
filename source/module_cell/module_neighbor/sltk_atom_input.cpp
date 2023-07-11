#include "sltk_atom_input.h"
#include "sltk_grid.h"
#include "module_base/memory.h"

//==========================================================
// define constructor and deconstructor
//==========================================================
Atom_input::Atom_input
(
	std::ofstream &ofs_in,
	const UnitCell &ucell,
	const int amount,
	const int ntype,
    const bool boundary_in,
	const double radius_in,
	const int &test_atom_in
)
		:d_amount(amount),
		d_amount_expand(amount),
		periodic_boundary(boundary_in),
		radius(radius_in),
		expand_flag(false),
		glayerX(1),
		glayerX_minus(0),
		glayerY(1),
		glayerY_minus(0),
		glayerZ(1),
		glayerZ_minus(0),
		x_min_expand(0),
		y_min_expand(0),
		z_min_expand(0),
		x_max_expand(0),
		y_max_expand(0),
		z_max_expand(0),
		d_current(0),
		test_atom_input(test_atom_in),

//----------------------------------------------------------
// WARNING :
// Please be very very careful!
// Here type = 0 ,and natom = -1 is a initial value,
// don't change it !!
//----------------------------------------------------------
		type(0),
		natom(-1)
{
	ModuleBase::TITLE("Atom_input", "Atom_input");

	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in, "ntype", ntype);
	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in, "Amount(atom number)", amount);
	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in, "Periodic_boundary", periodic_boundary);

//----------------------------------------------------------
// EXPLAIN : check searching raidus
//----------------------------------------------------------
	if(test_atom_input)ModuleBase::GlobalFunc::OUT(ofs_in, "Searching radius(lat0)", radius);

	if (radius < 0)
	{
		ModuleBase::WARNING_QUIT("atom_arrange::init", " search radius < 0,forbidden");
	}

	//===========================
	// Set three lattice vectors
	//===========================
	vec1[0] = ucell.latvec.e11;
	vec1[1] = ucell.latvec.e12;
	vec1[2] = ucell.latvec.e13;

	vec2[0] = ucell.latvec.e21;
	vec2[1] = ucell.latvec.e22;
	vec2[2] = ucell.latvec.e23;

	vec3[0] = ucell.latvec.e31;
	vec3[1] = ucell.latvec.e32;
	vec3[2] = ucell.latvec.e33;

	if(GlobalV::test_grid)
	{
		ofs_in << " Output lattice vectors now (unit:lat0):" << std::endl;
		ofs_in << " " << std::setw(5) << "Vec1" 
			<< std::setw(10) << vec1[0]
			<< std::setw(10) << vec1[1]
			<< std::setw(10) << vec1[2] << std::endl;
		ofs_in << " " << std::setw(5) << "Vec2" 
			<< std::setw(10) << vec2[0]
			<< std::setw(10) << vec2[1]
			<< std::setw(10) << vec2[2] << std::endl;
		ofs_in << " " << std::setw(5) << "Vec3" 
			<< std::setw(10) << vec3[0]
			<< std::setw(10) << vec3[1]
			<< std::setw(10) << vec3[2];
		ofs_in << std::endl;
	}

	//=============================================
	// Use Lattice vectors to generate cell_length
	//=============================================
	clength0 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]) ;
	clength1 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]) ;
	clength2 = sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]) ;

	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in,"CellLength(unit: lat0)",clength0,clength1,clength2);
	//==============================
	// set lattice constant
	//==============================
	lat_now = ucell.lat0;

	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in, "lat0_now (Bohr)", lat_now);

	// random selection, in order to estimate again.
	this->x_min = ucell.atoms[0].tau[0].x;
	this->y_min = ucell.atoms[0].tau[0].y;
	this->z_min = ucell.atoms[0].tau[0].z;
	this->x_max = ucell.atoms[0].tau[0].x;
	this->y_max = ucell.atoms[0].tau[0].y;
	this->z_max = ucell.atoms[0].tau[0].z;

	//calculate min & max value
	for (int i = 0;i < ntype;i++)
	{
		for (int j = 0;j < ucell.atoms[i].na;j++)
		{
			if (ucell.atoms[i].tau[j].x < x_min) this->x_min = ucell.atoms[i].tau[j].x;
			if (ucell.atoms[i].tau[j].y < y_min) this->y_min = ucell.atoms[i].tau[j].y;
			if (ucell.atoms[i].tau[j].z < z_min) this->z_min = ucell.atoms[i].tau[j].z;
			if (ucell.atoms[i].tau[j].x > x_max) this->x_max = ucell.atoms[i].tau[j].x;
			if (ucell.atoms[i].tau[j].y > y_max) this->y_max = ucell.atoms[i].tau[j].y;
			if (ucell.atoms[i].tau[j].z > z_max) this->z_max = ucell.atoms[i].tau[j].z;
		}
	}

	if(test_atom_input)
	{
		ofs_in << " Find the coordinate range of the input atom(unit:lat0)." << std::endl;
	}
	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in,"min_tau", x_min, y_min, z_min);
	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in,"max_tau", x_max, y_max, z_max);

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : Check_Expand_Condition(check if swe need to
// expand grid,and generate 6 MEMBER VARIABLE(number of
// layers for 6 dimension)
// initial value for "glayerX,Y,Z" : 1
// (if > 2 ,expand flag = 1)
// initial value for "glayerX,Y,Z_minus" : 0
// ( if > 1 ,expand flag = 1)
//----------------------------------------------------------

	this->Check_Expand_Condition(ucell);

	// for check
	//glayerX+=2;
	//glayerY+=2;
	//glayerZ+=2;
	//glayerX_minus-=2;
	//glayerY_minus-=2;
	//glayerZ_minus-=2;

	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in,"glayer+",glayerX,glayerY,glayerZ);
	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in,"glayer-",glayerX_minus,glayerY_minus,glayerZ_minus);

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : Expand_Grid (after expand grid,we set another
// 6 MEMBER VARIABLES : x(y,z)_min(max)_expand,record
// the peak cartesian value among all the cells
//----------------------------------------------------------

	if (this->expand_flag)
	{
		this->Expand_Grid(ucell, ntype);
	}

	if(GlobalV::test_grid) ModuleBase::GlobalFunc::OUT(ofs_in, "expand_flag", expand_flag);

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : calculate_cells
// Calculate how many cells we need in each direction.
//----------------------------------------------------------
	this->calculate_cells();
	if(test_atom_input) ModuleBase::GlobalFunc::OUT(ofs_in, "CellDim", cell_nx, cell_ny, cell_nz);
	return;
}

Atom_input::~Atom_input()
{
	if (expand_flag)
	{
		delete [] store_x;
		delete [] store_y;
		delete [] store_z;
		delete [] store_type;
		delete [] store_natom;
	}
}

//============================================
// !!!! May still have bug, be very careful!!
// should use the same algorithm to generate
// dxe, dye, dze in grid_meshcell.cpp.
//============================================
void Atom_input::Check_Expand_Condition(const UnitCell &ucell)
{
//	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

	if (!periodic_boundary) return;

	// mohan update 2011-04-14
	// d stands for direct coordinates.
	double dmaxX=ucell.atoms[0].taud[0].x;
	double dmaxY=ucell.atoms[0].taud[0].y;
	double dmaxZ=ucell.atoms[0].taud[0].z;
	
	double dminX=ucell.atoms[0].taud[0].x;
	double dminY=ucell.atoms[0].taud[0].y;
	double dminZ=ucell.atoms[0].taud[0].z;

	//calculate min & max value
	for (int i = 0;i < ucell.ntype;i++)
	{
		for (int j = 0;j < ucell.atoms[i].na;j++)
		{
			dminX = std::min( dminX, ucell.atoms[i].taud[j].x );
			dminY = std::min( dminY, ucell.atoms[i].taud[j].y );
			dminZ = std::min( dminZ, ucell.atoms[i].taud[j].z );

			dmaxX = std::max( dmaxX, ucell.atoms[i].taud[j].x );
			dmaxY = std::max( dmaxY, ucell.atoms[i].taud[j].y );
			dmaxZ = std::max( dmaxZ, ucell.atoms[i].taud[j].z );
		}
	}

	if(dminX<0.0)
	{
		std::cout << " dminX=" << dminX << std::endl;
		ModuleBase::WARNING_QUIT("Atom_input::Check_Expand_Condition","dminX<0.0");
	}
	if(dminY<0.0)
	{
		std::cout << " dminY=" << dminY << std::endl;
		ModuleBase::WARNING_QUIT("Atom_input::Check_Expand_Condition","dminY<0.0");
	}
	if(dminZ<0.0)
	{
		std::cout << " dminZ=" << dminZ << std::endl;
		ModuleBase::WARNING_QUIT("Atom_input::Check_Expand_Condition","dminZ<0.0");
	}


	if(test_atom_input)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Radius",radius);

/*2016-07-19, LiuXh
	// the unit of extent_1DX,Y,Z is lat0.
	// means still how far can be included now.
	double extent_1DX = glayerX * clength0 - dmaxX;
	while (radius > extent_1DX)
	{
		glayerX++;
		extent_1DX = glayerX * clength0 - dmaxX;
	}
	double extent_1DY = glayerY * clength1 - dmaxY;
	while (radius > extent_1DY)
	{
		glayerY++;
		extent_1DY = glayerY * clength1 - dmaxY;
	}
	double extent_1DZ = glayerZ * clength2 - dmaxZ;
	while (radius > extent_1DZ)
	{
		glayerZ++;
		extent_1DZ = glayerZ * clength2 - dmaxZ;
	}

	// in case the cell is not retangle.
	// mohan added 2009-10-23
	// if this is not added, it's a serious bug.
	glayerX++;
	glayerY++;
	glayerZ++;
	if(test_atom_input)
	{
		GlobalV::ofs_running << " Extend distance from the (maxX,maxY,maxZ) direct position in this unitcell: " << std::endl;
	}
	
	if(test_atom_input)OUT(GlobalV::ofs_running,"ExtentDim+",extent_1DX,extent_1DY,extent_1DZ);

	double extent_1DX_minus = glayerX_minus * clength0 + dminX;
	while (radius > extent_1DX_minus)
	{
		glayerX_minus++;
		extent_1DX_minus = glayerX_minus * clength0 + dminX;
	}
	double extent_1DY_minus = glayerY_minus * clength1 + dminY;
	while (radius > extent_1DY_minus)
	{
		glayerY_minus++;
		extent_1DY_minus = glayerY_minus * clength1 + dminY;
	}
	double extent_1DZ_minus = glayerZ_minus * clength2 + dminZ;
	while (radius > extent_1DZ_minus)
	{
		glayerZ_minus++;
		extent_1DZ_minus = glayerZ_minus * clength2 + dminZ;
	}

	// in case the cell is not retangle.
	// mohan added 2009-10-23
	// if this is not added, it's a serious bug.
	glayerX_minus++;
	glayerY_minus++;
	glayerZ_minus++;

	//glayerX_minus++;
	//glayerY_minus++;
	//glayerZ_minus++;
2016-07-19, LiuXh*/
	//Begin, 2016-07-19, LiuXh
	double a23_1 = ucell.latvec.e22*ucell.latvec.e33 - ucell.latvec.e23*ucell.latvec.e32;
	double a23_2 = ucell.latvec.e21*ucell.latvec.e33 - ucell.latvec.e23*ucell.latvec.e31;
	double a23_3 = ucell.latvec.e21*ucell.latvec.e32 - ucell.latvec.e22*ucell.latvec.e31;
	double a23_norm = sqrt(a23_1*a23_1 + a23_2*a23_2 + a23_3*a23_3);
	double extend_v = a23_norm * radius;
	double extend_d1 = extend_v/ucell.omega*ucell.lat0*ucell.lat0*ucell.lat0;
	int extend_d11 = static_cast<int>(extend_d1);
	//2016-09-05, LiuXh
	if(extend_d1 - extend_d11 > 0.0) extend_d11 += 1;

	double a31_1 = ucell.latvec.e32*ucell.latvec.e13 - ucell.latvec.e33*ucell.latvec.e12;
	double a31_2 = ucell.latvec.e31*ucell.latvec.e13 - ucell.latvec.e33*ucell.latvec.e11;
	double a31_3 = ucell.latvec.e31*ucell.latvec.e12 - ucell.latvec.e32*ucell.latvec.e11;
	double a31_norm = sqrt(a31_1*a31_1 + a31_2*a31_2 + a31_3*a31_3);
	double extend_d2 = a31_norm*radius/ucell.omega*ucell.lat0*ucell.lat0*ucell.lat0;
	int extend_d22 = static_cast<int>(extend_d2);
	//2016-09-05, LiuXh
	if(extend_d2 - extend_d22 > 0.0) extend_d22 += 1;

	double a12_1 = ucell.latvec.e12*ucell.latvec.e23 - ucell.latvec.e13*ucell.latvec.e22;
	double a12_2 = ucell.latvec.e11*ucell.latvec.e23 - ucell.latvec.e13*ucell.latvec.e21;
	double a12_3 = ucell.latvec.e11*ucell.latvec.e22 - ucell.latvec.e12*ucell.latvec.e21;
	double a12_norm = sqrt(a12_1*a12_1 + a12_2*a12_2 + a12_3*a12_3);
	double extend_d3 = a12_norm * radius/ucell.omega*ucell.lat0*ucell.lat0*ucell.lat0;
	int extend_d33 = static_cast<int>(extend_d3);
	//2016-09-05, LiuXh
	if(extend_d3 - extend_d33 > 0.0) extend_d33 += 1;

	glayerX = extend_d11 +1;
	glayerY = extend_d22 +1;
	glayerZ = extend_d33 +1;
	//Begin, 2016-09-05, LiuXh
	//glayerX_minus = extend_d11 +1;
	//glayerY_minus = extend_d22 +1;
	//glayerZ_minus = extend_d33 +1;
	glayerX_minus = extend_d11;
	glayerY_minus = extend_d22;
	glayerZ_minus = extend_d33;
	//End, 2016-09-05, LiuXh

	if(glayerX==1) glayerX++;
	if(glayerY==1) glayerY++;
	if(glayerZ==1) glayerZ++;
	if(glayerX_minus==1) glayerX_minus++;
	if(glayerY_minus==1) glayerY_minus++;
	if(glayerZ_minus==1) glayerZ_minus++;
	//End, 2016-07-19, LiuXh
/*	
	if(test_atom_input)
	{
		GlobalV::ofs_running << " Extend distance from the (minX,minY,minZ) direct position in this unitcell: " << std::endl;
	}

	if(test_atom_input)OUT(GlobalV::ofs_running,"ExtentDim-",extent_1DX_minus,extent_1DY_minus,extent_1DZ_minus);
*/
//----------------------------------------------------------
// EXPLAIN : if extent don't satisfty the searching
// requiment, we must expand one more layer
//----------------------------------------------------------

	if (glayerX > 2 || glayerY > 2 || glayerZ > 2)
	{
		this->expand_flag = true;
	}
	else if (glayerX_minus > 1 || glayerX_minus > 1 || glayerX_minus > 1)
	{
		this->expand_flag = true;
	}
	else
	{
		this->expand_flag = false;
	}
	return;
}

void Atom_input::Expand_Grid(const UnitCell &ucell, const int ntype)
{
	ModuleBase::TITLE("Atom_input", "Expand_Grid");

	if(test_atom_input)
	{
		GlobalV::ofs_running << " Be careful of thie grid adjacent searching program!" << std::endl;
		GlobalV::ofs_running << " Here I would like to say some gudlines:" << std::endl;
		GlobalV::ofs_running << " You are using Expand_Grid now, which means now you treat" << std::endl;
		GlobalV::ofs_running << " your 'unitcell' as a Cell Class which defined in grid class" << std::endl;
		GlobalV::ofs_running << " This Cell is diffenent from the 'Not expand' cell." << std::endl;
		GlobalV::ofs_running << " In most cases, it may not be a cubic, so please do it more carefully." << std::endl;
		GlobalV::ofs_running << " Good luck! " << std::endl;
	}

	double *x_old = new double[d_amount];
	double *y_old = new double[d_amount];
	double *z_old = new double[d_amount];

	int *type_old = new int[d_amount];
	int *natom_old = new int[d_amount];

	int ia = 0;
	for (int i = 0;i < ntype;i++)
	{
		for (int j = 0;j < ucell.atoms[i].na;j++)
		{
			x_old[ia] = ucell.atoms[i].tau[j].x;
			y_old[ia] = ucell.atoms[i].tau[j].y;
			z_old[ia] = ucell.atoms[i].tau[j].z;
			type_old[ia] = i;
			natom_old[ia] = j;
			ia++;
		}
	}

	// how many copys we need now.
	const int gcopy =
	    (glayerX + glayerX_minus) *
	    (glayerY + glayerY_minus) *
	    (glayerZ + glayerZ_minus) ;

	if(test_atom_input)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Grid_copy_times",gcopy);

	this->d_amount_expand = d_amount * gcopy;

	if(test_atom_input)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Atom_number_now",d_amount_expand);

	// store new atom positions.
	this->store_x = new double[d_amount_expand];
	this->store_y = new double[d_amount_expand];
	this->store_z = new double[d_amount_expand];

	ModuleBase::Memory::record("SLTK::Epd_Atom",sizeof(double) * d_amount_expand*3);

	// store which grid the atom is in.
	store_cell_x = new int[d_amount_expand];
	store_cell_y = new int[d_amount_expand];
	store_cell_z = new int[d_amount_expand];

	this->store_type = new int[d_amount_expand];
	this->store_natom = new int[d_amount_expand];

	ModuleBase::Memory::record("SLTK::Epd_atom_info",sizeof(int) * d_amount_expand*5);

	int ia_all = 0;

	for (int ix = -glayerX_minus; ix < glayerX; ix++)
	{
		for (int iy = -glayerY_minus; iy < glayerY; iy++)
		{
			for (int iz = -glayerZ_minus; iz < glayerZ; iz++)
			{
				for (int ia = 0; ia < d_amount; ia++)
				{
					store_x[ia_all] = x_old[ia] + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
					store_y[ia_all] = y_old[ia] + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
					store_z[ia_all] = z_old[ia] + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
					store_type[ia_all] = type_old[ia];
					store_natom[ia_all] = natom_old[ia];

					// notice '0' is not the origin unitcell.
					store_cell_x[ia_all] = ix + glayerX_minus;
					store_cell_y[ia_all] = iy + glayerY_minus;
					store_cell_z[ia_all] = iz + glayerZ_minus;

					if (test_atom_input > 1)
					{
						if (d_amount_expand < 1000)
						{
							GlobalV::ofs_running << "\n" << std::setw(6) << ia_all
							<< std::setw(10) << x_old[ia]
							<< std::setw(10) << y_old[ia]
							<< std::setw(10) << z_old[ia]
							<< std::setw(10) << store_x[ia_all]
							<< std::setw(10) << store_y[ia_all]
							<< std::setw(10) << store_z[ia_all]
							<< std::setw(6) << store_cell_x[ia_all]
							<< std::setw(6) << store_cell_y[ia_all]
							<< std::setw(6) << store_cell_z[ia_all];
						}
					}

					ia_all++;
				}
			}
		}
	}

	// mohan fix bug 2012-04-02 by valgrind.
	delete[] store_cell_x;
	delete[] store_cell_y;
	delete[] store_cell_z;

	assert(ia_all == d_amount_expand);

	// becareful! now the cell is not the origin cell,
	// it's unitcell in expand case! so don't add
	// like x_min ... to the x_min_expand.
	this->x_min_expand = //this->x_min
	    - glayerX_minus * vec1[0]
	    - glayerY_minus * vec2[0]
	    - glayerZ_minus * vec3[0];

	this->y_min_expand = //this->y_min
	    - glayerX_minus * vec1[1]
	    - glayerY_minus * vec2[1]
	    - glayerZ_minus * vec3[1];

	this->z_min_expand = //this->z_min
	    - glayerX_minus * vec1[2]
	    - glayerY_minus * vec2[2]
	    - glayerZ_minus * vec3[2];

	// in fact we don't use these 6 parameters at all.
	this->x_max_expand = //this->x_max 
				   + glayerX * vec1[0]
	               + glayerY * vec2[0]
	               + glayerZ * vec3[0];

	this->y_max_expand = //this->y_max 
				   + glayerX * vec1[1]
	               + glayerY * vec2[1]
	               + glayerZ * vec3[1];

	this->z_max_expand = //this->z_max 
	  			   + glayerX * vec1[2]
	               + glayerY * vec2[2]
	               + glayerZ * vec3[2];

	if(test_atom_input)
	{
		GlobalV::ofs_running << " New Xmin=" << x_min_expand
			<< " Ymin=" << y_min_expand
			<< " Zmin=" << z_min_expand << std::endl;

		GlobalV::ofs_running << " New Xmax=" << x_max_expand
			<< " Ymax=" << y_max_expand
			<< " Zmax=" << z_max_expand << std::endl;
	}

	delete[] x_old;
	delete[] y_old;
	delete[] z_old;
	delete[] type_old;
	delete[] natom_old;
	return;
}

void Atom_input::calculate_cells(void)
{
	ModuleBase::TITLE("Atom_input", "calculate_cells");
//----------------------------------------------------------
// EXPLAIN :
// Expand_Case : Simple , we already know the cell numbres,
// all the trouble is only to construct adjacentset using all
// the cells.
// Not_Expand_Case : Using searching radius to construct
// the cells ,  trouble here,but is the convenience of searching
// time , we then only need to search 27-adjacent cell for each cell.
//----------------------------------------------------------
	if (expand_flag)
	{
		cell_nx = glayerX + glayerX_minus ;
		cell_ny = glayerY + glayerY_minus ;
		cell_nz = glayerZ + glayerZ_minus ;
	}
	else
	{
		// maybe a bug, if we don't use direct
		// coordinates, mohan note 2011-04-14
		double real_nx, real_ny, real_nz;
		real_nx = (x_max - x_min) / radius;
		real_ny = (y_max - y_min) / radius;
		real_nz = (z_max - z_min) / radius;
		cell_nx = static_cast<int>(real_nx) + 1;
		cell_ny = static_cast<int>(real_ny) + 1;
		cell_nz = static_cast<int>(real_nz) + 1;
	}

	//================
	// Wrong !
	//================
//	if(int_nx != real_nx) this->cell_nx++;
//	if(int_ny != real_ny) this->cell_ny++;
//	if(int_nz != real_nz) this->cell_nz++;
	//=======================================
	// Not need because if int_nx = real_nx,
	// the position belong to the next cell
	//=======================================
	return;
}

void Atom_input::set_FAtom(const UnitCell &ucell, FAtom &a)const
{
//----------------------------------------------------------
// EXPLAIN : if expand grid , set from array
//----------------------------------------------------------
	if (expand_flag)
	{
		a.setX(store_x[d_current]);
		a.setY(store_y[d_current]);
		a.setZ(store_z[d_current]);
		a.setType(store_type[d_current]);
		a.setNatom(store_natom[d_current]);
		++ d_current;
	}

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : load_atom
//
// EXPLAIN : if not expand grid , set directly
//----------------------------------------------------------
	else
	{
		Load_atom(ucell

		);
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		a.setType(type);
		a.setNatom(natom);
//		GlobalV::ofs_running<<"\n x = "<<x;
//		GlobalV::ofs_running<<"\n y = "<<y;
//		GlobalV::ofs_running<<"\n z = "<<z;
//		GlobalV::ofs_running<<"\n Type = "<<type;
//		GlobalV::ofs_running<<"\n natom = "<<natom;
	}

	return;
}

void Atom_input::Load_atom(const UnitCell& ucell)const
{
//	ModuleBase::TITLE("Atom_input","load_atom");
	natom++;

	if (natom >= ucell.atoms[type].na)
	{
		type ++;
		natom = 0;
	}

	x = ucell.atoms[type].tau[natom].x;

	y = ucell.atoms[type].tau[natom].y;
	z = ucell.atoms[type].tau[natom].z;

//	std::cout<<" x = "<<ucell.atoms[type].tau[natom].x
//		<<" y = "<<ucell.atoms[type].tau[natom].y
//		<<" z = "<<ucell.atoms[type].tau[natom].z
//		<<" type = "<<type
//		<<" natom = "<<natom;
	return;
}

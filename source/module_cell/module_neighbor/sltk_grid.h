#ifndef GRID_H
#define GRID_H

#include <stdexcept>
#include <functional>
#include "sltk_util.h"
#include "sltk_atom.h"
#include "sltk_atom_input.h"

#include "module_cell/unitcell.h"
#include <unordered_map>
#include <tuple>

struct KeyHash {
    std::size_t operator()(const std::tuple<double, double, double>& key) const {
        auto h1 = std::hash<double>{}(std::get<0>(key));
        auto h2 = std::hash<double>{}(std::get<1>(key));
        auto h3 = std::hash<double>{}(std::get<2>(key));
        return h1 ^ (h2 << 1) ^ (h3 << 2); // a hash combine
    }
};

struct KeyEqual {
    bool operator()(const std::tuple<double, double, double>& lhs,
                    const std::tuple<double, double, double>& rhs) const {
        return std::get<0>(lhs) == std::get<0>(rhs) &&
               std::get<1>(lhs) == std::get<1>(rhs) &&
               std::get<2>(lhs) == std::get<2>(rhs);
    }
};

typedef std::unordered_map<std::tuple<double, double, double>, FAtom, KeyHash, KeyEqual> AtomMap;

struct CellSet
{
	AtomMap atom_map;
	int in_grid[3];
	CellSet();
};

//==========================================================
// CLASS NAME :
// Atom_input : defined elsewhere
//==========================================================

class Atom_input;

//==========================================================
// CLASS NAME :
// Grid :
//==========================================================

class Grid
{
public:

	// Constructors and destructor
	// Grid is Global class,so init it with constant number
	Grid():test_grid(0){};
	Grid(const int &test_grid_in);
	virtual ~Grid();

	void init(std::ofstream &ofs,
				const UnitCell &ucell, 
				const Atom_input &input);

	//2015-05-07
	void delete_vector(const Atom_input &input);


	//Static data
	static const char* const ERROR[3];

	//Data
	int natom;// Total atoms.
	bool pbc; // periodic boundary condition
	bool expand_flag;
	double sradius;// searching radius
	double d_minX;// origin of all cells
	double d_minY;
	double d_minZ;
	int dx;
	int dy;
	int dz;
	int layer;
	double cell_x_length;
	double cell_y_length;
	double cell_z_length;
	int true_cell_x;
	int true_cell_y;
	int true_cell_z;

	std::vector<std::vector<std::vector<CellSet>>> Cell; //dx , dy ,dz is cell number in each direction,respectly.
	void delete_Cell() //it will replace by container soon!
	{
		if (this->init_cell_flag)
		{
			for (int i = 0;i < this->dx;i++)
			{
				for (int j = 0;j < this->dy;j++)
				{
					this->Cell[i][j].clear();
				}
			}

			for (int i = 0;i < this->dx;i++)
			{
				this->Cell[i].clear();
			}

			this->Cell.clear();
			this->init_cell_flag = false;
		}
	}

	double grid_length[3];
	double vec1[3];
	double vec2[3];
	double vec3[3];
	double lat_now;
	bool init_cell_flag;
    //LiuXh add 2019-07-15
    const double& getD_minX(void) const {return d_minX;}
    const double& getD_minY(void) const {return d_minY;}
    const double& getD_minZ(void) const {return d_minZ;}

    const int& getCellX(void) const {return dx;}
    const int& getCellY(void) const {return dy;}
    const int& getCellZ(void) const {return dz;}

private:

	const int test_grid;
//==========================================================
// MEMBER FUNCTIONS :
// Three Main Steps:
// NAME : setMemberVariables (read in datas from Atom_input,
// 			init cells.)
// NAME : setBoundaryAdjacent( Consider different situations,
// 			if not_expand case : nature/periodic boundary
// 			condition , if expand_case)
//==========================================================
	void setMemberVariables(std::ofstream &ofs_in, 
							const Atom_input &input);

	void setBoundaryAdjacent(std::ofstream &ofs_in, 
								const Atom_input &input);

//==========================================================

	void In_Which_Cell(const UnitCell &ucell, int &a, int &b, int &c, const FAtom &atom)const;
	void Build_Hash_Table(const UnitCell &ucell, const Atom_input &input);

//==========================================================

	void Construct_Adjacent_expand(const int i, const int j, const int k);

	void Construct_Adjacent_expand_periodic(const int i,
	                                        const int j,
											const int k,
											FAtom& fatom);

	void Construct_Adjacent_begin();
	void Construct_Adjacent_nature(
	    const int i, const int j, const int k, FAtom & fatom1);
	void Construct_Adjacent_periodic(
	    const int i, const int j, const int k, FAtom & fatom1);
	void Construct_Adjacent_final(const int i, const int j, const int k, FAtom & fatom1,
 									const int i2, const int j2, const int k2, FAtom & fatom2);
};

#endif

#ifndef GRID_H
#define GRID_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"
#include "sltk_atom_input.h"
#include "sltk_util.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

typedef std::vector<std::vector<FAtom>> AtomMap;

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
    Grid() : test_grid(0){};
    Grid(const int& test_grid_in);
    virtual ~Grid();

    void init(std::ofstream& ofs, const UnitCell& ucell, const Atom_input& input);

    // 2015-05-07
    void delete_vector(int i, int j, int k);

    // Data
    bool pbc; // periodic boundary condition
    bool expand_flag;
    double sradius2; // searching radius squared
    double sradius;  // searching radius
    double d_minX;   // origin of all cells
    double d_minY;
    double d_minZ;
    int cell_nx;
    int cell_ny;
    int cell_nz;
    int layer;

    int true_cell_x;
    int true_cell_y;
    int true_cell_z;

    std::vector<std::vector<std::vector<CellSet>>> Cell; // dx , dy ,dz is cell number in each direction,respectly.
    void delete_Cell()                                   // it will replace by container soon!
    {
        if (this->init_cell_flag)
        {
            for (int i = 0; i < this->cell_nx; i++)
            {
                for (int j = 0; j < this->cell_ny; j++)
                {
                    this->Cell[i][j].clear();
                }
            }

            for (int i = 0; i < this->cell_nx; i++)
            {
                this->Cell[i].clear();
            }

            this->Cell.clear();
            this->init_cell_flag = false;
        }
    }
    bool init_cell_flag = false;
    // LiuXh add 2019-07-15
    double getD_minX() const
    {
        return d_minX;
    }
    double getD_minY() const
    {
        return d_minY;
    }
    double getD_minZ() const
    {
        return d_minZ;
    }

    int getCellX() const
    {
        return cell_nx;
    }
    int getCellY() const
    {
        return cell_ny;
    }
    int getCellZ() const
    {
        return cell_nz;
    }
    int getTrueCellX() const
    {
        return true_cell_x;
    }
    int getTrueCellY() const
    {
        return true_cell_y;
    }
    int getTrueCellZ() const
    {
        return true_cell_z;
    }

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
    void setMemberVariables(std::ofstream& ofs_in, const Atom_input& input);

    void setBoundaryAdjacent(std::ofstream& ofs_in, const Atom_input& input);

    //==========================================================
    void Build_Hash_Table(const UnitCell& ucell, const Atom_input& input);

    //==========================================================

    void Construct_Adjacent_expand(const int i, const int j, const int k);

    void Construct_Adjacent_expand_periodic(const int i, const int j, const int k, FAtom& fatom);

    void Construct_Adjacent_begin();
    void Construct_Adjacent_nature(const int i, const int j, const int k, FAtom& fatom1);
    void Construct_Adjacent_periodic(const int i, const int j, const int k, FAtom& fatom1);
    void Construct_Adjacent_final(const int i,
                                  const int j,
                                  const int k,
                                  FAtom& fatom1,
                                  const int i2,
                                  const int j2,
                                  const int k2,
                                  FAtom& fatom2);
};

#endif

#ifndef ESOLVER_H
#define ESOLVER_H

#include "../input.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../src_pw/energy.h"
#include "../module_base/matrix.h"
//--------------temporary----------------------------
#include "src_lcao/record_adj.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_hamilt.h"
//--------------\temporary----------------------------
//------It should be moved as fast as possible------

namespace ModuleESolver
{

class ESolver
{
// protected:
//     ModuleBase::matrix lattice_v;
public:
    ESolver(){
        classname = "ESolver";
    }
    virtual ~ESolver(){};

    //virtual void Init(Input_EnSolver &inp, matrix &lattice_v)=0
    virtual void Init(Input &inp, UnitCell_pseudo &cell)=0;

    // They shoud be add after atom class is refactored
    // virtual void UpdateLatAtom(ModuleBase::matrix &lat_in, Atom &atom_in);
    // virtual void UpdateLat(ModuleBase::matrix &lat_in);
    // virtual void UpdateAtom(Atom &atom_in);
   
    virtual void Run(int istep, UnitCell_pseudo& cell) = 0;
    
    //Deal with exx and other calculation than scf/md/relax: 
    // such as nscf, istate-charge or envelope
    virtual void othercalculation(const int istep) {};

    virtual void cal_Energy(energy& en) = 0;
    virtual void cal_Force(ModuleBase::matrix &force) = 0;
    virtual void cal_Stress(ModuleBase::matrix& stress) = 0;
    virtual void postprocess() {};

    //Print current classname.
    void printname();

    //temporarily
    //get iterstep used in current scf
    virtual int getniter(){return 0;}
    string classname;
};

void init_esolver(ESolver* &p_esolver, const string use_esol);
void clean_esolver(ESolver* &pesolver);

}

#endif

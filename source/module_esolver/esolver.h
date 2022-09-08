#ifndef ESOLVER_H
#define ESOLVER_H

#include "../input.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../src_pw/energy.h"
#include "../module_base/matrix.h"
//--------------temporary----------------------------
#include "module_psi/psi.h"

namespace ModuleESolver
{
    class ESolver
    {
        // protected:
        //     ModuleBase::matrix lattice_v;
    public:
        ESolver() {
            classname = "ESolver";
        }
        
        virtual ~ESolver() 
        {
            //--------------temporary----------------------------
            if(this->psi != nullptr)
            {
                delete psi;
            }
            if(this->psid != nullptr)
            {
                delete psid;
            }
        }

        //virtual void Init(Input_EnSolver &inp, matrix &lattice_v)=0
        virtual void Init(Input& inp, UnitCell_pseudo& cell) = 0;

        // They shoud be add after atom class is refactored
        // virtual void UpdateLatAtom(ModuleBase::matrix &lat_in, Atom &atom_in);
        // virtual void UpdateLat(ModuleBase::matrix &lat_in);
        // virtual void UpdateAtom(Atom &atom_in);

        virtual void Run(int istep, UnitCell_pseudo& cell) = 0;

        //Deal with exx and other calculation than scf/md/relax: 
        // such as nscf, istate-charge or envelope
        virtual void othercalculation(const int istep) {};

        virtual void cal_Energy(energy& en) = 0;
        virtual void cal_Force(ModuleBase::matrix& force) = 0;
        virtual void cal_Stress(ModuleBase::matrix& stress) = 0;
        virtual void postprocess() {};

        //Print current classname.
        void printname();

        //temporarily
        //get iterstep used in current scf
        virtual int getniter() { return 0; }
        string classname;

        //--------------temporary----------------------------
        // this is the interface of non-self-consistant calculation
        virtual void nscf() {};
        
        //wavefunction coefficients
        psi::Psi<std::complex<double>>* psi = nullptr;
        psi::Psi<double>* psid = nullptr;
    };

    string determine_type();
    void init_esolver(ESolver*& p_esolver);
    void clean_esolver(ESolver*& pesolver);

}

#endif

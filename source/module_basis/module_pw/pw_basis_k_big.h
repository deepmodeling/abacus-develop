#ifndef PW_BASIS_K_BIG_H
#define PW_BASIS_K_BIG_H
#include "module_base/constants.h"
#include "module_base/global_function.h"

// temporary class, because previous ABACUS consider big grid for fft grids 
// which are used for grid integration in LCAO.
// In fact, it is unnecessary. It will be moved after grid integration is refactored.
namespace ModulePW
{

class PW_Basis_K_Big: public PW_Basis_K
{
public:
    
    // combine [bx,by,bz] FFT grids into a big one
	// typical values are bx=2, by=2, bz=2
	// nbx=nx/bx, nby=ny/by, nbz=nz/bz, 
    // Note: this class can only use initgrids(lat0_in, latvec_in, PW_Basis_Big::nx, PW_Basis_Big::ny, PW_Basis_Big::nz)!!!
    PW_Basis_K_Big(){
        bx = 1;
        by = 1;
        bz = 1;
    }
    PW_Basis_K_Big(std::string device_, std::string precision_) : PW_Basis_K(device_, precision_) {}
    ~PW_Basis_K_Big(){};
    void setbxyz(const int bx_in, const int by_in, const int bz_in)
    {
        bx = bx_in;
        by = by_in;
        bz = bz_in;
    }
    int bx,by,bz;
    int nbx, nby, nbz;
    virtual void distribute_r()
    {
        bx = (bx == 0) ? 2 : bx;
        by = (by == 0) ? 2 : by;
        bz = (bz == 0) ? 2 : bz;
        this->nbx = this->nx / bx;
        this->nby = this->ny / by;
        this->nbz = this->nz / bz;
        delete[] this->numz; this->numz = new int[this->poolnproc];
        delete[] this->startz; this->startz = new int[this->poolnproc];
        ModuleBase::GlobalFunc::ZEROS(this->numz, this->poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->startz, this->poolnproc);

        int npbz = this->nbz / this->poolnproc;
        int modbz = this->nbz % this->poolnproc;
        this->startz[0] = 0;
        for(int ip = 0 ; ip < this->poolnproc ; ++ip)
        {
            this->numz[ip] = npbz*this->bz;
            if(ip < modbz)   this->numz[ip]+=this->bz;
            if(ip < this->poolnproc - 1)   this->startz[ip+1] = this->startz[ip] + numz[ip];
            if(ip == this->poolrank) 
            {
                this->nplane = numz[ip];
                this->startz_current = startz[ip];
            }
        }
        this->nrxx = this->numz[this->poolrank] * this->nxy;
        return;
    }

};
}
#endif
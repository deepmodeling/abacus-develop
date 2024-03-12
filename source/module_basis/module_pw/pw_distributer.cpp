#include "pw_basis.h"
#include "module_base/global_function.h"

namespace ModulePW
{
/// 
/// distribute real-space grids to different processors
/// Known: nx, ny, nz, poolnproc, poolrank
/// output: nrxx, startz, numz
/// 
void PW_Basis::distribute_r()
{
    delete[] this->numz; this->numz = new int[this->poolnproc];
    delete[] this->startz; this->startz = new int[this->poolnproc];
    ModuleBase::GlobalFunc::ZEROS(this->numz, this->poolnproc);
    ModuleBase::GlobalFunc::ZEROS(this->startz, this->poolnproc);

    int npz = this->nz / this->poolnproc;
    int modz = this->nz % this->poolnproc;
    this->startz[0] = 0;
    for(int ip = 0 ; ip < this->poolnproc ; ++ip)
    {
        this->numz[ip] = npz;
        if(ip < modz)   this->numz[ip]++;
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

}
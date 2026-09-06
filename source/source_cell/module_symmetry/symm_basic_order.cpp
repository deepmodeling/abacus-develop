#include "symmetry.h"
#include "source_base/mymath.h"

namespace ModuleSymmetry
{
// atom ordering for each atom type
// by a "weighted function" f
// (instead of ordering by x, y, z directly)
void Symmetry_Basic::atom_ordering_new(double *posi, const int natom, int *subindex) const
{
    //order the atomic positions inside a supercell by a unique ordering scheme
    subindex[0] = 0;

    if(natom == 1)
    {
        //if there is only one atom, it is not necessary to order
        return;
    }

    std::vector<double> tmpx(natom);
    std::vector<double> tmpy(natom);
    std::vector<double> tmpz(natom);
    for(int i=0; i<natom; i++)
    {
        tmpx[i] = posi[i*3];
        tmpy[i] = posi[i*3+1];
        tmpz[i] = posi[i*3+2];
    }
    double x_max = *max_element(tmpx.begin(),tmpx.end());
    double x_min = *min_element(tmpx.begin(),tmpx.end());
    double y_max = *max_element(tmpy.begin(),tmpy.end());
    double y_min = *min_element(tmpy.begin(),tmpy.end());
    double z_max = *max_element(tmpz.begin(),tmpz.end());
    double z_min = *min_element(tmpz.begin(),tmpz.end());

    double* weighted_func = new double[natom];

    //the first time: f(x, y, z)
    for(int i=0; i<natom; i++)
    {
        weighted_func[i]=1/epsilon/epsilon*tmpx[i]+1/epsilon*tmpy[i]+tmpz[i];
    }
    ModuleBase::heapsort(natom, weighted_func, subindex);
    this->order_atoms(posi, natom, subindex);
    for(int i=0; i<natom; i++)
    {
        tmpx[i] = posi[i*3];
        tmpy[i] = posi[i*3+1];
        tmpz[i] = posi[i*3+2];
    }

    //the second time: f(y, z) for fixed x
    for(int i=0; i<natom-1;)
    {
        int ix_right=i+1;    //right bound is no included
        while(ix_right<natom && equal(tmpx[ix_right],tmpx[i]))
        {
            ++ix_right;
        }

        int nxequal=ix_right-i;
        if(nxequal>1)    //need a new sort
        {
            subindex[0] = 0;
            for(int j=0; j<nxequal; ++j)
            {
                weighted_func[j]=1/epsilon*tmpy[i+j]+tmpz[i+j];
            }
            ModuleBase::heapsort(nxequal, weighted_func, subindex);
            this->order_atoms(&posi[i*3], nxequal, subindex);
        }
        i=ix_right;
    }

    delete[] weighted_func;
    return;
}

void Symmetry_Basic::test_atom_ordering(double *posi, const int natom, int *subindex) const
{
    //an interface to test a protected function
    this->atom_ordering_new(posi, natom, subindex);
}
}

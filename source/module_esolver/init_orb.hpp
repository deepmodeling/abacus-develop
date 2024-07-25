#ifndef INIT_ORB_HPP
#define INIT_ORB_HPP
#include "module_cell/unitcell.h"
#include "module_base/memory.h"
#include "module_basis/module_ao/ORB_read.h"
void inline init_orb(double& dr_uniform, 
                        std::vector<double>& rcuts,
                        UnitCell& ucell,
                        std::vector<std::vector<double>>& psi_u)
{
    // set the grid parameters
    dr_uniform=GlobalC::ORB.dr_uniform;
    
    const int nwmax=ucell.nwmax;
    const int ntype=ucell.ntype;
    
    rcuts=std::vector<double>(ntype);
    ModuleBase::Memory::record("rcuts", sizeof(double)*ntype*3);
    for(int T=0; T<ntype; T++)
    {
        rcuts[T]=GlobalC::ORB.Phi[T].getRcut();
    }
    
    const double max_cut = *std::max_element(rcuts.begin(), rcuts.end());
    const int nr_max = static_cast<int>(1/dr_uniform * max_cut) + 10;
    psi_u=std::vector<std::vector<double>>(ntype * nwmax);
    ModuleBase::Memory::record("psi_u", sizeof(double)*nwmax*ntype*2);
    
    Atom* atomx;
    const Numerical_Orbital_Lm* pointer;
    
    for (int i = 0; i < ntype; i++)
    {
        atomx = &ucell.atoms[i];
        for (int j = 0; j < nwmax; j++)
        {
            if (j < atomx->nw)
            {
                pointer = &GlobalC::ORB.Phi[i].PhiLN(atomx->iw2l[j],atomx->iw2n[j]);
                for (int k=0;k<pointer->psi_uniform.size();k++)
                {
                    psi_u[i*nwmax+j].push_back(pointer->psi_uniform[k]);
                    psi_u[i*nwmax+j].push_back(pointer->dpsi_uniform[k]);
                    // dpsi_u[i*nwmax+j].push_back(pointer->dpsi_uniform[k]);
                }
                // When the pointer points to nr_uniform, 
                // take up to 8 numbers backwards when calculating force, 
                // and set them to 0.
                for (int k=0;k<8;k++)
                {
                    psi_u[i*nwmax+j].push_back(0.0);
                }
            }
        }
    }
}
#endif
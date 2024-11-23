#include "gint_tools.h"
#include "module_base/memory.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/unitcell.h"

namespace Gint_Tools{

void init_orb(double& dr_uniform, 
                std::vector<double>& rcuts,
                UnitCell& ucell,
                const LCAO_Orbitals& orb,
                std::vector<std::vector<double>>& psi_u,
                std::vector<std::vector<double>>& dpsi_u,
                std::vector<std::vector<double>>& d2psi_u)
{
    //! set the grid parameters
    dr_uniform=orb.dr_uniform;

    assert(dr_uniform>0.0);
    
    const int nwmax=ucell.nwmax;
    const int ntype=ucell.ntype;

    assert(nwmax>0);
    assert(ntype>0);
    
    rcuts=std::vector<double>(ntype);
    ModuleBase::Memory::record("rcuts", sizeof(double)*ntype*3);
    
    for(int it=0; it<ntype; it++)
	{
		rcuts[it]=orb.Phi[it].getRcut();
	}
    
    const double max_cut = *std::max_element(rcuts.begin(), rcuts.end());
    const int nr_max = static_cast<int>(1/dr_uniform * max_cut) + 10;
    psi_u=std::vector<std::vector<double>>(ntype * nwmax);
    dpsi_u=std::vector<std::vector<double>>(ntype * nwmax);
    d2psi_u=std::vector<std::vector<double>>(ntype * nwmax);
    ModuleBase::Memory::record("psi_u", sizeof(double)*nwmax*ntype*3);
    
    Atom* atomx = nullptr;
    const Numerical_Orbital_Lm* pointer = nullptr;
    
    for (int i = 0; i < ntype; i++)
    {
        atomx = &ucell.atoms[i];
        for (int j = 0; j < nwmax; j++)
        {
            const int k=i*nwmax+j;
            if (j < atomx->nw)
            {
                pointer = &orb.Phi[i].PhiLN(atomx->iw2l[j],atomx->iw2n[j]);
                psi_u[k]=pointer->psi_uniform;
                dpsi_u[k]=pointer->dpsi_uniform;
                d2psi_u[k]=pointer->ddpsi_uniform;
            }
        }
    }
}// End of init_orb()

}// End of Gint_Tools

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
    // set the grid parameters
    dr_uniform=orb.dr_uniform;
    
    const int nwmax=ucell.nwmax;
    const int ntype=ucell.ntype;
    
    rcuts=std::vector<double>(ntype);
    ModuleBase::Memory::record("rcuts", sizeof(double)*ntype*3);
    for(int T=0; T<ntype; T++)
	{
		rcuts[T]=orb.Phi[T].getRcut();
	}
    
    const double max_cut = *std::max_element(rcuts.begin(), rcuts.end());
    const int nr_max = static_cast<int>(1/dr_uniform * max_cut) + 10;
    psi_u=std::vector<std::vector<double>>(ntype * nwmax);
    dpsi_u=std::vector<std::vector<double>>(ntype * nwmax);
    d2psi_u=std::vector<std::vector<double>>(ntype * nwmax);
    ModuleBase::Memory::record("psi_u", sizeof(double)*nwmax*ntype*3);
    
    Atom* atomx;
    const Numerical_Orbital_Lm* pointer;
    
    for (int i = 0; i < ntype; i++)
    {
        atomx = &ucell.atoms[i];
        for (int j = 0; j < nwmax; j++)
        {
            if (j < atomx->nw)
            {
                pointer = &orb.Phi[i].PhiLN(atomx->iw2l[j],atomx->iw2n[j]);
                psi_u[i*nwmax+j]=pointer->psi_uniform;
                dpsi_u[i*nwmax+j]=pointer->dpsi_uniform;
                d2psi_u[i*nwmax+j]=pointer->ddpsi_uniform;
            }
        }
    }
}
}// Gint_Tools

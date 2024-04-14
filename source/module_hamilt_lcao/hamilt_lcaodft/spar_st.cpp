#include <complex>
#include "spar_st.h"
#include "spar_dh.h"
#include "spar_hsr.h"

void sparse_format::cal_SR(
        const Parallel_Orbitals &pv,
		std::set<Abfs::Vector3_Order<int>> &all_R_coor,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &SR_sparse,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &SR_soc_sparse,
		Grid_Driver &grid,
		const double &sparse_thr, 
		hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("sparse_format","cal_SR");

    sparse_format::set_R_range(all_R_coor, grid);

    const int nspin = GlobalV::NSPIN;

    //cal_STN_R_sparse(current_spin, sparse_thr);
    if(nspin==1 || nspin==2)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao 
        = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);
        const int cspin = 0;
        sparse_format::cal_HContainer_d(pv, cspin, sparse_thr, *(p_ham_lcao->getSR()), SR_sparse);
    }
    else if(nspin==4)
    {
        hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao 
        = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);
        const int cspin = 0;
        sparse_format::cal_HContainer_cd(pv, cspin, sparse_thr, *(p_ham_lcao->getSR()), SR_soc_sparse);
    }

    return;
}


void sparse_format::cal_TR(
        const UnitCell &ucell,
        const Parallel_Orbitals &pv,
        LCAO_Matrix &lm,
	    Grid_Driver &grid,
		LCAO_gen_fixedH &gen_h,
		const double &sparse_thr)
{
    ModuleBase::TITLE("sparse_format","cal_TR");
    
    //need to rebuild T(R)
    lm.Hloc_fixedR.resize(lm.ParaV->nnr);
    lm.zeros_HSR('T');

    gen_h.build_ST_new('T', 0, ucell, GlobalC::ORB, GlobalC::UOT, &(GlobalC::GridD), lm.Hloc_fixedR.data());

    sparse_format::set_R_range(lm.all_R_coor, grid);

	sparse_format::cal_STN_R_for_T(
			ucell, 
			pv, 
			lm, 
			grid, 
			sparse_thr);

    return;
}


void sparse_format::cal_STN_R_for_T(
        const UnitCell &ucell,
        const Parallel_Orbitals &pv,
        LCAO_Matrix &lm,
		Grid_Driver &grid,
		const double &sparse_thr)
{
    ModuleBase::TITLE("sparse_format","cal_STN_R_for_T");

    const int nspin = GlobalV::NSPIN;

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double tmp=0.0;
    std::complex<double> tmpc=std::complex<double>(0.0,0.0);

    for(int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            grid.Find_atom(ucell, tau1, T1, I1);
            Atom* atom1 = &ucell.atoms[T1];
            const int start = ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < grid.getAdjacentNum()+1; ++ad)
            {
                const int T2 = grid.getType(ad);
                const int I2 = grid.getNatom(ad);
                Atom* atom2 = &ucell.atoms[T2];

                tau2 = grid.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) 
				{
					adj = true;
				}

                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < grid.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = grid.getType(ad0);

                        tau0 = grid.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = ucell.itiaiw2iwt(T2,I2,0);

					Abfs::Vector3_Order<int> dR(
							grid.getBox(ad).x, 
							grid.getBox(ad).y, 
							grid.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = pv.global2local_row(iw1_all);

						if(mu<0)
						{
							continue;
						}

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = pv.global2local_col(iw2_all);

							if(nu<0)
							{
								continue;
							}

                            if(nspin==1 || nspin==2)
                            {
                                tmp = lm.Hloc_fixedR[index];
                                if (std::abs(tmp) > sparse_thr)
                                {
                                    lm.TR_sparse[dR][iw1_all][iw2_all] = tmp;
                                }
                            }
                            else if(nspin==4)
                            {
                                tmpc = lm.Hloc_fixedR_soc[index];
                                if(std::abs(tmpc) > sparse_thr)
                                {
                                    lm.TR_soc_sparse[dR][iw1_all][iw2_all] = tmpc;
                                }
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    return;
}

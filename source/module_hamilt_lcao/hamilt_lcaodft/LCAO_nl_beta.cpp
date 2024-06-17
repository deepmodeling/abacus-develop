#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_base/timer.h"


#ifdef __MKL
#include <mkl_service.h>  // mkl_get_max_threads
#endif

namespace LCAO_domain
{

void build_Nonlocal_beta_new(
    LCAO_Matrix &lm,
    double* HSloc,
	const UnitCell &ucell,
	const LCAO_Orbitals& orb,
	const ORB_gen_tables& uot,
	Grid_Driver* GridD) //update by liuyu 2021-04-07
{
    ModuleBase::TITLE("LCAO_domain","b_NL_beta_new");
    ModuleBase::timer::tick ("LCAO_domain","b_NL_beta_new");

	const Parallel_Orbitals* pv = lm.ParaV;
	const int npol = GlobalV::NPOL;

#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif

    const std::vector<AdjacentAtomInfo> adjs_all = GridD->get_adjs(ucell);

#ifdef _OPENMP
    #pragma omp parallel
    {
        double* Nonlocal_thread;
        Nonlocal_thread = new double[pv->nloc];
        ModuleBase::GlobalFunc::ZEROS(Nonlocal_thread, pv->nloc);
        #pragma omp for schedule(dynamic)
#endif
        for(int iat=0; iat<ucell.nat; iat++)
        {
            const int T0 = ucell.iat2it[iat];
            const int I0 = ucell.iat2ia[iat];
            Atom* atom0 = &ucell.atoms[T0];

            //=======================================================
            //Step1:
            //saves <beta|psi>, where beta runs over L0,M0 on atom I0
            //and psi runs over atomic basis sets on the current core
            //=======================================================
            #ifdef _OPENMP
                std::vector<std::unordered_map<int,std::vector<double>>> nlm_tot_thread;
                nlm_tot_thread.resize(adjs_all[iat].adj_num + 1);
            #else 
                std::vector<std::unordered_map<int,std::vector<double>>> nlm_tot;
                nlm_tot.resize(adjs_all[iat].adj_num + 1);
            #endif 

            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            const double Rcut_Beta = ucell.infoNL.Beta[T0].get_rcut_max();

            //outermost loop : all adjacent atoms
            for(int ad_count=0; ad_count < adjs_all[iat].adj_num + 1; ad_count++)
            {
                const int T1 = adjs_all[iat].ntype[ad_count];
                const int I1 = adjs_all[iat].natom[ad_count];
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const double Rcut_AO1 = orb.Phi[T1].getRcut();
                const ModuleBase::Vector3<double> tau1 = adjs_all[iat].adjacent_tau[ad_count];
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw*npol;

                #ifdef _OPENMP
                    nlm_tot_thread[ad_count].clear();
                #else 
                    nlm_tot[ad_count].clear();
                #endif 

                //middle loop : atomic basis on current processor (either row or column)
                const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                if (dist1 > Rcut_Beta + Rcut_AO1)
                {
                    continue;
                }

                for(int iw1=0; iw1<nw1_tot; ++iw1)
                {
                    const int iw1_all = start1 + iw1;
                    const int iw1_local = pv->global2local_row(iw1_all);
                    const int iw2_local = pv->global2local_col(iw1_all);

                    if(iw1_local < 0 && iw2_local < 0) continue;

                    const int iw1_0 = iw1/npol;
                    std::vector<std::vector<double>> nlm;
                    //2D, but first dimension is only 1 here
                    //for force, the right hand side is the gradient
                    //and the first dimension is then 3
                    //inner loop : all projectors (L0,M0)


#ifdef USE_NEW_TWO_CENTER
                    //=================================================================
                    //          new two-center integral (temporary)
                    //=================================================================
                    int L1 = atom1->iw2l[ iw1_0 ];
                    int N1 = atom1->iw2n[ iw1_0 ];
                    int m1 = atom1->iw2m[ iw1_0 ];

                    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                    int M1 = (m1 % 2 == 0) ? -m1/2 : (m1+1)/2;

                    ModuleBase::Vector3<double> dtau = ucell.atoms[T0].tau[I0] - tau1;
                    uot.two_center_bundle->overlap_orb_beta->snap(
                            T1, L1, N1, M1, T0, dtau * ucell.lat0, false, nlm);
#else
                    uot.snap_psibeta_half(
                        orb,
                        ucell.infoNL,
                        nlm, tau1, T1,
                        atom1->iw2l[ iw1_0 ], // L1
                        atom1->iw2m[ iw1_0 ], // m1
                        atom1->iw2n[ iw1_0 ], // N1
                        ucell.atoms[T0].tau[I0], T0, 0); //R0,T0
#endif
                    //=================================================================
                    //          end of new two-center integral (temporary)
                    //=================================================================

                    #ifdef _OPENMP
                        nlm_tot_thread[ad_count].insert({iw1_all,nlm[0]});
                    #else 
                        nlm_tot[ad_count].insert({iw1_all,nlm[0]});
                    #endif 
                }//end iw
            }//end ad

            //=======================================================
            //Step2:
            //calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
            //and accumulate the value to Hloc_fixed(i,j)
            //=======================================================
            for(int ad1_count=0; ad1_count < adjs_all[iat].adj_num + 1; ad1_count++)
            {
                const int T1 = adjs_all[iat].ntype[ad1_count];
                const int I1 = adjs_all[iat].natom[ad1_count];
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = adjs_all[iat].adjacent_tau[ad1_count];
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw*npol;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                for (int ad2_count=0; ad2_count < adjs_all[iat].adj_num + 1; ad2_count++)
                {
                    const int T2 = adjs_all[iat].ntype[ad2_count];
                    const int I2 = adjs_all[iat].natom[ad2_count];
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = adjs_all[iat].adjacent_tau[ad2_count];
                    const Atom* atom2 = &ucell.atoms[T2];
                    const int nw2_tot = atom2->nw*npol;
                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                    const double dist2 = (tau2-tau0).norm() * ucell.lat0;

                    if (dist1 > Rcut_Beta + Rcut_AO1
                            || dist2 > Rcut_Beta + Rcut_AO2)
                    {
                        continue;
                    }

                    for(int iw1=0; iw1<nw1_tot; ++iw1)
                    {
                        const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_row(iw1_all);
                        if(iw1_local < 0) continue;
                        const int iw1_0 = iw1/npol;
                        for(int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_col(iw2_all);
                            if(iw2_local < 0) continue;
                            const int iw2_0 = iw2/npol;
                            #ifdef _OPENMP
                                std::vector<double> nlm1 = nlm_tot_thread[ad1_count][iw1_all];
                                std::vector<double> nlm2 = nlm_tot_thread[ad2_count][iw2_all];
                            #else 
                                std::vector<double> nlm1 = nlm_tot[ad1_count][iw1_all];
                                std::vector<double> nlm2 = nlm_tot[ad2_count][iw2_all];
                            #endif 

                            assert(nlm1.size()==nlm2.size());
                            #ifdef _OPENMP
                                double nlm_thread=0.0;
                            #else 
                                double nlm=0.0;
                            #endif
                            const int nproj = ucell.infoNL.nproj[T0];
                            int ib = 0;
                            for(int nb = 0; nb < nproj; nb++)
                            {
                                const int L0 = ucell.infoNL.Beta[T0].Proj[nb].getL();
                                for(int m=0;m<2*L0+1;m++)
                                {
                                    #ifdef _OPENMP
                                        nlm_thread += nlm1[ib]*nlm2[ib]*ucell.atoms[T0].ncpp.dion(nb,nb);
                                    #else 
                                        nlm += nlm1[ib]*nlm2[ib]*ucell.atoms[T0].ncpp.dion(nb,nb);
                                    #endif
                                    ib+=1;
                                }
                            }
                            assert(ib==nlm1.size());

                            const int ir = pv->global2local_row(iw1_all);
                            const int ic = pv->global2local_col(iw2_all);
                            long index=0;
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                index=ic*pv->nrow+ir;
                            }
                            else
                            {
                                index=ir*pv->ncol+ic;
                            }
#ifdef _OPENMP
                            Nonlocal_thread[index] += nlm_thread;
#else
                            lm.set_HSgamma(iw1_all, iw2_all, nlm, HSloc);
#endif
                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end iat

        #ifdef _OPENMP
            #pragma omp critical(cal_nonlocal)
            {
                for(int i=0; i<pv->nloc; i++)
                {
                    lm.Hloc_fixed[i] += Nonlocal_thread[i];
                }
            }
            delete[] Nonlocal_thread;
        #endif
#ifdef _OPENMP
    }
#endif

#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
	
    ModuleBase::timer::tick ("LCAO_domain","b_NL_beta_new");
	return;
}

}

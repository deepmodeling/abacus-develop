//wenfei 2022-1-11
//This file contains 2 subroutines:
//1. build_psialpha, which calculates the overlap
//between atomic basis and projector alpha : <psi_mu|alpha>
//which will be used in calculating pdm, gdmx, H_V_delta, F_delta;
//2. check_psialpha, which prints the results into .dat files
//for checking

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/vector3.h"
#include "module_base/timer.h"

void LCAO_Deepks::build_psialpha(const bool& calc_deri,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const ORB_gen_tables &UOT)
{
    ModuleBase::TITLE("LCAO_Deepks", "build_psialpha");
    ModuleBase::timer::tick ("LCAO_Deepks","build_psialpha");

    //cutoff for alpha is same for all types of atoms
    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    
    int job;
    if(!calc_deri)
    {
        job=0;
    }
    else
    {
        job=1;
    }

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            //iat: atom index on which |alpha> is located
            const int iat = ucell.itia2iat(T0,I0);
			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                this->nlm_save[iat].resize(GridD.getAdjacentNum()+1);
            }
            //outermost loop : find all adjacent atoms
            for (int ad=0; ad<GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;

                std::unordered_map<int,std::vector<std::vector<double>>> nlm_cur;
                if(GlobalV::GAMMA_ONLY_LOCAL)
                {
                    this->nlm_save[iat][ad].clear();
                }
                else
                {
                    nlm_cur.clear();
                }

				const double dist1 = (tau1-tau0).norm() * ucell.lat0;

				if (dist1 > Rcut_Alpha + Rcut_AO1)
				{
					continue;
				}

                //middle loop : all atomic basis on the adjacent atom ad
				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
                    const int iw1_local = pv->global2local_row(iw1_all);
                    const int iw2_local = pv->global2local_col(iw1_all);
					if(iw1_local < 0 && iw2_local < 0)continue;
					const int iw1_0 = iw1/GlobalV::NPOL;
					std::vector<std::vector<double>> nlm;
					//2D, dim 0 contains the overlap <psi|alpha>
                    //dim 1-3 contains the gradient of overlap

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
                    GlobalC::UOT.two_center_bundle->overlap_orb_alpha->snap(
                            T1, L1, N1, M1, 0, dtau * ucell.lat0, calc_deri, nlm);
#else
					//inner loop : all projectors (N,L,M)
					UOT.snap_psialpha_half(
                        orb,
						nlm, job, tau1, T1,
						atom1->iw2l[ iw1_0 ], // L1
						atom1->iw2m[ iw1_0 ], // m1
						atom1->iw2n[ iw1_0 ], // N1
						ucell.atoms[T0].tau[I0], T0, I0); //R0,T0
#endif
                    //=================================================================
                    //          end of new two-center integral (temporary)
                    //=================================================================

                    if(GlobalV::GAMMA_ONLY_LOCAL)
                    {
                        this->nlm_save[iat][ad].insert({iw1_all,nlm});
                    }
                    else
                    {
                        nlm_cur.insert({iw1_all,nlm});
                    }
				}//end iw

                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    const int ibt=ucell.itia2iat(T1, I1);
                    const int rx=GridD.getBox(ad).x;
                    const int ry=GridD.getBox(ad).y;
                    const int rz=GridD.getBox(ad).z;
                    key_tuple key_1(ibt,rx,ry,rz);
                    this->nlm_save_k[iat][key_1]=nlm_cur;
                }
			}//end ad
		}//end I0
	}//end T0

    ModuleBase::timer::tick ("LCAO_Deepks","build_psialpha");
	return;

}

void LCAO_Deepks::check_psialpha(const bool& calc_deri,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const ORB_gen_tables &UOT)
{
    ModuleBase::TITLE("LCAO_Deepks", "check_psialpha");
    ModuleBase::timer::tick ("LCAO_Deepks","check_psialpha");

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    //same for all types of atoms
    int job;
    if(!calc_deri)
    {
        job=0;
    }
    else
    {
        job=1;
    }

    std::ofstream ofs("psialpha.dat");
    std::ofstream ofs_x("dpsialpha_x.dat");
    std::ofstream ofs_y("dpsialpha_y.dat");
    std::ofstream ofs_z("dpsialpha_z.dat");

    ofs<<std::setprecision(10);
    ofs_x<<std::setprecision(10);
    ofs_y<<std::setprecision(10);
    ofs_z<<std::setprecision(10);

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
			//=======================================================
            //Step 1 : 
			//saves <alpha|psi>, where alpha runs over all projectors
			//and psi runs over atomic basis sets on the current core
			//=======================================================

			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            ofs << "iat : " << iat << std::endl;
            ofs_x << "iat : " << iat << std::endl;
            ofs_y << "iat : " << iat << std::endl;
            ofs_z << "iat : " << iat << std::endl;

            for (int ad=0; ad<GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;

				const double dist1 = (tau1-tau0).norm() * ucell.lat0;

				if (dist1 > Rcut_Alpha + Rcut_AO1)
				{
					continue;
				}

                int ibt, rx, ry, rz;
                if(GlobalV::GAMMA_ONLY_LOCAL)
                {
                    ofs << "ad : " << ad << " " << dist1 << std::endl;
                    ofs_x << "ad : " << ad << " " << dist1 << std::endl;
                    ofs_y << "ad : " << ad << " " << dist1 << std::endl;
                    ofs_z << "ad : " << ad << " " << dist1 << std::endl;
                }
                else
                {
                    ibt=ucell.itia2iat(T1, I1);
                    rx=GridD.getBox(ad).x;
                    ry=GridD.getBox(ad).y;
                    rz=GridD.getBox(ad).z;
                    ofs << "key : " << ibt << " " << rx << " " << ry << " " << rz << std::endl;
                    ofs_x << "key : " << ibt << " " << rx << " " << ry << " " << rz << std::endl; 
                    ofs_y << "key : " << ibt << " " << rx << " " << ry << " " << rz << std::endl; 
                    ofs_z << "key : " << ibt << " " << rx << " " << ry << " " << rz << std::endl; 
                }

				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
                    ofs << "iw : " << iw1_all << std::endl;
                    ofs_x << "iw : " << iw1_all << std::endl;
                    ofs_y << "iw : " << iw1_all << std::endl;
                    ofs_z << "iw : " << iw1_all << std::endl;
                    const int iw1_local = pv->global2local_row(iw1_all);
                    const int iw2_local = pv->global2local_col(iw1_all);
					if(iw1_local < 0 && iw2_local < 0)continue;
					const int iw1_0 = iw1/GlobalV::NPOL;

                    std::vector<std::vector<double>> nlm;

                    if(GlobalV::GAMMA_ONLY_LOCAL)
                    {
                        nlm = this->nlm_save[iat][ad][iw1_all];
                    }
                    else
                    {
                        key_tuple key_1(ibt,rx,ry,rz);
                        nlm = this->nlm_save_k[iat][key_1][iw1_all];
                    }
                    
                    for(int ind=0;ind<nlm[0].size();ind++)
                    {
                        ofs << nlm[0][ind] << " ";
                        if(ind%6 == 5) ofs << "\n";
                        if(calc_deri)
                        {
                            ofs_x << nlm[1][ind] << " ";
                            if(ind%6 == 5) ofs_x << "\n";
                            ofs_y << nlm[2][ind] << " ";
                            if(ind%6 == 5) ofs_y << "\n";
                            ofs_z << nlm[3][ind] << " ";
                            if(ind%6 == 5) ofs_z << "\n";
                        }
                    }
                    ofs << std::endl;
                    if(calc_deri)
                    {
                        ofs_x << std::endl;
                        ofs_y << std::endl;
                        ofs_z << std::endl;
                    }
				}//end iw
			}//end ad
		}//end I0
	}//end T0

    ModuleBase::timer::tick ("LCAO_Deepks","check_psialpha");
	return;

}
#endif

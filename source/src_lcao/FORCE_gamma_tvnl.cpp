#include "FORCE_gamma.h"
#include "../src_pw/global.h"


void Force_LCAO_gamma::cal_ftvnl_dphi(
	ModuleBase::matrix& dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& ftvnl_dphi, 
	ModuleBase::matrix& stvnl_dphi)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_ftvnl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi");
    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        const int iat = GlobalC::ucell.iwt2iat[i];
        for(int j=0; j<GlobalV::NLOCAL; j++)
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];

            if (mu >= 0 && nu >= 0 )
            {
                const int index = mu * GlobalC::ParaO.ncol + nu;
                //contribution from deriv of AO's in T+VNL term
                
                double sum = 0.0;
                for(int is=0; is<GlobalV::NSPIN; ++is)
                {
                    sum += dm2d(is, index);
                }
                sum *= 2.0;

				if(isforce)
				{
					ftvnl_dphi(iat,0) += sum * GlobalC::LM.DHloc_fixed_x[index];
					ftvnl_dphi(iat,1) += sum * GlobalC::LM.DHloc_fixed_y[index];
					ftvnl_dphi(iat,2) += sum * GlobalC::LM.DHloc_fixed_z[index];
				}
                if(isstress)
                {
                    stvnl_dphi(0,0) += sum/2.0 * GlobalC::LM.DHloc_fixed_11[index];
                    stvnl_dphi(0,1) += sum/2.0 * GlobalC::LM.DHloc_fixed_12[index];
                    stvnl_dphi(0,2) += sum/2.0 * GlobalC::LM.DHloc_fixed_13[index];
                    stvnl_dphi(1,1) += sum/2.0 * GlobalC::LM.DHloc_fixed_22[index];
                    stvnl_dphi(1,2) += sum/2.0 * GlobalC::LM.DHloc_fixed_23[index];
                    stvnl_dphi(2,2) += sum/2.0 * GlobalC::LM.DHloc_fixed_33[index];   
                }
            }
        }
    }
    if(isstress){
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) stvnl_dphi(j,i) = stvnl_dphi(i,j);
				stvnl_dphi(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi");
    return;
}



void Force_LCAO_gamma::cal_fvnl_dbeta(
	ModuleBase::matrix& dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& fvnl_dbeta, 
	ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvnl_dbeta");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta");

	double r0[3];
	double r1[3];

    for(int iat=0; iat<GlobalC::ucell.nat; iat++)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];
        const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it].tau[ia];
        //find ajacent atom of atom ia
        //GlobalC::GridD.Find_atom( GlobalC::ucell.atoms[it].tau[ia] );
		GlobalC::GridD.Find_atom(GlobalC::ucell, GlobalC::ucell.atoms[it].tau[ia] ,it, ia);
		const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();

        //FOLLOWING ARE CONTRIBUTIONS FROM
        //VNL DUE TO PROJECTOR'S DISPLACEMENT
        for (int ad1 =0 ; ad1 < GlobalC::GridD.getAdjacentNum()+1; ad1++)
        {
            const int T1 = GlobalC::GridD.getType (ad1);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int I1 = GlobalC::GridD.getNatom (ad1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau (ad1);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            for (int ad2 =0 ; ad2 < GlobalC::GridD.getAdjacentNum()+1; ad2++)
            {
                const int T2 = GlobalC::GridD.getType (ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int I2 = GlobalC::GridD.getNatom (ad2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau (ad2);
                const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;
				if(isstress)
                {
                    r1[0] = ( tau1.x - tau0.x) ;
                    r1[1] = ( tau1.y - tau0.y) ;
                    r1[2] = ( tau1.z - tau0.z) ;
                    r0[0] = ( tau2.x - tau0.x) ;
                    r0[1] = ( tau2.y - tau0.y) ;
                    r0[2] = ( tau2.z - tau0.z) ;
                }

                if (dist1 > Rcut_Beta + Rcut_AO1
                        || dist2 > Rcut_Beta + Rcut_AO2)
                {
                    continue;
                }

                for (int jj = 0; jj < GlobalC::ucell.atoms[T1].nw; jj++)
                {
                    const int iw1_all = start1 + jj;
                    const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
                    if(mu<0) continue;
                    for (int kk = 0; kk < GlobalC::ucell.atoms[T2].nw; kk++)
                    {
                        const int iw2_all = start2 + kk;
                        const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
                        if(nu<0) continue;
                    
                        double nlm[3] = {0,0,0};
                                
                        GlobalC::UOT.snap_psibeta(
							GlobalC::ORB,
                            GlobalC::ucell.infoNL,
                            nlm, 1,
                            tau1, T1,
                            atom1->iw2l[jj], // L2
                            atom1->iw2m[jj], // m2
                            atom1->iw2n[jj], // N2
                            tau2, T2,
                            atom2->iw2l[kk], // L1
                            atom2->iw2m[kk], // m1
                            atom2->iw2n[kk], // n1
                            tau0, it, GlobalC::ucell.atoms[it].dion, GlobalV::NSPIN,
							GlobalC::ucell.atoms[it].d_so,
							GlobalC::ucell.atoms[it].non_zero_count_soc[0], // index stands for spin
							GlobalC::ucell.atoms[it].index1_soc[0],
							GlobalC::ucell.atoms[it].index2_soc[0],
							GlobalC::ucell.atoms[it].nproj_soc
							); // mohan  add 2021-05-07

                        double nlm1[3] = {0,0,0};

                        if(isstress) 
						{
								GlobalC::UOT.snap_psibeta(
								GlobalC::ORB,
                                GlobalC::ucell.infoNL,
                                nlm1, 1,
                                tau2, T2,
                                atom2->iw2l[kk], // L2
                                atom2->iw2m[kk], // m2
                                atom2->iw2n[kk], // N2
                                tau1, T1,
                                atom1->iw2l[jj], // L1
                                atom1->iw2m[jj], // m1
                                atom1->iw2n[jj], // n1
                                tau0, it, GlobalC::ucell.atoms[it].dion, GlobalV::NSPIN,
								GlobalC::ucell.atoms[it].d_so,
								GlobalC::ucell.atoms[it].non_zero_count_soc[0], // index stands for spin
								GlobalC::ucell.atoms[it].index1_soc[0],
								GlobalC::ucell.atoms[it].index2_soc[0],
								GlobalC::ucell.atoms[it].nproj_soc);
						}

                        const int index = mu * GlobalC::ParaO.ncol + nu;

                        // dbeta is minus, that's consistent.
                        // only one projector for each atom force.

                        double sum = 0.0;
                        for(int is=0; is<GlobalV::NSPIN; ++is)
                        {
                            sum += dm2d(is,index);
                        }
                        sum *= 2.0;
	
						if(isforce)
						{
							fvnl_dbeta(iat,0) -= sum * nlm[0];
							fvnl_dbeta(iat,1) -= sum * nlm[1];
							fvnl_dbeta(iat,2) -= sum * nlm[2];
						}

                        if(isstress) 
                        {
                            for(int ipol=0;ipol<3;ipol++)
							{
								// mohan update 2021-03-19
                                svnl_dbeta(0,ipol) += sum/2.0 * (nlm[0] * r0[ipol] + nlm1[0] * r1[ipol]);
                                svnl_dbeta(1,ipol) += sum/2.0 * (nlm[1] * r0[ipol] + nlm1[1] * r1[ipol]);
                                svnl_dbeta(2,ipol) += sum/2.0 * (nlm[2] * r0[ipol] + nlm1[2] * r1[ipol]);
                            }
                        }
                    }//!kk
                }//!ad2
            }//!jj
        }//!ad1
    }//!iat

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                svnl_dbeta(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta");
    return;
}


void Force_LCAO_gamma::cal_ftvnl_dphi(
	const std::vector<ModuleBase::matrix> &dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& ftvnl_dphi, 
	ModuleBase::matrix& stvnl_dphi)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_ftvnl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi");

    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        const int iat = GlobalC::ucell.iwt2iat[i];
        for(int j=0; j<GlobalV::NLOCAL; j++)
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];

            if (mu >= 0 && nu >= 0 )
            {
                const int index = mu * GlobalC::ParaO.ncol + nu;
                //contribution from deriv of AO's in T+VNL term

                double sum = 0.0;
                for(int is=0; is<GlobalV::NSPIN; ++is)
                {
                    sum += dm2d[is](nu, mu);
                }
                sum *= 2.0;

                if(isforce)
				{
					ftvnl_dphi(iat,0) += sum * GlobalC::LM.DHloc_fixed_x[index];
					ftvnl_dphi(iat,1) += sum * GlobalC::LM.DHloc_fixed_y[index];
					ftvnl_dphi(iat,2) += sum * GlobalC::LM.DHloc_fixed_z[index];
				}
                if(isstress)
                {
                    stvnl_dphi(0,0) += sum/2.0 * GlobalC::LM.DHloc_fixed_11[index];
                    stvnl_dphi(0,1) += sum/2.0 * GlobalC::LM.DHloc_fixed_12[index];
                    stvnl_dphi(0,2) += sum/2.0 * GlobalC::LM.DHloc_fixed_13[index];
                    stvnl_dphi(1,1) += sum/2.0 * GlobalC::LM.DHloc_fixed_22[index];
                    stvnl_dphi(1,2) += sum/2.0 * GlobalC::LM.DHloc_fixed_23[index];
                    stvnl_dphi(2,2) += sum/2.0 * GlobalC::LM.DHloc_fixed_33[index];   
                }
            }
        }
    }
    if(isstress){
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) stvnl_dphi(j,i) = stvnl_dphi(i,j);
				stvnl_dphi(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi");
    return;
}


void Force_LCAO_gamma::cal_fvnl_dbeta(
	const std::vector<ModuleBase::matrix> &dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& fvnl_dbeta, 
	ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvnl_dbeta");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta");

    for(int iat=0; iat<GlobalC::ucell.nat; iat++)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];
        const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it].tau[ia];
        //find ajacent atom of atom ia
        //GlobalC::GridD.Find_atom( GlobalC::ucell.atoms[it].tau[ia] );
        GlobalC::GridD.Find_atom(GlobalC::ucell, GlobalC::ucell.atoms[it].tau[ia] ,it, ia);

        //FOLLOWING ARE CONTRIBUTIONS FROM
        //VNL DUE TO PROJECTOR'S DISPLACEMENT
        for (int ad1 =0 ; ad1 < GlobalC::GridD.getAdjacentNum()+1; ad1++)
        {
            const int T1 = GlobalC::GridD.getType (ad1);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int I1 = GlobalC::GridD.getNatom (ad1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau (ad1);

            for (int ad2 =0 ; ad2 < GlobalC::GridD.getAdjacentNum()+1; ad2++)
            {
                const int T2 = GlobalC::GridD.getType (ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int I2 = GlobalC::GridD.getNatom (ad2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau (ad2);

                const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();
                const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;
                double r0[3];
				double r1[3];
				if(isstress)
                {
                    r1[0] = ( tau1.x - tau0.x) ;
                    r1[1] = ( tau1.y - tau0.y) ;
                    r1[2] = ( tau1.z - tau0.z) ;
                    r0[0] = ( tau2.x - tau0.x) ;
                    r0[1] = ( tau2.y - tau0.y) ;
                    r0[2] = ( tau2.z - tau0.z) ;
                }

                if (dist1 > Rcut_Beta + Rcut_AO1
                        || dist2 > Rcut_Beta + Rcut_AO2)
                {
                    continue;
                }

                for (int jj = 0; jj < GlobalC::ucell.atoms[T1].nw; jj++)
                {
                    const int iw1_all = start1 + jj;
                    const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
                    if(mu<0) continue;
                    for (int kk = 0; kk < GlobalC::ucell.atoms[T2].nw; kk++)
                    {
                        const int iw2_all = start2 + kk;
                        const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
                        if(nu<0) continue;

                        double nlm[3] = {0,0,0};

                        GlobalC::UOT.snap_psibeta(
										GlobalC::ORB,
                                        GlobalC::ucell.infoNL,
                                        nlm, 1,
                                        tau1, T1,
                                        atom1->iw2l[jj], // L2
                                        atom1->iw2m[jj], // m2
                                        atom1->iw2n[jj], // N2
                                        tau2, T2,
                                        atom2->iw2l[kk], // L1
                                        atom2->iw2m[kk], // m1
                                        atom2->iw2n[kk], // n1
										tau0, it, GlobalC::ucell.atoms[it].dion, GlobalV::NSPIN,
										GlobalC::ucell.atoms[it].d_so,
										GlobalC::ucell.atoms[it].non_zero_count_soc[0], // index stands for spin
										GlobalC::ucell.atoms[it].index1_soc[0],
										GlobalC::ucell.atoms[it].index2_soc[0],
										GlobalC::ucell.atoms[it].nproj_soc
								); // mohan  add 2021-05-07

                        double nlm1[3] = {0,0,0};
                        if(isstress) GlobalC::UOT.snap_psibeta(
												   GlobalC::ORB,
                                                   GlobalC::ucell.infoNL,
                                                   nlm1, 1,
                                                   tau2, T2,
                                                   atom2->iw2l[kk], // L2
                                                   atom2->iw2m[kk], // m2
                                                   atom2->iw2n[kk], // N2
                                                   tau1, T1,
                                                   atom1->iw2l[jj], // L1
                                                   atom1->iw2m[jj], // m1
                                                   atom1->iw2n[jj], // n1
                                                   tau0, it, GlobalC::ucell.atoms[it].dion, GlobalV::NSPIN,
												   GlobalC::ucell.atoms[it].d_so,
												   GlobalC::ucell.atoms[it].non_zero_count_soc[0], // index stands for spin
												   GlobalC::ucell.atoms[it].index1_soc[0],
												   GlobalC::ucell.atoms[it].index2_soc[0],
												   GlobalC::ucell.atoms[it].nproj_soc
								); // mohan  add 2021-05-07

                        //const int index = mu * GlobalC::ParaO.ncol + nu;

                        // dbeta is minus, that's consistent.
                        // only one projector for each atom force.

                        double sum = 0.0;
                        for(int is=0; is<GlobalV::NSPIN; ++is)
                        {
                            //sum += dm2d[is][index];
                            sum += dm2d[is](nu, mu);
                        }
                        sum *= 2.0;

                        if(isforce)
						{
							fvnl_dbeta(iat,0) -= sum * nlm[0];
							fvnl_dbeta(iat,1) -= sum * nlm[1];
							fvnl_dbeta(iat,2) -= sum * nlm[2];
						}

                        if(isstress) 
                        {
                            for(int ipol=0;ipol<3;ipol++)
							{
                                svnl_dbeta(0,ipol) -= sum/2.0 * (nlm[0] * r0[ipol] + nlm1[0] * r1[ipol])* -1;
                                svnl_dbeta(1,ipol) -= sum/2.0 * (nlm[1] * r0[ipol] + nlm1[1] * r1[ipol])* -1;
                                svnl_dbeta(2,ipol) -= sum/2.0 * (nlm[2] * r0[ipol] + nlm1[2] * r1[ipol])* -1;
                            }
                        }
                    }//!kk
                }//!ad2
            }//!jj
        }//!ad1
    }//!iat
    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                svnl_dbeta(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta");
    return;
}

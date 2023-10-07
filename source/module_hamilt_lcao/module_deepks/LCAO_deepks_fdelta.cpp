//wenfei 2022-1-11
//This file contains subroutines for calculating F_delta,
//which is defind as sum_mu,nu rho_mu,nu d/dX (<chi_mu|alpha>V(D)<alpha|chi_nu>)

//There are 3 subroutines in this file:
//1. cal_f_delta_gamma, which is used for gamma point calculation
//2. cal_f_delta_k, which is used for multi-k calculation
//3. check_f_delta, which prints F_delta into F_delta.dat for checking

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/vector3.h"
#include "module_base/timer.h"
#include "module_base/constants.h"

void stress_fill( 
    const double& lat0_, 
    const double& omega_,
    ModuleBase::matrix& stress_matrix)
{
    assert(omega_>0.0);
    double weight = lat0_ / omega_ ;
    for(int i=0;i<3;++i)
    {
        for(int j=0;j<3;++j)
        {
            if(j>i) stress_matrix(j,i) = stress_matrix(i,j);
            stress_matrix(i,j) *= weight ;
        }
    }
}


//force for gamma only calculations
//Pulay and HF terms are calculated together
void LCAO_Deepks::cal_f_delta_gamma(const std::vector<std::vector<double>>& dm,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const bool isstress, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_f_delta_gamma");
    this->F_delta.zero_out();

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    int nrow = this->pv->nrow;
    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int ibt = ucell.itia2iat(T2,I2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &ucell.atoms[T2];
                    const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    
                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                    const double dist2 = (tau2-tau0).norm() * ucell.lat0;

                    if (dist1 > Rcut_Alpha + Rcut_AO1
                            || dist2 > Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }					

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

                    for (int iw1=0; iw1<nw1_tot; ++iw1)
                    {
                        const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_col(iw1_all);
                        if(iw1_local < 0)continue;

                        for (int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_row(iw2_all);
                            if(iw2_local < 0)continue;

                            double nlm[3]={0,0,0};
                            double nlm_t[3] = {0,0,0}; //for stress
                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2;
                            nlm2.resize(3);

                            for(int dim=0;dim<3;dim++)
                            {
                                nlm2[dim] = this->nlm_save[iat][ad2][iw2_all][dim+1];
                            }

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < nm; ++m1)
                                    {
                                        for (int m2 = 0; m2 < nm; ++m2)
                                        {
                                            for(int dim=0;dim<3;dim++)
                                            {                                            
                                                nlm[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());

                            for(int is = 0; is < GlobalV::NSPIN; is ++)
                            {
                                // HF term is minus, only one projector for each atom force.
                                this->F_delta(iat, 0) -= 2 * dm[is][iw1_local * nrow + iw2_local] * nlm[0];
                                this->F_delta(iat, 1) -= 2 * dm[is][iw1_local * nrow + iw2_local] * nlm[1];
                                this->F_delta(iat, 2) -= 2 * dm[is][iw1_local * nrow + iw2_local] * nlm[2];

                                // Pulay term is plus, only one projector for each atom force.
                                this->F_delta(ibt, 0) += 2 * dm[is][iw1_local * nrow + iw2_local] * nlm[0];
                                this->F_delta(ibt, 1) += 2 * dm[is][iw1_local * nrow + iw2_local] * nlm[1];
                                this->F_delta(ibt, 2) += 2 * dm[is][iw1_local * nrow + iw2_local] * nlm[2];
                            }

                            if(isstress)
                            {
                                nlm1 = this->nlm_save[iat][ad2][iw2_all][0];
                                for(int i=0;i<3;i++)
                                {
                                    nlm2[i] = this->nlm_save[iat][ad1][iw1_all][i+1];
                                }

                                assert(nlm1.size()==nlm2[0].size());                                

                                int ib=0;
                                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                for(int dim=0;dim<3;dim++)
                                                {                                            
                                                    nlm_t[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                                }
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                                assert(ib==nlm1.size());

                                for(int ipol=0;ipol<3;ipol++)
                                {
                                    for(int jpol=ipol;jpol<3;jpol++)
                                    {
                                        for(int is = 0; is < GlobalV::NSPIN; is ++)
                                        {
                                            //svnl_dalpha(ipol, jpol) += dm[is](iw1_local, iw2_local) * (nlm[jpol] * r0[ipol] + nlm_t[jpol] * r1[ipol]);
                                            svnl_dalpha(ipol, jpol) += dm[is][iw1_local * nrow + iw2_local] * (nlm[jpol] * r0[ipol] + nlm_t[jpol] * r1[ipol]);
                                        }
                                    }
                                }
                            }
                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end I0
    }//end T0

    if(isstress)
    {
        stress_fill(ucell.lat0, ucell.omega, svnl_dalpha);
    }

    return;
}

//force for multi-k calculations
//Pulay and HF terms are calculated together

void LCAO_Deepks::cal_f_delta_k(const std::vector<std::vector<std::complex<double>>>& dm,/**<[in] density matrix*/
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>> &kvec_d,
    const bool isstress, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_f_delta_hf_k_new");
    ModuleBase::timer::tick("LCAO_Deepks","cal_f_delta_hf_k_new");
    this->F_delta.zero_out();

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    int nrow = this->pv->nrow;
    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);
	

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z);

                for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2,I2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &ucell.atoms[T2];
                    const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);
                    
                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                    const double dist2 = (tau2-tau0).norm() * ucell.lat0;

                    if (dist1 > Rcut_Alpha + Rcut_AO1
                            || dist2 > Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }

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

                    for (int iw1=0; iw1<nw1_tot; ++iw1)
                    {
                        const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_col(iw1_all);
                        if(iw1_local < 0)continue;

                        for (int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_row(iw2_all);
                            if(iw2_local < 0)continue;
                            double dm_current;
                            std::complex<double> tmp = 0.0;
                            for(int ik=0;ik<nks;ik++)
                            {
                                const double arg = - ( kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                //tmp += dm[ik](iw1_local, iw2_local) * kphase;
                                tmp += dm[ik][iw1_local * nrow + iw2_local] * kphase;
                            }
                            dm_current=tmp.real();

                            double nlm[3]={0,0,0};
                            double nlm_t[3] = {0,0,0}; //for stress
                            key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                            key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
                            std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2;
                            nlm2.resize(3);
                            for(int dim=0;dim<3;dim++)
                            {
                                nlm2[dim] = this->nlm_save_k[iat][key_2][iw2_all][dim+1];
                            }

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < nm; ++m1)
                                    {
                                        for (int m2 = 0; m2 < nm; ++m2)
                                        {
                                            for(int dim=0;dim<3;dim++)
                                            {                                            
                                                nlm[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());

                            // Pulay term is plus
                            this->F_delta(ibt2, 0) += 2.0 * dm_current * nlm[0];
                            this->F_delta(ibt2, 1) += 2.0 * dm_current * nlm[1];
                            this->F_delta(ibt2, 2) += 2.0 * dm_current * nlm[2];

                            // HF term is minus, only one projector for each atom force.
                            this->F_delta(iat, 0) -= 2.0 * dm_current * nlm[0];
                            this->F_delta(iat, 1) -= 2.0 * dm_current * nlm[1];
                            this->F_delta(iat, 2) -= 2.0 * dm_current * nlm[2];

                            if(isstress)
                            {
                                nlm1 = this->nlm_save_k[iat][key_2][iw2_all][0];
                                for(int i=0;i<3;i++)
                                {
                                    nlm2[i] = this->nlm_save_k[iat][key_1][iw1_all][i+1];
                                }

                                assert(nlm1.size()==nlm2[0].size());                                

                                int ib=0;
                                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                for(int dim=0;dim<3;dim++)
                                                {                                            
                                                    nlm_t[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                                }
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                                assert(ib==nlm1.size());
                                   
                                for(int ipol=0;ipol<3;ipol++)
                                {
                                    svnl_dalpha(0,ipol) -= dm_current * (nlm[0] * r0[ipol] + nlm_t[0] * r1[ipol])* -1.0;
                                    svnl_dalpha(1,ipol) -= dm_current * (nlm[1] * r0[ipol] + nlm_t[1] * r1[ipol])* -1.0;
                                    svnl_dalpha(2,ipol) -= dm_current * (nlm[2] * r0[ipol] + nlm_t[2] * r1[ipol])* -1.0;
                                }

                            }
                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end I0
    }//end T0

    if(isstress)
    {
        stress_fill(ucell.lat0, ucell.omega, svnl_dalpha);
    }
    ModuleBase::timer::tick("LCAO_Deepks","cal_f_delta_hf_k_new");
    return;
}

//prints F_delta into F_delta.dat
void LCAO_Deepks::check_f_delta(const int nat, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Deepks", "check_F_delta");

    std::ofstream ofs("F_delta.dat");
    ofs<<std::setprecision(10);

    for (int iat=0; iat<nat; iat++)
    {
        ofs << F_delta(iat,0) << " " << F_delta(iat,1) << " " << F_delta(iat,2) << std::endl;
    }

    std::ofstream ofs1("stress_delta.dat");
    ofs1<<std::setprecision(10);
    for (int ipol=0; ipol<3; ipol++)
    {
        for (int jpol=0; jpol<3; jpol++)
        {
            ofs1 << svnl_dalpha(ipol,jpol) << " ";
        }
        ofs1 << std::endl;
    }    
    return;
}

#endif
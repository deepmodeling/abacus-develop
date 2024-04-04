#include "sparse_format.h"

void sparse_format::cal_dH(
		LCAO_Matrix &lm,
		LCAO_gen_fixedH &gen_h, 
		const int &current_spin, 
		const double &sparse_threshold,
		Gint_k &gint_k)
{
    ModuleBase::TITLE("sparse_format","cal_dH");

    sparse_format::set_R_range(lm);

    const int nnr = lm.ParaV->nnr;
    lm.DHloc_fixedR_x = new double[nnr];
    lm.DHloc_fixedR_y = new double[nnr];
    lm.DHloc_fixedR_z = new double[nnr];

    ModuleBase::GlobalFunc::ZEROS(lm.DHloc_fixedR_x, lm.ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(lm.DHloc_fixedR_y, lm.ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(lm.DHloc_fixedR_z, lm.ParaV->nloc);
    // cal dT=<phi|kin|dphi> in LCAO
    // cal T + VNL(P1) in LCAO basis
    if(GlobalV::CAL_STRESS)
	{
        GlobalV::CAL_STRESS = false;

        gen_h.build_ST_new('T', true, GlobalC::ucell, lm.Hloc_fixedR.data());

        GlobalV::CAL_STRESS = true;
    }
    else
    {
        gen_h.build_ST_new('T', true, GlobalC::ucell, lm.Hloc_fixedR.data());
    }
    gen_h.build_Nonlocal_mu_new (lm.Hloc_fixed.data(), true);
    
    sparse_format::cal_dSTN_R(lm, current_spin, sparse_threshold);

    delete[] lm.DHloc_fixedR_x;
    delete[] lm.DHloc_fixedR_y;
    delete[] lm.DHloc_fixedR_z;

    gint_k.cal_dvlocal_R_sparseMatrix(current_spin, sparse_threshold, &lm);
}


void sparse_format::set_R_range(LCAO_Matrix &lm)
{
    int R_minX = int(GlobalC::GridD.getD_minX());
    int R_minY = int(GlobalC::GridD.getD_minY());
    int R_minZ = int(GlobalC::GridD.getD_minZ());

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    for(int ix = 0; ix < R_x; ix++)
    {
        for(int iy = 0; iy < R_y; iy++)
        {
            for(int iz = 0; iz < R_z; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix+R_minX, iy+R_minY, iz+R_minZ);
                lm.all_R_coor.insert(temp_R);
            }
        }
    }

    return;
}


void sparse_format::cal_dSTN_R(
		LCAO_Matrix &lm,
		const int &current_spin, 
		const double &sparse_threshold)
{
    ModuleBase::TITLE("sparse_format","cal_dSTN_R");

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double temp_value_double;
    std::complex<double> temp_value_complex;

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

					Abfs::Vector3_Order<int> dR(
							GlobalC::GridD.getBox(ad).x, 
							GlobalC::GridD.getBox(ad).y, 
							GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = lm.ParaV->global2local_row(iw1_all);

						if(mu<0)
						{
							continue;
						}

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = lm.ParaV->global2local_col(iw2_all);

							if(nu<0)
							{
								continue;
							}

                            if(GlobalV::NSPIN!=4)
                            {
                                temp_value_double = lm.DHloc_fixedR_x[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    lm.dHRx_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = lm.DHloc_fixedR_y[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    lm.dHRy_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = lm.DHloc_fixedR_z[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    lm.dHRz_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                ModuleBase::WARNING_QUIT("cal_dSTN_R","nspin=4 with SOC is not supported yet.");
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

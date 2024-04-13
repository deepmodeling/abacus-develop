#include "spar_dh.h"

void sparse_format::cal_dH(
		LCAO_Matrix &lm,
		Grid_Driver &grid,
		LCAO_gen_fixedH &gen_h, 
		const int &current_spin, 
		const double &sparse_thr,
		Gint_k &gint_k)
{
    ModuleBase::TITLE("sparse_format","cal_dH");

    sparse_format::set_R_range(lm.all_R_coor, grid);

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

        gen_h.build_ST_new('T', true, GlobalC::ucell, GlobalC::ORB, GlobalC::UOT, &(GlobalC::GridD), lm.Hloc_fixedR.data());

        GlobalV::CAL_STRESS = true;
    }
    else
    {
        gen_h.build_ST_new('T', true, GlobalC::ucell, GlobalC::ORB, GlobalC::UOT, &(GlobalC::GridD), lm.Hloc_fixedR.data());
    }
    gen_h.build_Nonlocal_mu_new (lm.Hloc_fixed.data(), true, GlobalC::ucell, GlobalC::ORB, GlobalC::UOT, &(GlobalC::GridD));
    
    sparse_format::cal_dSTN_R(lm, grid, current_spin, sparse_thr);

    delete[] lm.DHloc_fixedR_x;
    delete[] lm.DHloc_fixedR_y;
    delete[] lm.DHloc_fixedR_z;

    gint_k.cal_dvlocal_R_sparseMatrix(current_spin, sparse_thr, &lm);

    return;
}


void sparse_format::set_R_range(
        std::set<Abfs::Vector3_Order<int>> &all_R_coor,
		Grid_Driver &grid)
{
    const int RminX = int(grid.getD_minX());
    const int RminY = int(grid.getD_minY());
    const int RminZ = int(grid.getD_minZ());

    const int Rx = grid.getCellX();
    const int Ry = grid.getCellY();
    const int Rz = grid.getCellZ();

    for(int ix = 0; ix < Rx; ix++)
    {
        for(int iy = 0; iy < Ry; iy++)
        {
            for(int iz = 0; iz < Rz; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix+RminX, iy+RminY, iz+RminZ);
                all_R_coor.insert(temp_R);
            }
        }
    }

    return;
}


void sparse_format::cal_dSTN_R(
		LCAO_Matrix &lm,
		Grid_Driver &grid,
		const int &current_spin, 
		const double &sparse_thr)
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
            grid.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < grid.getAdjacentNum()+1; ++ad)
            {
                const int T2 = grid.getType(ad);
                const int I2 = grid.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = grid.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
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
							grid.getBox(ad).x, 
							grid.getBox(ad).y, 
							grid.getBox(ad).z);

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
                                if (std::abs(temp_value_double) > sparse_thr)
                                {
                                    lm.dHRx_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = lm.DHloc_fixedR_y[index];
                                if (std::abs(temp_value_double) > sparse_thr)
                                {
                                    lm.dHRy_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = lm.DHloc_fixedR_z[index];
                                if (std::abs(temp_value_double) > sparse_thr)
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

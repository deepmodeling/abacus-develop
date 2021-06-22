#include "../src_pw/global.h"
#include "LCAO_hamilt.h"
#include "build_st_pw.h"
#include "../src_global/sltk_atom_arrange.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "dftu.h"

LCAO_Hamilt::LCAO_Hamilt()
{ 
    init_s = false;
}

LCAO_Hamilt::~LCAO_Hamilt()
{
	if(test_deconstructor)
	{
		cout << " ~LCAO_Hamilt()" << endl;	
	}
}

//--------------------------------------------
// 'calculate_STNR_gamma' or 
// 'calculate_STNR_k' functions are called
//--------------------------------------------
void LCAO_Hamilt::set_lcao_matrices(void)
{
	TITLE("LCAO_Hamilt","set_lcao_matrices");
	timer::tick("LCAO_Hamilt","set_lcao_matrices",'E');

	if(GAMMA_ONLY_LOCAL)
	{
		// mohan add 2012-03-29
		// calculate the grid integration of 'Vl' matrix for gamma algorithms.
		this->GG.prepare(ucell.latvec, ucell.lat0);
	
		// calulate the 'S', 'T' and 'Vnl' matrix for gamma algorithms.
	    this->calculate_STNR_gamma();	

	}
	else // multiple k-points
	{
		// calculate the 'S', 'T' and 'Vnl' matrix for k-points algorithms.
		this->calculate_STNR_k();

		// calculate the grid integration of 'Vl' matrix for l-points algorithms.
		this->GK.init(pw.nbx, pw.nby, pw.nbzp, pw.nbzp_start, pw.ncxyz);

	}

	// initial the overlap matrix is done.	
    this->init_s = true;
	//cout << " init_s=" << init_s << endl; //delete 2015-09-06, xiaohui
//	OUT(ofs_running,"init_s",init_s);

	timer::tick("LCAO_Hamilt","set_lcao_matrices",'E');
	return;
}

void LCAO_Hamilt::calculate_Hgamma( const int &ik )				// Peize Lin add ik 2016-12-03
{
	TITLE("LCAO_Hamilt","calculate_Hgamma");
	timer::tick("LCAO_Hamilt","cal_Hgamma",'F');

	// Set the matrix 'H' to zero.
	LM.zeros_HSgamma('H'); // 3 stands for Hloc.

	bool local_pw = false;

	if(local_pw)
	{
		cout << "\n Call build_H in plane wave basis!" << endl;
	 	// Use plane wave basis to calculate 'Vl' matrix. 
		Build_ST_pw bsp;
		// 0 stands for, 0 stands for k point.
		bsp.set_local(0);
	}
	else
	{
		time_t time_vlocal_start = time(NULL);
	
		// calculate the 'Vl' matrix using gamma-algorithms.
		if(VL_IN_H)
		{	
			this->GG.cal_vlocal(pot.vr_eff1);

			// Peize Lin add 2016-12-03
			if( 5==xcf.iexch_now && 0==xcf.igcx_now )				// HF
			{
				exx_lcao.add_Hexx(ik,1);
			}
			else if( 6==xcf.iexch_now && 8==xcf.igcx_now )			// PBE0
			{
				exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha);
			}
			else if( 9==xcf.iexch_now && 12==xcf.igcx_now )			// HSE
			{
				exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha);
			}
		}

		time_t time_vlocal_end = time(NULL);
		OUT_TIME("vlocal integration",time_vlocal_start,time_vlocal_end);
	}

	//add T+VNL+Vl matrix.
	LM.update_Hloc();


	//test
	if(NURSE)
	{
		LM.print_HSgamma('S'); // S
		LM.print_HSgamma('T');
		LM.print_HSgamma('H');
	//	WARNING_QUIT("LCAO_Hamilt::calculate_Hgamma","print the H,S matrix");
//		QUIT();
	}


	timer::tick("LCAO_Hamilt","cal_Hgamma",'F');
	return;
}



void LCAO_Hamilt::calculate_STNR_gamma(void)
{
	TITLE("LCAO_Hamilt","calculate_fixed");

	OUT(ofs_running,"gamma_only_local",GAMMA_ONLY_LOCAL);

	// must be done after "setup_this_ion_iter"
	// because some basic parameters should be initialized
	// in UHM.GG.init();

	LM.zeros_HSgamma('S');    	

	this->genH.calculate_S_no();	

	//LM.print_HSgamma('S');

	//-------------------------------------
	// test using plane wave calculations.
	// all the matrixs are stored in LM.
	// LM.allocate_HS_k(ParaO.nloc);
	// Build_ST_pw bsp;
	// bsp.set_ST(0, 'S');
	// LM.print_HSk('S','R',1.0e-5);
	//-------------------------------------

	// set T and Vnl matrix to zero.
	// 2 stands for LM.Hloc_fixed matrix.
	LM.zeros_HSgamma('T'); 

	//add nonlocal pseudopotential matrix element
	time_t time_vnl_start = time(NULL);
	if(VNL_IN_H)
	{
		genH.calculate_NL_no();
	}
	time_t time_vnl_end = time(NULL);

//	OUT(ofs_running, "Time to calculate <psi|Vnl|psi>", std::difftime(time_vnl_end, time_vnl_start));
	
	//add kinetic energy matrix element
	time_t time_t_start = time(NULL);
	if(T_IN_H)
	{
		genH.calculate_T_no();
//		LM.print_HSgamma('T');
	}
	time_t time_t_end = time(NULL);

	//	ofs_running << " T+Vnl matrix" << endl;
	//LM.print_HSgamma('T');

	OUT_TIME("kinetical matrix",time_t_start, time_t_end);
	OUT_TIME("vnl matrix",time_vnl_start, time_vnl_end);

	return;
}


#include "LCAO_nnr.h"
// be called in LOOP_elec::cal_bands(). 
void LCAO_Hamilt::calculate_Hk(const int &ik)
{
	TITLE("LCAO_Hamilt","calculate_Hk");
	timer::tick("LCAO_Hamilt","calculate_Hk",'F');

	// whether you want to calculate the local potential
	// or not, you need to set this matrix to 0.
	LM.zeros_HSk('H');

	if(VL_IN_H)
	{
		//-------------------------
		// set the local potential
		// in plane wave basis.
		//-------------------------
//		Build_ST_pw bsp;
//		bsp.set_local(ik);	
//		LM.print_HSk('H','C',1.0e-5);

		//--------------------------
		// set the local potential
		// in LCAO basis.
		//--------------------------
		LM.zeros_HSR('H', LNNR.nnr);

		if(NSPIN!=4) 
		{
			this->GK.folding_vl_k(ik);
		}
		else 
		{
			this->GK.folding_vl_k_nc(ik);
		}

		// Peize Lin add 2016-12-03
		if( 5==xcf.iexch_now && 0==xcf.igcx_now )				// HF
		{
			exx_lcao.add_Hexx(ik,1);
		}
		else if( 6==xcf.iexch_now && 8==xcf.igcx_now )			// PBE0
		{
			exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha);
		}
		else if( 9==xcf.iexch_now && 12==xcf.igcx_now )			// HSE
		{
			exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha);
		}
	}


	//-----------------------------------------
	// folding matrix here: S(k) (SlocR->Sloc2)
	// folding matrix here: T(k)+Vnl(k)
	// (Hloc_fixed->Hloc_fixed2)
	//-----------------------------------------
	LM.zeros_HSk('S');
	LM.zeros_HSk('T');
//	cout << " after folding Hfixed k." << endl;
	LNNR.folding_fixedH(ik);

	//------------------------------------------
	// Add T(k)+Vnl(k)+Vlocal(k)
	// (Hloc2 += Hloc_fixed2), (complex matrix)
	//------------------------------------------
//	cout << " Folding matrix here." << endl;
	LM.update_Hloc2();

/*
	if(NURSE)
	{
		LM.print_HSk('H','R',1.0e-5);
//		LM.print_HSk('S','R',1.0e-5);
	}
	*/
	
	timer::tick("LCAO_Hamilt","calculate_Hk",'F');
	return;
}

// only need to do the first time.
// available for all k points.
void LCAO_Hamilt::calculate_STNR_k(void)
{
    TITLE("Hamilt_Linear","calculate_STBR_k");

	//--------------------------------------------
	// set S(R) to zero.
	// the total value of S(R) in this processor
	// is LNNR.nnr.
	// and store in LM.SlocR.
	//--------------------------------------------
	LM.zeros_HSR('S', LNNR.nnr);
    this->genH.calculate_S_no();	

	//------------------------------
	// set T(R) and Vnl(R) to zero.
	// and then calculate it
	// and store in LM.Hloc_fixedR.
	//------------------------------
	LM.zeros_HSR('T', LNNR.nnr);
	


	if(T_IN_H)
	{
		this->genH.calculate_T_no();	
	}


	if(VNL_IN_H)
	{
		this->genH.calculate_NL_no();
	}


	return;

	//-----------------------------------
	// this part is used for checking
	// the consistent between LCAO
	// and plane wave basis.
	//-----------------------------------
	
	// check in plane wave basis.	
	Build_ST_pw bsp;
	for(int ik=0; ik<kv.nks; ik++)
	{
		cout << " ik=" << ik << " ------------------------------------------" << endl;
	
		//----------------------------------------------
		// Check the matrix in plane wave basis.
		//----------------------------------------------
		LM.zeros_HSk('S');
		LM.zeros_HSk('T');

		bsp.set_ST(ik, 'S');
		bsp.set_ST(ik, 'T');

		cout << " --> PW S" << endl;
		LM.print_HSk('S','R',1.0e-5);
		cout << " --> PW T" << endl;
		LM.print_HSk('T','R',1.0e-5);

		string fn = "Sloc2pw.dat";
		LM.output_HSk('S', fn);
		
		//------------------------------------------
		// folding the SlocR and Hloc_fixedR matrix
		// into Sloc2 and Hloc_fixed2 matrix.
		//------------------------------------------
		LM.zeros_HSk('S');
		LM.zeros_HSk('T');
		LNNR.folding_fixedH(ik);
		cout << " --> LCAO S" << endl;
		LM.print_HSk('S','R',1.0e-5);	
		cout << " --> LCAO T+Vnl" << endl;
		LM.print_HSk('T','R',1.0e-5);	

		string fn2 = "Sloc2lcao.dat";
		LM.output_HSk('S',fn2);

		//----------------
		// test gamma Vnl	
		//----------------
//		GAMMA_ONLY_LOCAL = true;
//		LM.allocate_HS_gamma(ParaO.nloc);
//		LM.zeros_HSgamma('H');
//		UHM.genH.calculate_NL_no( nstart );
//		GAMMA_ONLY_LOCAL = false;
//		cout << " Correct LCAO Vnl " << endl;
//		LM.print_HSgamma('H');		
//		UHM.genH.calculate_NL_no( nstart );
//		GAMMA_ONLY_LOCAL = false;
//		cout << " Correct LCAO Vnl " << endl;
//		LM.print_HSgamma('H');		
		
	}
	
    return;
}


void LCAO_Hamilt::calculate_STN_R(void)
{
    TITLE("LCAO_Hamilt","calculate_STN_R");

    //int iat = 0;
    int index = 0;
    Vector3<double> dtau, tau1, tau2;
    Vector3<double> dtau1, dtau2, tau0;

    LM.allocate_Hloc_fixedR_tr();
    LM.allocate_HR_tr();
    LM.allocate_SlocR_tr();

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    for(int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GridD.Find_atom(tau1);
            GridD.Find_atom(ucell, tau1, T1, I1);
            Atom* atom1 = &ucell.atoms[T1];
            const int start = ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GridD.getType(ad);
                const int I2 = GridD.getNatom(ad);
                Atom* atom2 = &ucell.atoms[T2];

                tau2 = GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GridD.getType(ad0);
                        //const int I0 = GridD.getNatom(ad0);
                        //const int iat0 = ucell.itia2iat(T0, I0);
                        //const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
                        double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

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

                    Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z);
                    R_x = (int) (dR.x - R_minX);
                    R_y = (int) (dR.y - R_minY);
                    R_z = (int) (dR.z - R_minZ);

                    for(int ii=0; ii<atom1->nw*NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = ParaO.trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = ParaO.trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            int iic;
                            if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic=mu+nu*ParaO.nrow;
                            }
                            else
                            {
                                iic=mu*ParaO.ncol+nu;
                            }

                            if(NSPIN!=4)
                            {
                                LM.SlocR_tr[R_x][R_y][R_z][iic] = LM.SlocR[index];
                                LM.Hloc_fixedR_tr[R_x][R_y][R_z][iic] = LM.Hloc_fixedR[index];
                            }
                            else
                            {
                                LM.SlocR_tr_soc[R_x][R_y][R_z][iic] = LM.SlocR_soc[index];
                                LM.Hloc_fixedR_tr_soc[R_x][R_y][R_z][iic] = LM.Hloc_fixedR_soc[index];
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

void LCAO_Hamilt::calculate_STN_R_sparse(const double &sparse_threshold)
{
    TITLE("LCAO_Hamilt","calculate_STN_R_sparse");

    //int iat = 0;
    int index = 0;
    Vector3<double> dtau, tau1, tau2;
    Vector3<double> dtau1, dtau2, tau0;

    LM.allocate_HS_R_sparse();

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    for(int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GridD.Find_atom(tau1);
            GridD.Find_atom(ucell, tau1, T1, I1);
            Atom* atom1 = &ucell.atoms[T1];
            const int start = ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GridD.getType(ad);
                const int I2 = GridD.getNatom(ad);
                Atom* atom2 = &ucell.atoms[T2];

                tau2 = GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GridD.getType(ad0);
                        //const int I0 = GridD.getNatom(ad0);
                        //const int iat0 = ucell.itia2iat(T0, I0);
                        //const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
                        double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

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

                    Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z);
                    R_x = (int) (dR.x - R_minX);
                    R_y = (int) (dR.y - R_minY);
                    R_z = (int) (dR.z - R_minZ);

                    for(int ii=0; ii<atom1->nw*NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = ParaO.trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = ParaO.trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            if(NSPIN!=4)
                            {
								double temp_value = LM.SlocR[index];
								if (abs(temp_value) > sparse_threshold)
								{
									LM.SR_sparse[R_x][R_y][R_z][iw1_all].insert(pair<size_t, double>(iw2_all, temp_value));
								}

								temp_value = LM.Hloc_fixedR[index];
								if (abs(temp_value) > sparse_threshold)
								{
									LM.HR_sparse[R_x][R_y][R_z][iw1_all].insert(pair<size_t, double>(iw2_all, temp_value));
								}
                            }
                            else
                            {
								complex<double> temp_value = LM.SlocR_soc[index];
								if(abs(temp_value) > sparse_threshold)
								{
									LM.SR_soc_sparse[R_x][R_y][R_z][iw1_all].insert(pair<size_t, complex<double>>(iw2_all, temp_value));
								}

								temp_value = LM.Hloc_fixedR_soc[index];
								if(abs(temp_value) > sparse_threshold)
								{
									LM.HR_soc_sparse[R_x][R_y][R_z][iw1_all].insert(pair<size_t, complex<double>>(iw2_all, temp_value));
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


void LCAO_Hamilt::calculate_HSR_sparse(const int &current_spin, const double &sparse_threshold)
{
	TITLE("LCAO_Hamilt","calculate_HSR_sparse");

	calculate_STN_R_sparse(sparse_threshold);

	GK.cal_vlocal_R_sparseMatrix(current_spin, sparse_threshold);

	if (INPUT.dft_plus_u)
	{
		if (NSPIN == 4)
		{
			calculat_HR_dftu_soc_sparse(current_spin, sparse_threshold);
		}
		else
		{
			calculat_HR_dftu_sparse(current_spin, sparse_threshold);
		}
	}

}

void LCAO_Hamilt::calculat_HR_dftu_sparse(const int &current_spin, const double &sparse_threshold)
{
	TITLE("LCAO_Hamilt","calculat_HR_dftu_sparse");
	timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

	int R_x = GridD.getCellX();
    int R_y = GridD.getCellY();
    int R_z = GridD.getCellZ();

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

	double *HR_tmp = new double[ParaO.nloc];
	double *SR_tmp = new double[ParaO.nloc];

	int ir;
	int ic;
	int iic;

	for(int ix=0; ix<R_x; ix++)
    {
        for(int iy=0; iy<R_y; iy++)
        {
            for(int iz=0; iz<R_z; iz++)
            {	
				map<size_t, map<size_t, double>> &temp_HR_sparse = LM.HR_sparse[ix][iy][iz];
				map<size_t, map<size_t, double>> &temp_SR_sparse = LM.SR_sparse[ix][iy][iz];

				ZEROS(HR_tmp, ParaO.nloc);
				ZEROS(SR_tmp, ParaO.nloc);
				
				for (auto &iter : temp_SR_sparse)
				{
					ir = ParaO.trace_loc_row[iter.first];
					for (auto &value : iter.second)
					{
						ic = ParaO.trace_loc_col[value.first];
						if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
						{
							iic = ir + ic * ParaO.nrow;
						}
						else
						{
							iic = ir * ParaO.ncol + ic;
						}
						SR_tmp[iic] = value.second;
					}
				}

				dftu.cal_eff_pot_mat_R_double(current_spin, SR_tmp, HR_tmp);

				for (int i = 0; i < NLOCAL; ++i)
				{
					ir = ParaO.trace_loc_row[i];
					if (ir >= 0)
					{
						for (int j = 0; j < NLOCAL; ++j)
						{
							ic = ParaO.trace_loc_col[j];
							if (ic >= 0)
							{
								if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
								{
									iic = ir + ic * ParaO.nrow;
								}
								else
								{
									iic = ir * ParaO.ncol + ic;
								}

								if (abs(HR_tmp[iic]) > sparse_threshold)
								{
									double &value = temp_HR_sparse[i][j];
									value += HR_tmp[iic];
									if (abs(value) < sparse_threshold)
									{
										temp_HR_sparse[i].erase(j);
									}
								}
							}
						}
					}
				}

			}
		}
	}

	delete[] HR_tmp;
	delete[] SR_tmp;
	HR_tmp = nullptr;
	SR_tmp = nullptr;

	timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

}

void LCAO_Hamilt::calculat_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold)
{
	TITLE("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");
	timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

	int R_x = GridD.getCellX();
    int R_y = GridD.getCellY();
    int R_z = GridD.getCellZ();

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

	complex<double> *HR_soc_tmp = new complex<double>[ParaO.nloc];
	complex<double> *SR_soc_tmp = new complex<double>[ParaO.nloc];

	int ir;
	int ic;
	int iic;

	for(int ix=0; ix<R_x; ix++)
    {
        for(int iy=0; iy<R_y; iy++)
        {
            for(int iz=0; iz<R_z; iz++)
            {
				map<size_t, map<size_t, complex<double>>> &temp_HR_soc_sparse = LM.HR_soc_sparse[ix][iy][iz];
				map<size_t, map<size_t, complex<double>>> &temp_SR_soc_sparse = LM.SR_soc_sparse[ix][iy][iz];

				ZEROS(HR_soc_tmp, ParaO.nloc);
				ZEROS(SR_soc_tmp, ParaO.nloc);
				
				for (auto &iter : temp_SR_soc_sparse)
				{
					ir = ParaO.trace_loc_row[iter.first];
					for (auto &value : iter.second)
					{
						ic = ParaO.trace_loc_col[value.first];
						if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
						{
							iic = ir + ic * ParaO.nrow;
						}
						else
						{
							iic = ir * ParaO.ncol + ic;
						}
						SR_soc_tmp[iic] = value.second;
					}
				}

				dftu.cal_eff_pot_mat_R_complex_double(current_spin, SR_soc_tmp, HR_soc_tmp);

				for (int i = 0; i < NLOCAL; ++i)
				{
					ir = ParaO.trace_loc_row[i];
					if (ir >= 0)
					{
						for (int j = 0; j < NLOCAL; ++j)
						{
							ic = ParaO.trace_loc_col[j];
							if (ic >= 0)
							{
								if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
								{
									iic = ir + ic * ParaO.nrow;
								}
								else
								{
									iic = ir * ParaO.ncol + ic;
								}

								if (abs(HR_soc_tmp[iic]) > sparse_threshold)
								{
									complex<double> &value = temp_HR_soc_sparse[i][j];
									value += HR_soc_tmp[iic];
									if (abs(value) < sparse_threshold)
									{
										temp_HR_soc_sparse[i].erase(j);
									}
								}
							}
						}
					}
				}

			}
		}
	}

	delete[] HR_soc_tmp;
	delete[] SR_soc_tmp;
	HR_soc_tmp = nullptr;
	SR_soc_tmp = nullptr;
	
	timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

}

void LCAO_Hamilt::destroy_all_HSR_sparse(void)
{
	LM.destroy_HS_R_sparse();
}

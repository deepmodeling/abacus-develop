#include "../src_pw/global.h"
#include "LCAO_hamilt.h"
#include "build_st_pw.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "dftu.h"
#ifdef __DEEPKS
#include "LCAO_descriptor.h"	//caoyu add 2021-07-26
#endif

LCAO_Hamilt::LCAO_Hamilt()
{ 
    init_s = false;
}

LCAO_Hamilt::~LCAO_Hamilt()
{
	if(GlobalV::test_deconstructor)
	{
		std::cout << " ~LCAO_Hamilt()" << std::endl;	
	}
}

//--------------------------------------------
// 'calculate_STNR_gamma' or 
// 'calculate_STNR_k' functions are called
//--------------------------------------------
void LCAO_Hamilt::set_lcao_matrices(void)
{
	ModuleBase::TITLE("LCAO_Hamilt","set_lcao_matrices");
	ModuleBase::timer::tick("LCAO_Hamilt","set_lcao_matrices");

	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		// mohan add 2012-03-29
		// calculate the grid integration of 'Vl' matrix for gamma algorithms.
		this->GG.prepare(GlobalC::ucell.latvec, GlobalC::ucell.lat0);
	
		// calulate the 'S', 'T' and 'Vnl' matrix for gamma algorithms.
	    this->calculate_STNR_gamma();	

	}
	else // multiple k-points
	{
		// calculate the 'S', 'T' and 'Vnl' matrix for k-points algorithms.
		this->calculate_STNR_k();

		// calculate the grid integration of 'Vl' matrix for l-points algorithms.
		this->GK.init(GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbzp, GlobalC::pw.nbzp_start, GlobalC::pw.ncxyz);

	}

	// initial the overlap matrix is done.	
    this->init_s = true;
	//std::cout << " init_s=" << init_s << std::endl; //delete 2015-09-06, xiaohui
//	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"init_s",init_s);

	ModuleBase::timer::tick("LCAO_Hamilt","set_lcao_matrices");
	return;
}

void LCAO_Hamilt::calculate_Hgamma( const int &ik )				// Peize Lin add ik 2016-12-03
{
	ModuleBase::TITLE("LCAO_Hamilt","calculate_Hgamma");
	ModuleBase::timer::tick("LCAO_Hamilt","cal_Hgamma");

	// Set the matrix 'H' to zero.
	GlobalC::LM.zeros_HSgamma('H'); // 3 stands for Hloc.

	bool local_pw = false;

	if(local_pw)
	{
		std::cout << "\n Call build_H in plane wave basis!" << std::endl;
	 	// Use plane wave basis to calculate 'Vl' matrix. 
		Build_ST_pw bsp;
		// 0 stands for, 0 stands for k point.
		bsp.set_local(0);
	}
	else
	{
		time_t time_vlocal_start = time(NULL);
	
		// calculate the 'Vl' matrix using gamma-algorithms.
		if(GlobalV::VL_IN_H)
		{	
			this->GG.cal_vlocal(GlobalC::pot.vr_eff1);

			// Peize Lin add 2016-12-03
			if( 5==GlobalC::xcf.iexch_now && 0==GlobalC::xcf.igcx_now )				// HF
			{
				GlobalC::exx_lcao.add_Hexx(ik,1);
			}
			else if( 6==GlobalC::xcf.iexch_now && 8==GlobalC::xcf.igcx_now )			// PBE0
			{
				GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha);
			}
			else if( 9==GlobalC::xcf.iexch_now && 12==GlobalC::xcf.igcx_now )			// HSE
			{
				GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha);
			}
		}

		time_t time_vlocal_end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("vlocal integration",time_vlocal_start,time_vlocal_end);
	}
	
#ifdef __DEEPKS	//caoyu add 2021-07-26 for DeePKS

	if (INPUT.deepks_scf)
	{
		ld.cal_v_delta(LOC.wfc_dm_2d.dm_gamma[0]);
		ld.add_v_delta();
	}
	
#endif
	
	//add T+VNL+Vl matrix.
	GlobalC::LM.update_Hloc();


	//test
	if(GlobalV::NURSE)
	{
		GlobalC::LM.print_HSgamma('S'); // S
		GlobalC::LM.print_HSgamma('T');
		GlobalC::LM.print_HSgamma('H');
	//	ModuleBase::WARNING_QUIT("LCAO_Hamilt::calculate_Hgamma","print the H,S matrix");
//		ModuleBase::QUIT();
	}


	ModuleBase::timer::tick("LCAO_Hamilt","cal_Hgamma");
	return;
}



void LCAO_Hamilt::calculate_STNR_gamma(void)
{
	ModuleBase::TITLE("LCAO_Hamilt","calculate_fixed");

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"gamma_only_local",GlobalV::GAMMA_ONLY_LOCAL);

	// must be done after "setup_this_ion_iter"
	// because some basic parameters should be initialized
	// in GlobalC::UHM.GG.init();

	GlobalC::LM.zeros_HSgamma('S');    	

	this->genH.calculate_S_no();	

	//GlobalC::LM.print_HSgamma('S');

	//-------------------------------------
	// test using plane wave calculations.
	// all the matrixs are stored in GlobalC::LM.
	// GlobalC::LM.allocate_HS_k(GlobalC::ParaO.nloc);
	// Build_ST_pw bsp;
	// bsp.set_ST(0, 'S');
	// GlobalC::LM.print_HSk('S','R',1.0e-5);
	//-------------------------------------

	// set T and Vnl matrix to zero.
	// 2 stands for GlobalC::LM.Hloc_fixed matrix.
	GlobalC::LM.zeros_HSgamma('T'); 

	//add nonlocal pseudopotential matrix element
	time_t time_vnl_start = time(NULL);
	if(GlobalV::VNL_IN_H)
	{
		genH.calculate_NL_no();
	}
	time_t time_vnl_end = time(NULL);

//	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Time to calculate <psi|Vnl|psi>", std::difftime(time_vnl_end, time_vnl_start));
	
	//add kinetic energy matrix element
	time_t time_t_start = time(NULL);
	if(GlobalV::T_IN_H)
	{
		genH.calculate_T_no();
//		GlobalC::LM.print_HSgamma('T');
	}
	time_t time_t_end = time(NULL);

	//	GlobalV::ofs_running << " T+Vnl matrix" << std::endl;
	//GlobalC::LM.print_HSgamma('T');

	ModuleBase::GlobalFunc::OUT_TIME("kinetical matrix",time_t_start, time_t_end);
	ModuleBase::GlobalFunc::OUT_TIME("vnl matrix",time_vnl_start, time_vnl_end);

	return;
}


#include "LCAO_nnr.h"
// be called in LOOP_elec::cal_bands(). 
void LCAO_Hamilt::calculate_Hk(const int &ik)
{
	ModuleBase::TITLE("LCAO_Hamilt","calculate_Hk");
	ModuleBase::timer::tick("LCAO_Hamilt","calculate_Hk");

	// whether you want to calculate the local potential
	// or not, you need to set this matrix to 0.
	GlobalC::LM.zeros_HSk('H');

	if(GlobalV::VL_IN_H)
	{
		//-------------------------
		// set the local potential
		// in plane wave basis.
		//-------------------------
//		Build_ST_pw bsp;
//		bsp.set_local(ik);	
//		GlobalC::LM.print_HSk('H','C',1.0e-5);

		//--------------------------
		// set the local potential
		// in LCAO basis.
		//--------------------------
		GlobalC::LM.zeros_HSR('H', GlobalC::LNNR.nnr);

		if(GlobalV::NSPIN!=4) 
		{
			this->GK.folding_vl_k(ik);
		}
		else 
		{
			this->GK.folding_vl_k_nc(ik);
		}

		// Peize Lin add 2016-12-03
		if( 5==GlobalC::xcf.iexch_now && 0==GlobalC::xcf.igcx_now )				// HF
		{
			GlobalC::exx_lcao.add_Hexx(ik,1);
		}
		else if( 6==GlobalC::xcf.iexch_now && 8==GlobalC::xcf.igcx_now )			// PBE0
		{
			GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha);
		}
		else if( 9==GlobalC::xcf.iexch_now && 12==GlobalC::xcf.igcx_now )			// HSE
		{
			GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha);
		}
	}


	//-----------------------------------------
	// folding matrix here: S(k) (SlocR->Sloc2)
	// folding matrix here: T(k)+Vnl(k)
	// (Hloc_fixed->Hloc_fixed2)
	//-----------------------------------------
	GlobalC::LM.zeros_HSk('S');
	GlobalC::LM.zeros_HSk('T');
//	std::cout << " after folding Hfixed k." << std::endl;
	GlobalC::LNNR.folding_fixedH(ik);

	//------------------------------------------
	// Add T(k)+Vnl(k)+Vlocal(k)
	// (Hloc2 += Hloc_fixed2), (std::complex matrix)
	//------------------------------------------
//	std::cout << " Folding matrix here." << std::endl;
	GlobalC::LM.update_Hloc2();

/*
	if(GlobalV::NURSE)
	{
		GlobalC::LM.print_HSk('H','R',1.0e-5);
//		GlobalC::LM.print_HSk('S','R',1.0e-5);
	}
	*/
	
	ModuleBase::timer::tick("LCAO_Hamilt","calculate_Hk");
	return;
}

// only need to do the first time.
// available for all k points.
void LCAO_Hamilt::calculate_STNR_k(void)
{
    ModuleBase::TITLE("Hamilt_Linear","calculate_STBR_k");

	//--------------------------------------------
	// set S(R) to zero.
	// the total value of S(R) in this processor
	// is GlobalC::LNNR.nnr.
	// and store in GlobalC::LM.SlocR.
	//--------------------------------------------
	GlobalC::LM.zeros_HSR('S', GlobalC::LNNR.nnr);
    this->genH.calculate_S_no();	

	//------------------------------
	// set T(R) and Vnl(R) to zero.
	// and then calculate it
	// and store in GlobalC::LM.Hloc_fixedR.
	//------------------------------
	GlobalC::LM.zeros_HSR('T', GlobalC::LNNR.nnr);
	


	if(GlobalV::T_IN_H)
	{
		this->genH.calculate_T_no();	
	}


	if(GlobalV::VNL_IN_H)
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
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		std::cout << " ik=" << ik << " ------------------------------------------" << std::endl;
	
		//----------------------------------------------
		// Check the matrix in plane wave basis.
		//----------------------------------------------
		GlobalC::LM.zeros_HSk('S');
		GlobalC::LM.zeros_HSk('T');

		bsp.set_ST(ik, 'S');
		bsp.set_ST(ik, 'T');

		std::cout << " --> PW S" << std::endl;
		GlobalC::LM.print_HSk('S','R',1.0e-5);
		std::cout << " --> PW T" << std::endl;
		GlobalC::LM.print_HSk('T','R',1.0e-5);

		std::string fn = "Sloc2pw.dat";
		GlobalC::LM.output_HSk('S', fn);
		
		//------------------------------------------
		// folding the SlocR and Hloc_fixedR matrix
		// into Sloc2 and Hloc_fixed2 matrix.
		//------------------------------------------
		GlobalC::LM.zeros_HSk('S');
		GlobalC::LM.zeros_HSk('T');
		GlobalC::LNNR.folding_fixedH(ik);
		std::cout << " --> LCAO S" << std::endl;
		GlobalC::LM.print_HSk('S','R',1.0e-5);	
		std::cout << " --> LCAO T+Vnl" << std::endl;
		GlobalC::LM.print_HSk('T','R',1.0e-5);	

		std::string fn2 = "Sloc2lcao.dat";
		GlobalC::LM.output_HSk('S',fn2);

		//----------------
		// test gamma Vnl	
		//----------------
//		GlobalV::GAMMA_ONLY_LOCAL = true;
//		GlobalC::LM.allocate_HS_gamma(GlobalC::ParaO.nloc);
//		GlobalC::LM.zeros_HSgamma('H');
//		GlobalC::UHM.genH.calculate_NL_no( nstart );
//		GlobalV::GAMMA_ONLY_LOCAL = false;
//		std::cout << " Correct LCAO Vnl " << std::endl;
//		GlobalC::LM.print_HSgamma('H');		
//		GlobalC::UHM.genH.calculate_NL_no( nstart );
//		GlobalV::GAMMA_ONLY_LOCAL = false;
//		std::cout << " Correct LCAO Vnl " << std::endl;
//		GlobalC::LM.print_HSgamma('H');		
		
	}
	
    return;
}


void LCAO_Hamilt::calculate_STN_R(void)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R");

    //int iat = 0;
    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    GlobalC::LM.allocate_Hloc_fixedR_tr();
    GlobalC::LM.allocate_HR_tr();
    GlobalC::LM.allocate_SlocR_tr();

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GlobalC::GridD.Find_atom(tau1);
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
                        //const int I0 = GlobalC::GridD.getNatom(ad0);
                        //const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        //const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

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

                    ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);
                    R_x = (int) (dR.x - R_minX);
                    R_y = (int) (dR.y - R_minY);
                    R_z = (int) (dR.z - R_minZ);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            int iic;
                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic=mu+nu*GlobalC::ParaO.nrow;
                            }
                            else
                            {
                                iic=mu*GlobalC::ParaO.ncol+nu;
                            }

                            if(GlobalV::NSPIN!=4)
                            {
                                GlobalC::LM.SlocR_tr[R_x][R_y][R_z][iic] = GlobalC::LM.SlocR[index];
                                GlobalC::LM.Hloc_fixedR_tr[R_x][R_y][R_z][iic] = GlobalC::LM.Hloc_fixedR[index];
                            }
                            else
                            {
                                GlobalC::LM.SlocR_tr_soc[R_x][R_y][R_z][iic] = GlobalC::LM.SlocR_soc[index];
                                GlobalC::LM.Hloc_fixedR_tr_soc[R_x][R_y][R_z][iic] = GlobalC::LM.Hloc_fixedR_soc[index];
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
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R_sparse");

    //int iat = 0;
    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    GlobalC::LM.allocate_HS_R_sparse();

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GlobalC::GridD.Find_atom(tau1);
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
                        //const int I0 = GlobalC::GridD.getNatom(ad0);
                        //const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        //const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

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

                    ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);
                    R_x = (int) (dR.x - R_minX);
                    R_y = (int) (dR.y - R_minY);
                    R_z = (int) (dR.z - R_minZ);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
								double temp_value = GlobalC::LM.SlocR[index];
								if (abs(temp_value) > sparse_threshold)
								{
									GlobalC::LM.SR_sparse[R_x][R_y][R_z][iw1_all].insert(std::pair<size_t, double>(iw2_all, temp_value));
								}

								temp_value = GlobalC::LM.Hloc_fixedR[index];
								if (abs(temp_value) > sparse_threshold)
								{
									GlobalC::LM.HR_sparse[R_x][R_y][R_z][iw1_all].insert(std::pair<size_t, double>(iw2_all, temp_value));
								}
                            }
                            else
                            {
								std::complex<double> temp_value = GlobalC::LM.SlocR_soc[index];
								if(abs(temp_value) > sparse_threshold)
								{
									GlobalC::LM.SR_soc_sparse[R_x][R_y][R_z][iw1_all].insert(std::pair<size_t, std::complex<double>>(iw2_all, temp_value));
								}

								temp_value = GlobalC::LM.Hloc_fixedR_soc[index];
								if(abs(temp_value) > sparse_threshold)
								{
									GlobalC::LM.HR_soc_sparse[R_x][R_y][R_z][iw1_all].insert(std::pair<size_t, std::complex<double>>(iw2_all, temp_value));
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
	ModuleBase::TITLE("LCAO_Hamilt","calculate_HSR_sparse");

	calculate_STN_R_sparse(sparse_threshold);

	GK.cal_vlocal_R_sparseMatrix(current_spin, sparse_threshold);

	if (INPUT.dft_plus_u)
	{
		if (GlobalV::NSPIN == 4)
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
	ModuleBase::TITLE("LCAO_Hamilt","calculat_HR_dftu_sparse");
	ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

	int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

	double *HR_tmp = new double[GlobalC::ParaO.nloc];
	double *SR_tmp = new double[GlobalC::ParaO.nloc];

	int ir;
	int ic;
	int iic;

	for(int ix=0; ix<R_x; ix++)
    {
        for(int iy=0; iy<R_y; iy++)
        {
            for(int iz=0; iz<R_z; iz++)
            {	
				std::map<size_t, std::map<size_t, double>> &temp_HR_sparse = GlobalC::LM.HR_sparse[ix][iy][iz];
				std::map<size_t, std::map<size_t, double>> &temp_SR_sparse = GlobalC::LM.SR_sparse[ix][iy][iz];

				ModuleBase::GlobalFunc::ZEROS(HR_tmp, GlobalC::ParaO.nloc);
				ModuleBase::GlobalFunc::ZEROS(SR_tmp, GlobalC::ParaO.nloc);

				for (auto &iter : temp_SR_sparse)
				{
					ir = GlobalC::ParaO.trace_loc_row[iter.first];
					for (auto &value : iter.second)
					{
						ic = GlobalC::ParaO.trace_loc_col[value.first];
						if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
						{
							iic = ir + ic * GlobalC::ParaO.nrow;
						}
						else
						{
							iic = ir * GlobalC::ParaO.ncol + ic;
						}
						SR_tmp[iic] = value.second;
					}
				}

				GlobalC::dftu.cal_eff_pot_mat_R_double(current_spin, SR_tmp, HR_tmp);

				for (int i = 0; i < GlobalV::NLOCAL; ++i)
				{
					ir = GlobalC::ParaO.trace_loc_row[i];
					if (ir >= 0)
					{
						for (int j = 0; j < GlobalV::NLOCAL; ++j)
						{
							ic = GlobalC::ParaO.trace_loc_col[j];
							if (ic >= 0)
							{
								if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
								{
									iic = ir + ic * GlobalC::ParaO.nrow;
								}
								else
								{
									iic = ir * GlobalC::ParaO.ncol + ic;
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

	ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

}

void LCAO_Hamilt::calculat_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold)
{
	ModuleBase::TITLE("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");
	ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

	int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

	std::complex<double> *HR_soc_tmp = new std::complex<double>[GlobalC::ParaO.nloc];
	std::complex<double> *SR_soc_tmp = new std::complex<double>[GlobalC::ParaO.nloc];

	int ir;
	int ic;
	int iic;

	for(int ix=0; ix<R_x; ix++)
    {
        for(int iy=0; iy<R_y; iy++)
        {
            for(int iz=0; iz<R_z; iz++)
            {
				std::map<size_t, std::map<size_t, std::complex<double>>> &temp_HR_soc_sparse = GlobalC::LM.HR_soc_sparse[ix][iy][iz];
				std::map<size_t, std::map<size_t, std::complex<double>>> &temp_SR_soc_sparse = GlobalC::LM.SR_soc_sparse[ix][iy][iz];

				ModuleBase::GlobalFunc::ZEROS(HR_soc_tmp, GlobalC::ParaO.nloc);
				ModuleBase::GlobalFunc::ZEROS(SR_soc_tmp, GlobalC::ParaO.nloc);

				for (auto &iter : temp_SR_soc_sparse)
				{
					ir = GlobalC::ParaO.trace_loc_row[iter.first];
					for (auto &value : iter.second)
					{
						ic = GlobalC::ParaO.trace_loc_col[value.first];
						if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
						{
							iic = ir + ic * GlobalC::ParaO.nrow;
						}
						else
						{
							iic = ir * GlobalC::ParaO.ncol + ic;
						}
						SR_soc_tmp[iic] = value.second;
					}
				}

				GlobalC::dftu.cal_eff_pot_mat_R_complex_double(current_spin, SR_soc_tmp, HR_soc_tmp);

				for (int i = 0; i < GlobalV::NLOCAL; ++i)
				{
					ir = GlobalC::ParaO.trace_loc_row[i];
					if (ir >= 0)
					{
						for (int j = 0; j < GlobalV::NLOCAL; ++j)
						{
							ic = GlobalC::ParaO.trace_loc_col[j];
							if (ic >= 0)
							{
								if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
								{
									iic = ir + ic * GlobalC::ParaO.nrow;
								}
								else
								{
									iic = ir * GlobalC::ParaO.ncol + ic;
								}

								if (abs(HR_soc_tmp[iic]) > sparse_threshold)
								{
									std::complex<double> &value = temp_HR_soc_sparse[i][j];
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

	ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

}

void LCAO_Hamilt::destroy_all_HSR_sparse(void)
{
	GlobalC::LM.destroy_HS_R_sparse();
}
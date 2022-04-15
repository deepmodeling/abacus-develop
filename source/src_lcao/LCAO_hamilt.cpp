#include "../src_pw/global.h"
#include "LCAO_hamilt.h"
#include "build_st_pw.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "dftu.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_xc/xc_functional.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"	//caoyu add 2021-07-26
#endif
#include "../module_base/timer.h"

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
        this->GG.prepare(GlobalC::ucell.latvec, GlobalC::ucell.lat0, this->LM);
    
        // calulate the 'S', 'T' and 'Vnl' matrix for gamma algorithms.
        this->calculate_STNR_gamma();	

    }
    else // multiple k-points
    {
        // calculate the 'S', 'T' and 'Vnl' matrix for k-points algorithms.
        this->calculate_STNR_k();

        // calculate the grid integration of 'Vl' matrix for l-points algorithms.
        this->GK.init(GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbzp, GlobalC::pw.nbzp_start, GlobalC::pw.ncxyz, this->LM);

    }

    // initial the overlap matrix is done.	
    this->init_s = true;
    //std::cout << " init_s=" << init_s << std::endl; //delete 2015-09-06, xiaohui
//	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"init_s",init_s);

    ModuleBase::timer::tick("LCAO_Hamilt","set_lcao_matrices");
    return;
}

void LCAO_Hamilt::calculate_Hgamma( const int &ik , vector<ModuleBase::matrix> dm_gamma)				// Peize Lin add ik 2016-12-03
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_Hgamma");
    ModuleBase::timer::tick("LCAO_Hamilt","cal_Hgamma");

    // Set the matrix 'H' to zero.
    this->LM->zeros_HSgamma('H'); // 3 stands for Hloc.

    bool local_pw = false;

    if(local_pw)
    {
        std::cout << "\n Call build_H in plane wave basis!" << std::endl;
         // Use plane wave basis to calculate 'Vl' matrix. 
        Build_ST_pw bsp(this->LM);
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

        #ifdef __MPI //liyuanbo 2022/2/23
            // Peize Lin add 2016-12-03
            if(XC_Functional::get_func_type()==4)
            {
                if( Exx_Global::Hybrid_Type::HF == GlobalC::exx_lcao.info.hybrid_type ) //HF
                {
                    GlobalC::exx_lcao.add_Hexx(ik,1, *this->LM);
                }
                else if( Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type )			// PBE0
                {
                    GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha, *this->LM);
                }
                else if( Exx_Global::Hybrid_Type::HSE  == GlobalC::exx_lcao.info.hybrid_type )			// HSE
                {
                    GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha, *this->LM);
                }
            }
        #endif
        }

        time_t time_vlocal_end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("vlocal integration",time_vlocal_start,time_vlocal_end);
    }
    
#ifdef __DEEPKS	//caoyu add 2021-07-26 for DeePKS

	if (GlobalV::deepks_scf)
    {
        const Parallel_Orbitals* pv = this->LM->ParaV;
        GlobalC::ld.cal_projected_DM(dm_gamma[0],
            GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->trace_loc_row,
            pv->trace_loc_col);
    	GlobalC::ld.cal_descriptor();        
		GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
		GlobalC::ld.add_v_delta(GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->trace_loc_row,
            pv->trace_loc_col,
            pv->nrow,
            pv->ncol);
        for(int iic=0;iic<pv->nloc;iic++)
        {
            this->LM->Hloc[iic] += GlobalC::ld.H_V_delta[iic];
        }
	}
	
#endif
	
	//add T+VNL+Vl matrix.
	this->LM->update_Hloc();

	//test
	if(GlobalV::NURSE)
	{
		this->LM->print_HSgamma('S'); // S
		this->LM->print_HSgamma('T');
		this->LM->print_HSgamma('H');
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

    this->LM->zeros_HSgamma('S');    	

    this->genH.calculate_S_no(this->LM->Sloc.data());	

    //this->LM->print_HSgamma('S');

    //-------------------------------------
    // test using plane wave calculations.
    // all the matrixs are stored in this->LM->
    // this->LM->allocate_HS_k(this->LM->ParaV->nloc);
    // Build_ST_pw bsp;
    // bsp.set_ST(0, 'S');
    // this->LM->print_HSk('S','R',1.0e-5);
    //-------------------------------------

    // set T and Vnl matrix to zero.
    // 2 stands for this->LM->Hloc_fixed matrix.
    this->LM->zeros_HSgamma('T'); 

    //add nonlocal pseudopotential matrix element
    time_t time_vnl_start = time(NULL);
    if(GlobalV::VNL_IN_H)
    {
        genH.calculate_NL_no(this->LM->Hloc_fixed.data());
    }
    time_t time_vnl_end = time(NULL);

//	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Time to calculate <psi|Vnl|psi>", std::difftime(time_vnl_end, time_vnl_start));
    
    //add kinetic energy matrix element
    time_t time_t_start = time(NULL);
    if(GlobalV::T_IN_H)
    {
        genH.calculate_T_no(this->LM->Hloc_fixed.data());
//		this->LM->print_HSgamma('T');
    }
    time_t time_t_end = time(NULL);

    //	GlobalV::ofs_running << " T+Vnl matrix" << std::endl;
    //this->LM->print_HSgamma('T');

    ModuleBase::GlobalFunc::OUT_TIME("kinetical matrix",time_t_start, time_t_end);
    ModuleBase::GlobalFunc::OUT_TIME("vnl matrix",time_vnl_start, time_vnl_end);

    return;
}


// be called in LOOP_elec::cal_bands(). 
void LCAO_Hamilt::calculate_Hk(const int &ik)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_Hk");
    ModuleBase::timer::tick("LCAO_Hamilt","calculate_Hk");

    // whether you want to calculate the local potential
    // or not, you need to set this matrix to 0.
    this->LM->zeros_HSk('H');

    if(GlobalV::VL_IN_H)
    {
        //-------------------------
        // set the local potential
        // in plane wave basis.
        //-------------------------
//		Build_ST_pw bsp;
//		bsp.set_local(ik);	
//		this->LM->print_HSk('H','C',1.0e-5);

        //--------------------------
        // set the local potential
        // in LCAO basis.
        //--------------------------

        if(GlobalV::NSPIN!=4) 
        {
            this->GK.folding_vl_k(ik);
        }
        else 
        {
            this->GK.folding_vl_k_nc(ik);
        }

    #ifdef __MPI //liyuanbo 2022/2/23
        // Peize Lin add 2016-12-03
        if(XC_Functional::get_func_type()==4)
        {
            if( Exx_Global::Hybrid_Type::HF  == GlobalC::exx_lcao.info.hybrid_type )				// HF
            {
                GlobalC::exx_lcao.add_Hexx(ik,1, *this->LM);
            }
            else if( Exx_Global::Hybrid_Type::PBE0  == GlobalC::exx_lcao.info.hybrid_type )			// PBE0
            {
                GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha, *this->LM);
            }
            else if( Exx_Global::Hybrid_Type::HSE  == GlobalC::exx_lcao.info.hybrid_type )			// HSE
            {
                GlobalC::exx_lcao.add_Hexx(ik,GlobalC::exx_global.info.hybrid_alpha, *this->LM);
            }
        }
    #endif
    }


    //-----------------------------------------
    // folding matrix here: S(k) (SlocR->Sloc2)
    // folding matrix here: T(k)+Vnl(k)
    // (Hloc_fixed->Hloc_fixed2)
    //-----------------------------------------
    this->LM->zeros_HSk('S');
    this->LM->zeros_HSk('T');
//	std::cout << " after folding Hfixed k." << std::endl;
    this->LM->folding_fixedH(ik);

    //------------------------------------------
    // Add T(k)+Vnl(k)+Vlocal(k)
    // (Hloc2 += Hloc_fixed2), (std::complex matrix)
    //------------------------------------------
//	std::cout << " Folding matrix here." << std::endl;
	this->LM->update_Hloc2(ik);
/*
    if(GlobalV::NURSE)
    {
        this->LM->print_HSk('H','R',1.0e-5);
//		this->LM->print_HSk('S','R',1.0e-5);
    }
    */
    
    ModuleBase::timer::tick("LCAO_Hamilt","calculate_Hk");
    return;
}

// only need to do the first time.
// available for all k points.
void LCAO_Hamilt::calculate_STNR_k(void)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STNR_k");

    //--------------------------------------------
    // set S(R) to zero.
    // the total value of S(R) in this processor
    // is GlobalC::LNNR.nnr.
    // and store in this->LM->SlocR.
    //--------------------------------------------
    this->LM->zeros_HSR('S');
    this->genH.calculate_S_no(this->LM->SlocR.data());	

    //------------------------------
    // set T(R) and Vnl(R) to zero.
    // and then calculate it
    // and store in this->LM->Hloc_fixedR.
    //------------------------------
    this->LM->zeros_HSR('T');
    


    if(GlobalV::T_IN_H)
    {
        this->genH.calculate_T_no(this->LM->Hloc_fixedR.data());	
    }


    if(GlobalV::VNL_IN_H)
    {
        this->genH.calculate_NL_no(this->LM->Hloc_fixedR.data());
    }


    return;

    //-----------------------------------
    // this part is used for checking
    // the consistent between LCAO
    // and plane wave basis.
    //-----------------------------------
    
    // check in plane wave basis.	
    Build_ST_pw bsp(this->LM);
    for(int ik=0; ik<GlobalC::kv.nks; ik++)
    {
        std::cout << " ik=" << ik << " ------------------------------------------" << std::endl;
    
        //----------------------------------------------
        // Check the matrix in plane wave basis.
        //----------------------------------------------
        this->LM->zeros_HSk('S');
        this->LM->zeros_HSk('T');

        bsp.set_ST(ik, 'S');
        bsp.set_ST(ik, 'T');

        std::cout << " --> PW S" << std::endl;
        this->LM->print_HSk('S','R',1.0e-5);
        std::cout << " --> PW T" << std::endl;
        this->LM->print_HSk('T','R',1.0e-5);

        std::string fn = "Sloc2pw.dat";
        this->LM->output_HSk('S', fn);
        
        //------------------------------------------
        // folding the SlocR and Hloc_fixedR matrix
        // into Sloc2 and Hloc_fixed2 matrix.
        //------------------------------------------
        this->LM->zeros_HSk('S');
        this->LM->zeros_HSk('T');
        this->LM->folding_fixedH(ik);
        std::cout << " --> LCAO S" << std::endl;
        this->LM->print_HSk('S','R',1.0e-5);	
        std::cout << " --> LCAO T+Vnl" << std::endl;
        this->LM->print_HSk('T','R',1.0e-5);	

        std::string fn2 = "Sloc2lcao.dat";
        this->LM->output_HSk('S',fn2);

        //----------------
        // test gamma Vnl	
        //----------------
//		GlobalV::GAMMA_ONLY_LOCAL = true;
//		this->LM->allocate_HS_gamma(this->LM->ParaV->nloc);
//		this->LM->zeros_HSgamma('H');
//		GlobalC::UHM.genH.calculate_NL_no( nstart );
//		GlobalV::GAMMA_ONLY_LOCAL = false;
//		std::cout << " Correct LCAO Vnl " << std::endl;
//		this->LM->print_HSgamma('H');		
//		GlobalC::UHM.genH.calculate_NL_no( nstart );
//		GlobalV::GAMMA_ONLY_LOCAL = false;
//		std::cout << " Correct LCAO Vnl " << std::endl;
//		this->LM->print_HSgamma('H');		
        
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

    this->LM->allocate_Hloc_fixedR_tr();
    this->LM->allocate_HR_tr();
    this->LM->allocate_SlocR_tr();

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
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            int iic;
                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic=mu+nu*this->LM->ParaV->nrow;
                            }
                            else
                            {
                                iic=mu*this->LM->ParaV->ncol+nu;
                            }

                            if(GlobalV::NSPIN!=4)
                            {
                                this->LM->SlocR_tr[R_x][R_y][R_z][iic] = this->LM->SlocR[index];
                                this->LM->Hloc_fixedR_tr[R_x][R_y][R_z][iic] = this->LM->Hloc_fixedR[index];
                            }
                            else
                            {
                                this->LM->SlocR_tr_soc[R_x][R_y][R_z][iic] = this->LM->SlocR_soc[index];
                                this->LM->Hloc_fixedR_tr_soc[R_x][R_y][R_z][iic] = this->LM->Hloc_fixedR_soc[index];
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

void LCAO_Hamilt::set_R_range_sparse()
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
                this->LM->all_R_coor.insert(temp_R);
            }
        }
    }

    return;
}

void LCAO_Hamilt::calculate_STN_R_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R_sparse");

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

                    Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                if (current_spin == 0)
                                {
                                    temp_value_double = this->LM->SlocR[index];
                                    if (std::abs(temp_value_double) > sparse_threshold)
                                    {
                                        this->LM->SR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                                    }
                                }

                                temp_value_double = this->LM->Hloc_fixedR[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->HR_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                temp_value_complex = this->LM->SlocR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->SR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
                                }

                                temp_value_complex = this->LM->Hloc_fixedR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->HR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
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


void LCAO_Hamilt::calculate_STN_R_sparse_for_S(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R_sparse_for_S");

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

                    Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                temp_value_double = this->LM->SlocR[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->SR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                temp_value_complex = this->LM->SlocR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->SR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
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

    set_R_range_sparse();

    calculate_STN_R_sparse(current_spin, sparse_threshold);

    GK.cal_vlocal_R_sparseMatrix(current_spin, sparse_threshold);

    if (INPUT.dft_plus_u)
    {
        if (GlobalV::NSPIN != 4)
        {
            calculat_HR_dftu_sparse(current_spin, sparse_threshold);
        }
        else
        {
            calculat_HR_dftu_soc_sparse(current_spin, sparse_threshold);
        }
    }

#ifdef __MPI
    if (GlobalC::exx_global.info.hybrid_type==Exx_Global::Hybrid_Type::HF
        || GlobalC::exx_global.info.hybrid_type==Exx_Global::Hybrid_Type::PBE0
        || GlobalC::exx_global.info.hybrid_type==Exx_Global::Hybrid_Type::HSE)
    {
        calculate_HR_exx_sparse(current_spin, sparse_threshold);
    }
#endif

    clear_zero_elements(current_spin, sparse_threshold);

}

void LCAO_Hamilt::calculate_SR_sparse(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_SR_sparse");
    set_R_range_sparse();
    calculate_STN_R_sparse_for_S(sparse_threshold);
}

void LCAO_Hamilt::calculat_HR_dftu_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculat_HR_dftu_sparse");
    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

    int total_R_num = this->LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        auto iter = this->LM->SR_sparse.find(R_coor);
        if (iter != this->LM->SR_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);

    double *HR_tmp = new double[this->LM->ParaV->nloc];
    double *SR_tmp = new double[this->LM->ParaV->nloc];

    int ir;
    int ic;
    int iic;
    auto &temp_HR_sparse = this->LM->HR_sparse[current_spin];

    count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        if (nonzero_num[count] != 0)
        {
            ModuleBase::GlobalFunc::ZEROS(HR_tmp, this->LM->ParaV->nloc);
            ModuleBase::GlobalFunc::ZEROS(SR_tmp, this->LM->ParaV->nloc);

            auto iter = this->LM->SR_sparse.find(R_coor);
            if (iter != this->LM->SR_sparse.end())
            {
                for (auto &row_loop : iter->second)
                {
                    ir = this->LM->ParaV->trace_loc_row[row_loop.first];
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = this->LM->ParaV->trace_loc_col[col_loop.first];
                        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic = ir + ic * this->LM->ParaV->nrow;
                        }
                        else
                        {
                            iic = ir * this->LM->ParaV->ncol + ic;
                        }
                        SR_tmp[iic] = col_loop.second;
                    }
                }
            }

            GlobalC::dftu.cal_eff_pot_mat_R_double(current_spin, SR_tmp, HR_tmp);

            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                ir = this->LM->ParaV->trace_loc_row[i];
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = this->LM->ParaV->trace_loc_col[j];
                        if (ic >= 0)
                        {
                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic = ir + ic * this->LM->ParaV->nrow;
                            }
                            else
                            {
                                iic = ir * this->LM->ParaV->ncol + ic;
                            }

                            if (std::abs(HR_tmp[iic]) > sparse_threshold)
                            {
                                double &value = temp_HR_sparse[R_coor][i][j];
                                value += HR_tmp[iic];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    temp_HR_sparse[R_coor][i].erase(j);
                                }
                            }
                        }
                    }
                }
            }

        }
        
        count++;
    }

    delete[] nonzero_num;
    delete[] HR_tmp;
    delete[] SR_tmp;
    nonzero_num = nullptr;
    HR_tmp = nullptr;
    SR_tmp = nullptr;

    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

}

void LCAO_Hamilt::calculat_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");
    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

    int total_R_num = this->LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        auto iter = this->LM->SR_soc_sparse.find(R_coor);
        if (iter != this->LM->SR_soc_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);

    std::complex<double> *HR_soc_tmp = new std::complex<double>[this->LM->ParaV->nloc];
    std::complex<double> *SR_soc_tmp = new std::complex<double>[this->LM->ParaV->nloc];

    int ir;
    int ic;
    int iic;

    count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        if (nonzero_num[count] != 0)
        {
            ModuleBase::GlobalFunc::ZEROS(HR_soc_tmp, this->LM->ParaV->nloc);
            ModuleBase::GlobalFunc::ZEROS(SR_soc_tmp, this->LM->ParaV->nloc);

            auto iter = this->LM->SR_soc_sparse.find(R_coor);
            if (iter != this->LM->SR_soc_sparse.end())
            {
                for (auto &row_loop : iter->second)
                {
                    ir = this->LM->ParaV->trace_loc_row[row_loop.first];
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = this->LM->ParaV->trace_loc_col[col_loop.first];
                        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic = ir + ic * this->LM->ParaV->nrow;
                        }
                        else
                        {
                            iic = ir * this->LM->ParaV->ncol + ic;
                        }
                        SR_soc_tmp[iic] = col_loop.second;
                    }
                }
            }

            GlobalC::dftu.cal_eff_pot_mat_R_complex_double(current_spin, SR_soc_tmp, HR_soc_tmp);

            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                ir = this->LM->ParaV->trace_loc_row[i];
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = this->LM->ParaV->trace_loc_col[j];
                        if (ic >= 0)
                        {
                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic = ir + ic * this->LM->ParaV->nrow;
                            }
                            else
                            {
                                iic = ir * this->LM->ParaV->ncol + ic;
                            }

                            if (std::abs(HR_soc_tmp[iic]) > sparse_threshold)
                            {
                                std::complex<double> &value = this->LM->HR_soc_sparse[R_coor][i][j];
                                value += HR_soc_tmp[iic];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    this->LM->HR_soc_sparse[R_coor][i].erase(j);
                                }
                            }
                        }
                    }
                }
            }

        }
        
        count++;
    }

    delete[] nonzero_num;
    delete[] HR_soc_tmp;
    delete[] SR_soc_tmp;
    nonzero_num = nullptr;
    HR_soc_tmp = nullptr;
    SR_soc_tmp = nullptr;

    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

}

#ifdef __MPI
// Peize Lin add 2021.11.16
void LCAO_Hamilt::calculate_HR_exx_sparse(const int &current_spin, const double &sparse_threshold)
{
	ModuleBase::TITLE("LCAO_Hamilt","calculate_HR_exx_sparse");
	ModuleBase::timer::tick("LCAO_Hamilt","calculate_HR_exx_sparse");	

	const Abfs::Vector3_Order<int> Rs_period(GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
	if(Rs_period.x<=0 || Rs_period.y<=0 || Rs_period.z<=0)
		throw std::invalid_argument("Rs_period = ("+ModuleBase::GlobalFunc::TO_STRING(Rs_period.x)+","+ModuleBase::GlobalFunc::TO_STRING(Rs_period.y)+","+ModuleBase::GlobalFunc::TO_STRING(Rs_period.z)+").\n"
			+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	const std::vector<Abfs::Vector3_Order<int>> Rs = Abfs::get_Born_von_Karmen_boxes( Rs_period );


	const int ik_begin = (GlobalV::NSPIN==2) ? (current_spin*GlobalC::kv.nks/2) : 0;
	const int ik_end = (GlobalV::NSPIN==2) ? ((current_spin+1)*GlobalC::kv.nks/2) : GlobalC::kv.nks;
	for(const Abfs::Vector3_Order<int> &R : Rs)
	{
		ModuleBase::matrix HexxR;
		for(int ik=ik_begin; ik<ik_end; ++ik)
		{
			ModuleBase::matrix HexxR_tmp;
			if(GlobalV::GAMMA_ONLY_LOCAL)
				HexxR_tmp = GlobalC::exx_global.info.hybrid_alpha
					* GlobalC::exx_lcao.Hexx_para.HK_Gamma_m2D[ik];
			else
				HexxR_tmp = GlobalC::exx_global.info.hybrid_alpha
					* (GlobalC::exx_lcao.Hexx_para.HK_K_m2D[ik]
					* std::exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (R*GlobalC::ucell.latvec)) )).real();

			if(HexxR.c)
				HexxR += HexxR_tmp;
			else
				HexxR = std::move(HexxR_tmp);
		}

		for(int iwt1_local=0; iwt1_local<HexxR.nr; ++iwt1_local)
		{
			const int iwt1_global = (GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
				? this->LM->ParaV->MatrixInfo.col_set[iwt1_local]
				: this->LM->ParaV->MatrixInfo.row_set[iwt1_local];
			for(int iwt2_local=0; iwt2_local<HexxR.nc; ++iwt2_local)
			{
				const int iwt2_global = (GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
					? this->LM->ParaV->MatrixInfo.row_set[iwt2_local]
					: this->LM->ParaV->MatrixInfo.col_set[iwt2_local];
				if(std::abs(HexxR(iwt1_local,iwt2_local)) > sparse_threshold)
				{
					if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
					{
						auto &HR_sparse_ptr = this->LM->HR_sparse[current_spin][R][iwt1_global];
						auto &HR_sparse = HR_sparse_ptr[iwt2_global];
						HR_sparse += HexxR(iwt1_local,iwt2_local);
						if(std::abs(HR_sparse) < sparse_threshold)
							HR_sparse_ptr.erase(iwt2_global);
					}
					else
					{
						auto &HR_sparse_ptr = this->LM->HR_soc_sparse[R][iwt1_global];
						auto &HR_sparse = HR_sparse_ptr[iwt2_global];
						HR_sparse += HexxR(iwt1_local,iwt2_local);
						if(std::abs(HR_sparse) < sparse_threshold)
							HR_sparse_ptr.erase(iwt2_global);
					}
				}
			}
		}
	}

    // In the future it should be changed to mpi communication, since some Hexx(R) of R in Rs may be zeros
    this->LM->all_R_coor.insert(Rs.begin(),Rs.end());
    
	ModuleBase::timer::tick("LCAO_Hamilt","calculate_HR_exx_sparse");	
}
#endif

// in case there are elements smaller than the threshold
void LCAO_Hamilt::clear_zero_elements(const int &current_spin, const double &sparse_threshold)
{
    if(GlobalV::NSPIN != 4)
    {
        for (auto &R_loop : this->LM->HR_sparse[current_spin])
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

        for (auto &R_loop : this->LM->SR_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

    }
    else
    {
        for (auto &R_loop : this->LM->HR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

        for (auto &R_loop : this->LM->SR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

    }

}

void LCAO_Hamilt::destroy_all_HSR_sparse(void)
{
	this->LM->destroy_HS_R_sparse();
}

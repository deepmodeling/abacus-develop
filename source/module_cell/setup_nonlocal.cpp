#include "setup_nonlocal.h"
#include "module_base/parallel_common.h"

#ifdef __LCAO
//#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
// mohan add 2013-08-02
// In order to get rid of the read in file .NONLOCAL.

InfoNonlocal::InfoNonlocal()
{
    this->Beta = new Numerical_Nonlocal[1];
	this->nproj = nullptr;
    this->nprojmax = 0;
    this->rcutmax_Beta = 0.0;
}
InfoNonlocal::~InfoNonlocal()
{
    delete[] Beta;
	delete[] nproj;
}

#include "../module_base/complexmatrix.h"
void InfoNonlocal::Set_NonLocal(
    const int &it, 
    Atom* atom, 
    int &n_projectors,
    const int& kmesh,
    const double& dk,
    const double& dr_uniform,
	std::ofstream &log)
{
	ModuleBase::TITLE("InfoNonlocal","Set_NonLocal");

	// set a pointer
	//Atom* atom = &GlobalC::ucell.atoms[it];

	// get the number of non-local projectors
	n_projectors = atom->ncpp.nbeta;

	const int nh = atom->ncpp.nh;//zhengdy-soc

	// set the nonlocal projector objects
	Numerical_Nonlocal_Lm* tmpBeta_lm = new Numerical_Nonlocal_Lm[n_projectors];

	ModuleBase::ComplexMatrix coefficient_D_nc_in(nh*2, nh*2);//zhengdy-soc

		int lmaxkb = - 1;
		for (int ibeta = 0; ibeta < atom->ncpp.nbeta; ibeta++)
		{
			lmaxkb = std::max( lmaxkb, atom->ncpp.lll[ibeta]);
		}
		Soc soc;
		if(atom->ncpp.has_so)
		{
			soc.rot_ylm(lmaxkb);
			soc.fcoef.create(1, atom->ncpp.nh, atom->ncpp.nh);
		}

		int ip1=0;
		for(int p1 = 0; p1<n_projectors; p1++)//nbeta
		{
			const int lnow = atom->ncpp.lll[p1];
			
			const int l1 = atom->ncpp.lll[p1];
			const double j1 = atom->ncpp.jjj[p1];
			for(int m1=0; m1<2*l1+1; m1++)
			{
				int ip2=0;
				for(int p2 = 0; p2<n_projectors; p2++)
				{
					const int l2 = atom->ncpp.lll[p2];
					const double j2 = atom->ncpp.jjj[p2];
					for(int m2 = 0;m2<2*l2+1;m2++)
					{
						if(l1==l2 && fabs(j1-j2)<1e-7)
						{
							for(int is1=0;is1<2;is1++)
							{
								for(int is2=0;is2<2;is2++)
								{
									if(atom->ncpp.has_so)
									{
										soc.set_fcoef(l1, l2,
											is1, is2,
											m1, m2,
											j1, j2,
											0, ip1, ip2);
										
										coefficient_D_nc_in(ip1 + nh*is1, ip2 + nh*is2) = atom->ncpp.dion(p1,p2) 
										* soc.fcoef(0, is1, is2, ip1, ip2);
										if(p1 != p2) 
										{
											soc.fcoef(0, is1, is2, ip1, ip2) = std::complex<double>(0.0,0.0);
										}
									}
									else
									{
										if(is1==is2 && ip1==ip2)
										{
											coefficient_D_nc_in(ip1 + nh*is1, ip2 + nh*is2) = atom->ncpp.dion(p1,p2);
										}
									}
								}// end is2
							}// end is1
						}// end l1==l2
						ip2++;
					}// end m2
				}// end p2
				assert(ip2==nh);
				ip1++;
			}// end m1

			// only keep the nonzero part.
			int cut_mesh = atom->ncpp.mesh; 
			for(int ir=atom->ncpp.mesh-1; ir>=0; --ir)
			{
				if( std::abs( atom->ncpp.betar(p1,ir) ) > 1.0e-10 )
				{
					cut_mesh = ir; 
					break;
				}
			}
			if(cut_mesh %2 == 0) 
			{
				++cut_mesh;
			}

			double* beta_r = new double[cut_mesh];
			ModuleBase::GlobalFunc::ZEROS(beta_r, cut_mesh);
			for(int ir=0; ir<cut_mesh; ++ir)
			{
				beta_r[ir] = atom->ncpp.betar(p1,ir);
			}

			tmpBeta_lm[p1].set_NL_proj(
					atom->label,
					it, //type
					lnow, // angular momentum L
					cut_mesh, // number of radial mesh
					atom->ncpp.rab,
					atom->ncpp.r, // radial mesh value (a.u.)
					beta_r,
				    kmesh,
					dk,
					dr_uniform); // delta k mesh in reciprocal space

			if(GlobalV::out_element_info)tmpBeta_lm[p1].plot(GlobalV::MY_RANK);

			delete[] beta_r;
		}

		assert(ip1==nh);

		this->Beta[it].set_type_info(
			it, 
			atom->label, 
			atom->ncpp.pp_type, 
			atom->ncpp.lmax, 
			n_projectors, 
			tmpBeta_lm);//zhengdy-soc 2018-09-10

		// mohan add 2021-05-07
		atom->ncpp.set_d_so(coefficient_D_nc_in,n_projectors,nh,1);

	delete[] tmpBeta_lm;

	log << " SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS" << std::endl;
	return;
}


void InfoNonlocal::Read_NonLocal(
	const int &it, 
	Atom* atom,
	int &n_projectors,
	const int &my_rank,
    const int& kmesh,
    const double& dk,
    const double& dr_uniform,
    const std::string& nonlocalFile)
{
	ModuleBase::TITLE("InfoNonlocal","Read_NonLocal");

	std::ifstream ifs;

	// mohan add 2010-09-08.
	// check if the non-local pseudopotential file exist.
	bool open = false;	
	if(my_rank==0)
	{
		ifs.open( nonlocalFile.c_str() );
		if(ifs)
		{
			open = true;
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(open);
#endif
	if(!open)
	{
		std::cout << " Non-local File : " << nonlocalFile << std::endl;
		ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal","Can not find the NONLOCAL file.");
	}
	else
	{
//		GlobalV::ofs_running << " open nonLocal pseudopotential file: " << nonlocal_file[it] << std::endl;
	}


	std::string label;
	std::string ps_type;

	// maximal lmax allowed in this calculation
	int nlmax = 0;

	if(my_rank==0)
	{
		if(ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<HEADER>"))
		{
			ModuleBase::GlobalFunc::READ_VALUE(ifs, label);
			ModuleBase::GlobalFunc::READ_VALUE(ifs, ps_type);
			if(ps_type != "NC")
			{
				ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal","Only available for NC nonlocal pseudopotential");
			}
			ModuleBase::GlobalFunc::READ_VALUE(ifs, nlmax);
//			std::cout << " " << label << " " << ps_type << " " << nlmax << std::endl; 
			assert(nlmax >= -1);
			ModuleBase::GlobalFunc::SCAN_END(ifs,"</HEADER>");
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_string(label);
	Parallel_Common::bcast_string(ps_type);
	Parallel_Common::bcast_int(nlmax);
#endif

	//mohan add 2012-06-09
	if( nlmax != -1 )
	{
		bool find_lmax = false;
		for(int ic=0; ic<atom->ncpp.nbeta; ic++)
		{
			if( nlmax == atom->ncpp.lll[ic] )
			{
				//			std::cout << " nlmax = " << nlmax << std::endl;
				//			std::cout << " lchi = " << atom->lll[ic] << std::endl;
				find_lmax = true;
				break;
			}
		}

		if( !find_lmax )
		{
			std::cout << " For element " << label << std::endl;
			std::cout << " Max L Read in from NONLOCAL = " << nlmax << std::endl;
			for(int ib=0; ib<atom->ncpp.nbeta; ++ib)
			{
				std::cout << " Max L Read in from pseudopotential file = " << atom->ncpp.lll[ib] << std::endl;
			}
			ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal","nlmax != atom->lll");
		}
	}


//	OUT(GlobalV::ofs_running,"Type",it);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"label",label);
//	OUT(GlobalV::ofs_running,"ps_type",ps_type);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nlmax",nlmax);

	//-------------------------------------------
	// if each L has projectors more than once,
	// this needed to be modified.	
	//-------------------------------------------
	int nproj_allowed = nlmax+1;
	ModuleBase::matrix coefficient_D_in(nproj_allowed, nproj_allowed);
	ModuleBase::ComplexMatrix coefficient_D_nc_in(nproj_allowed*2, nproj_allowed*2);

//	OUT(GlobalV::ofs_running,"nproj_allowed",nproj_allowed);

	if(my_rank==0)
	{
		if(ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<DIJ>"))
		{
			//--------------------------------------
			// this parameter is very important!!!
			//--------------------------------------
			ModuleBase::GlobalFunc::READ_VALUE(ifs, n_projectors);
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"n_projectors",n_projectors);
			
			for (int p1 = 0; p1 < n_projectors; p1++)
        	{
            	for (int p2 = 0; p2 < n_projectors; p2++)
            	{
                	int L1_read, L2_read;
                	
					ifs >> L1_read >> L2_read;
					
					assert(L1_read <= nlmax);
					assert(L2_read <= nlmax);
                	
					ifs >> coefficient_D_in(L1_read, L2_read);
					
//					GlobalV::ofs_running << " L1=" << L1_read << " L2=" << L2_read << " Coef=" << coefficient_D_in(L1_read,L2_read) << std::endl;
            	}
        	}
			ModuleBase::GlobalFunc::SCAN_END(ifs,"</DIJ>");
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_int(n_projectors); // mohan add 2010-12-20
//	Parallel_Common::bcast_double(coefficient_D_in.c, coefficient_D_in.nr * coefficient_D_in.nc);
#endif

	Numerical_Nonlocal_Lm* tmpBeta_lm = new Numerical_Nonlocal_Lm[n_projectors];
	int* LfromBeta = new int[n_projectors];
	ModuleBase::GlobalFunc::ZEROS(LfromBeta, n_projectors);

	for(int p1 = 0; p1<n_projectors; p1++)
	{
		int meshr_ps = 0;
		if(my_rank==0)
		{
			if(ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_BETA>", 0))
			{
				int iproj;
				ModuleBase::GlobalFunc::READ_VALUE(ifs, iproj);
				if(iproj!=p1)
				{
					std::cout << " iproj=" << iproj << " p1=" << p1 << std::endl;
					ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal","Check non-local projector index.");
				}
				
				ModuleBase::GlobalFunc::READ_VALUE(ifs, LfromBeta[p1]);
				assert( LfromBeta[p1] >= 0 );
				assert( LfromBeta[p1] <= nlmax );

				ModuleBase::GlobalFunc::READ_VALUE(ifs, meshr_ps);
				if(meshr_ps%2==0)
				{
					std::cout << " meshr_ps = " << meshr_ps << std::endl;
					ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal","meshr_ps must be odd!");
				}
			}
			else
			{
				ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal","<PP_BETA> doesn't match!");
			}
		}// end my_rank==0

//		OUT(GlobalV::ofs_running,"meshr_ps",meshr_ps);

#ifdef __MPI
		Parallel_Common::bcast_int(meshr_ps);
		Parallel_Common::bcast_int(LfromBeta[p1]);
#endif		

		double* radial_ps = new double[meshr_ps];
		double* rab_ps = new double[meshr_ps];
		double* beta_r = new double[meshr_ps];
		ModuleBase::GlobalFunc::ZEROS(radial_ps, meshr_ps);
		ModuleBase::GlobalFunc::ZEROS(rab_ps, meshr_ps);
		ModuleBase::GlobalFunc::ZEROS(beta_r, meshr_ps);

		if(my_rank==0)
		{
            for (int ir = 0; ir < meshr_ps; ir++)
            {
               	ifs >> radial_ps[ir];
               	ifs >> beta_r[ir];
               	ifs >> rab_ps[ir];
            }
		}

#ifdef __MPI
		Parallel_Common::bcast_double(radial_ps, meshr_ps);
		Parallel_Common::bcast_double(beta_r, meshr_ps);
		Parallel_Common::bcast_double(rab_ps, meshr_ps);
#endif
			
//		OUT(GlobalV::ofs_running,"radial_ps max",radial_ps[meshr_ps-1]);

//		std::cout << this->kmesh << std::endl;
        tmpBeta_lm[p1].set_NL_proj(
        		label,
                it, //type
                LfromBeta[p1], //angular momentum L
                meshr_ps, // number of radial mesh
                rab_ps,
                radial_ps,// radial mesh value(a.u.)
                beta_r,
                kmesh,
                dk,
				dr_uniform); // delta k mesh in reciprocal space

		if(GlobalV::out_element_info)tmpBeta_lm[p1].plot(my_rank);

		delete[] radial_ps;
		delete[] rab_ps;
		delete[] beta_r;
		
		if(my_rank==0)
		{
			ModuleBase::GlobalFunc::SCAN_END(ifs,"</PP_BETA>");
		}
	}// end projectors.
	
	this->Beta[it].set_type_info(
		it, 
		label, 
		ps_type, 
		nlmax, 
		n_projectors, 
		tmpBeta_lm);
		
	ifs.close();

	delete[] LfromBeta;
	delete[] tmpBeta_lm;

	return;
}

void InfoNonlocal::setupNonlocal(
    const int& ntype,
    Atom* atoms,
    std::ofstream &log,
    LCAO_Orbitals &orb
)
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//~~~~~~~~~~~~~~~~~~~~~~   2    ~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read in non-local projector for each atom type.
	// In fact this should be improved,
	// the non-local projector should be transferred
	// from .UPF file directly.
	// mohan note 2011-03-04
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		delete[] this->Beta;
		this->Beta = new Numerical_Nonlocal[ntype];

		delete[] this->nproj;
		this->nproj = new int[ntype];
		ModuleBase::GlobalFunc::ZEROS(this->nproj, ntype);
		
		this->nprojmax = 0;
		

		// if true: read in the nonlocal file from file.
		// if false: get nonlocal information from .upf or .vwr directly
		bool readin_nonlocal = false;

		for(int it=0; it<ntype; it++)
		{
			Atom* atom = &atoms[it];
			if(readin_nonlocal)
			{
				this->Read_NonLocal(
					it, 
					atom, 
					this->nproj[it], 
					GlobalV::MY_RANK, 
					orb.get_kmesh(),
					orb.get_dk(),
					orb.get_dr_uniform(),
					orb.orbital_file[it] );	
			}
			else
			{
				this->Set_NonLocal(
					it, 
					atom, 
					this->nproj[it],
					orb.get_kmesh(),
					orb.get_dk(),
					orb.get_dr_uniform(),
					log);
			}
			this->nprojmax = std::max(this->nprojmax, this->nproj[it]);
			//caoyu add 2021-05-24 to reconstruct atom_arrange::set_sr_NL
			this->rcutmax_Beta = std::max(this->rcutmax_Beta, this->Beta[it].get_rcut_max());
		}

		log << " max number of nonlocal projetors among all species is " << this->nprojmax << std::endl; 
	}
    return;
}

#endif
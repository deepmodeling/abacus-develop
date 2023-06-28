#include "unk_overlap_lcao.h"

#include "ctime"
#include "module_base/scalapack_connector.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

unkOverlap_lcao::unkOverlap_lcao()
{
    allocate_flag = false;
    /*
    const int kpoints_number = kv.nkstot;
    lcao_wfc_global = new std::complex<double>**[kpoints_number];
    for(int ik = 0; ik < kpoints_number; ik++)
    {
        lcao_wfc_global[ik] = new std::complex<double>*[GlobalV::NBANDS];
        for(int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            lcao_wfc_global[ik][ib] = new std::complex<double>[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(lcao_wfc_global[ik][ib], GlobalV::NLOCAL);
        }
    }

    cal_tag = new int*[GlobalV::NLOCAL];
    for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
    {
        cal_tag[iw] = new int[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(cal_tag[iw],GlobalV::NLOCAL);
    }
    */
    // GlobalV::ofs_running << "this is unkOverlap_lcao()" << std::endl;
}

unkOverlap_lcao::~unkOverlap_lcao()
{
    if (allocate_flag)
    {
        for (int ik = 0; ik < this->kpoints_number; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                delete lcao_wfc_global[ik][ib];
            }
            delete lcao_wfc_global[ik];
        }
        delete lcao_wfc_global;

        for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
        {
            delete cal_tag[iw];
        }
        delete cal_tag;
    }

    // GlobalV::ofs_running << "this is ~unkOverlap_lcao()" << std::endl;
}

void unkOverlap_lcao::init(const Grid_Technique& gt, std::complex<double>*** wfc_k_grid, const int nkstot)
{
    // std::cout << "unkOverlap_lcao::init start" << std::endl;

    int Lmax_used, Lmax;

    MOT.allocate(GlobalC::ORB.get_ntype(), // number of atom types
                 GlobalC::ORB.get_lmax(),  // max L used to calculate overlap
                 GlobalC::ORB.get_kmesh(), // kpoints, for integration in k space
                 GlobalC::ORB.get_Rmax(),  // max value of radial table
                 GlobalC::ORB.get_dR(),    // delta R, for making radial table
                 GlobalC::ORB.get_dk());   // delta k, for integration in k space

    int exx_lmax = 0;
#ifdef __EXX
    exx_lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif
    MOT.init_Table_Spherical_Bessel(2, 3, Lmax_used, Lmax, exx_lmax, GlobalC::ORB, GlobalC::ucell.infoNL.Beta);

    ModuleBase::Ylm::set_coefficients();

    MGT.init_Gaunt_CH(Lmax);
    MGT.init_Gaunt(Lmax);

    const int T = 0;                                                    // any selected element type
    orb_r.set_orbital_info(GlobalC::ORB.Phi[T].PhiLN(0, 0).getLabel(),  // atom label
                           T,                                           // atom type
                           1,                                           // angular momentum L
                           1,                                           // number of orbitals of this L , just N
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getNr(),     // number of radial mesh
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getRab(),    // the mesh interval in radial mesh
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getRadial(), // radial mesh value(a.u.)
                           Numerical_Orbital_Lm::Psi_Type::Psi,
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getRadial(), // radial wave function
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getNk(),
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getDk(),
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getDruniform(),
                           false,
                           true,
                           GlobalV::CAL_FORCE);

    // array initialization
    allocate_flag = true;
    this->kpoints_number = nkstot;
    if (allocate_flag)
    {
        lcao_wfc_global = new std::complex<double>**[kpoints_number];
		for(int ik = 0; ik < kpoints_number; ik++)
		{
			lcao_wfc_global[ik] = new std::complex<double>*[GlobalV::NBANDS];
			for(int ib = 0; ib < GlobalV::NBANDS; ib++)
			{
				lcao_wfc_global[ik][ib] = new std::complex<double>[GlobalV::NLOCAL];
				ModuleBase::GlobalFunc::ZEROS(lcao_wfc_global[ik][ib], GlobalV::NLOCAL);
			}
		}
	
		cal_tag = new int*[GlobalV::NLOCAL];
		for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
		{
			cal_tag[iw] = new int[GlobalV::NLOCAL];
			ModuleBase::GlobalFunc::ZEROS(cal_tag[iw],GlobalV::NLOCAL);
		}
    }

    // translate: get the atomic orbital coefficients of each cpu core
    for (int ik = 0; ik < kpoints_number; ik++)
    {
        get_lcao_wfc_global_ik(gt, lcao_wfc_global[ik], wfc_k_grid[ik]);
    }

#ifdef __MPI
    // parallel scheme
    int nproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	const int total_term = GlobalV::NLOCAL * GlobalV::NLOCAL;
	const int remain = total_term % nproc;
	int local_term = total_term / nproc;
	if(myrank < remain)
	{
		local_term++;
	}
	int start;
	for(int rank = 0; rank < nproc; rank++)
	{
		if(rank == myrank)
		{
			if(myrank < remain) 
			{
				start = myrank * local_term;
			}
			else
			{
				start = myrank * local_term + remain;
			}
		}
	}
#else 
	int start=0;
	int local_term=GlobalV::NLOCAL * GlobalV::NLOCAL;
#endif
	int count = -1;
	for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
		{
			count++;
			if(count >= start && count < (start + local_term))
			{
				cal_tag[iw1][iw2] = 1;
			}			
		}
	}
	
	for(int TA = 0; TA < GlobalC::ucell.ntype; TA++)
	{
		for (int TB = 0;  TB < GlobalC::ucell.ntype; TB++)
		{
			for (int LA=0; LA <= GlobalC::ORB.Phi[TA].getLmax() ; LA++)
			{
				for (int NA = 0; NA < GlobalC::ORB.Phi[TA].getNchi(LA); ++NA)
				{
					for (int LB = 0; LB <= GlobalC::ORB.Phi[TB].getLmax(); ++LB)
					{
						for (int NB = 0; NB < GlobalC::ORB.Phi[TB].getNchi(LB); ++NB)
						{
							center2_orb11[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb11(
									GlobalC::ORB.Phi[TA].PhiLN(LA,NA),								
									GlobalC::ORB.Phi[TB].PhiLN(LB,NB),
									MOT, MGT)));
						}
					}
				}
			}
		}
	}
	
	for(int TA = 0; TA < GlobalC::ucell.ntype; TA++)
	{
		for (int TB = 0;  TB < GlobalC::ucell.ntype; TB++)
		{
			for (int LA=0; LA <= GlobalC::ORB.Phi[TA].getLmax() ; LA++)
			{
				for (int NA = 0; NA < GlobalC::ORB.Phi[TA].getNchi(LA); ++NA)
				{
					for (int LB = 0; LB <= GlobalC::ORB.Phi[TB].getLmax(); ++LB)
					{
						for (int NB = 0; NB < GlobalC::ORB.Phi[TB].getNchi(LB); ++NB)
						{
							center2_orb21_r[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb21(
									GlobalC::ORB.Phi[TA].PhiLN(LA,NA),	
									orb_r,									
									GlobalC::ORB.Phi[TB].PhiLN(LB,NB),
									MOT, MGT)));
						}
					}
				}
			}
		}
	}
	
	for( auto &co1 : center2_orb11 )
		for( auto &co2 : co1.second )
			for( auto &co3 : co2.second )
				for( auto &co4 : co3.second )
					for( auto &co5 : co4.second )
						for( auto &co6 : co5.second )
							co6.second.init_radial_table();
						
	for( auto &co1 : center2_orb21_r )
		for( auto &co2 : co1.second )
			for( auto &co3 : co2.second )
				for( auto &co4 : co3.second )
					for( auto &co5 : co4.second )
						for( auto &co6 : co5.second )
							co6.second.init_radial_table();
	
	//std::cout << "unkOverlap_lcao::init end" << std::endl; 
	return;
}

int unkOverlap_lcao::iw2it(int iw)
{
    int ic, type;
    ic = 0;
    for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i=0; i<(2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           type = it;
                        }
                        ic++;
					}
                }
			}
        }
	}
    return type;
}

int unkOverlap_lcao::iw2ia(int iw)
{
    int ic, na;
    ic = 0;
    for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia<  GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           na = ia;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return na;
}

int unkOverlap_lcao::iw2iL(int iw)
{
	int ic, iL;
    ic = 0;
    for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia<  GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           iL = L;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return iL;
}

int unkOverlap_lcao::iw2iN(int iw)
{
	int ic, iN;
    ic = 0;
    for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia<  GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           iN = N;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return iN;
}

int unkOverlap_lcao::iw2im(int iw)
{
	int ic, im;
    ic = 0;
    for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia<  GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           im = i;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return im;
}

// search for the nearest neighbor atoms
void unkOverlap_lcao::cal_R_number()
{
    // The number of overlaps between atomic orbitals 1 and atomic orbitals 2,
    // or the number of R, is empty when there is no overlap
    orb1_orb2_R.resize(GlobalV::NLOCAL);
    for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
    {
		orb1_orb2_R[iw].resize(GlobalV::NLOCAL);
	}
	
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				Atom* atom2 = &GlobalC::ucell.atoms[T2];
				const double R_direct_x = (double)GlobalC::GridD.getBox(ad).x;
				const double R_direct_y = (double)GlobalC::GridD.getBox(ad).y;
				const double R_direct_z = (double)GlobalC::GridD.getBox(ad).z;
			
				tau2 = GlobalC::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
				if(distance < rcut - 1.0e-15)
				{
                    // translate: the unit of R_car is GlobalC::ucell.lat0
                    ModuleBase::Vector3<double> R_car = R_direct_x * GlobalC::ucell.a1 + R_direct_y * GlobalC::ucell.a2
                                                        + R_direct_z * GlobalC::ucell.a3;

                    for (int iw1 = 0; iw1 < atom1->nw; iw1++)
                    {
						int orb_index_in_NLOCAL_1 = GlobalC::ucell.itiaiw2iwt( T1, I1, iw1 );
						for(int iw2 = 0; iw2 < atom2->nw; iw2++)
						{
							int orb_index_in_NLOCAL_2 = GlobalC::ucell.itiaiw2iwt( T2, I2, iw2 );
							orb1_orb2_R[orb_index_in_NLOCAL_1][orb_index_in_NLOCAL_2].push_back(R_car);
						} // end iw2
						
					} // end iw1
					
				}
				
			} //  end ad
			
		} // end I1
		
	} // end T1

	return;
}

void unkOverlap_lcao::cal_orb_overlap()
{
	//std::cout << "the cal_orb_overlap is start" << std::endl;
	psi_psi.resize(GlobalV::NLOCAL);
	psi_r_psi.resize(GlobalV::NLOCAL);
	for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
	{
		psi_psi[iw].resize(GlobalV::NLOCAL);
		psi_r_psi[iw].resize(GlobalV::NLOCAL);
	}

	ModuleBase::Vector3<double> origin_point(0.0,0.0,0.0); 

	for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
		{
			//if ( !pv.in_this_processor(iw1,iw2) ) continue;

            // iw1 and iw2 never have overlap
            if (orb1_orb2_R[iw1][iw2].empty())
                continue;

            int atomType1 = iw2it(iw1);
            int ia1 = iw2ia(iw1);
            int N1 = iw2iN(iw1);
            int L1 = iw2iL(iw1);
            int m1 = iw2im(iw1);
            int atomType2 = iw2it(iw2);  int ia2 = iw2ia(iw2);  int N2 = iw2iN(iw2);  int L2 = iw2iL(iw2);  int m2 = iw2im(iw2);
			
			for(int iR = 0; iR < orb1_orb2_R[iw1][iw2].size(); iR++)
			{
				ModuleBase::Vector3<double> r_distance = ( GlobalC::ucell.atoms[atomType2].tau[ia2] - GlobalC::ucell.atoms[atomType1].tau[ia1] + orb1_orb2_R[iw1][iw2][iR] ) * GlobalC::ucell.lat0;
				psi_psi[iw1][iw2].push_back(center2_orb11[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, m2 ));
				
				double overlap_x = -1 * sqrt(ModuleBase::FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 1, m2 ); // m = 1
				double overlap_y = -1 * sqrt(ModuleBase::FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 2, m2 ); // m = -1
				double overlap_z =      sqrt(ModuleBase::FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 0, m2 ); // m =0
				ModuleBase::Vector3<double> overlap( overlap_x,overlap_y,overlap_z );
				
				psi_r_psi[iw1][iw2].push_back(overlap);
			}
		}
	}
	
	//std::cout << "the cal_orb_overlap is end" << std::endl;
	return;
}

// dk 's unit is GlobalC::ucell.tpiba
std::complex<double> unkOverlap_lcao::unkdotp_LCAO(const int ik_L,
                                                   const int ik_R,
                                                   const int iband_L,
                                                   const int iband_R,
                                                   const ModuleBase::Vector3<double> dk,
                                                   const K_Vectors& kv)
{	
	//std::cout << "unkdotp_LCAO start" << std::endl;

	std::complex<double> result(0.0,0.0);
	
	for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
		{
			//if ( !pv.in_this_processor(iw1,iw2) ) continue;
			if( !cal_tag[iw1][iw2] ) 
			{
				//GlobalV::ofs_running << "the no calculate iw1 and iw2 is " << iw1 << "," << iw2 << std::endl;
				continue;
			}
			
			//GlobalV::ofs_running << "the calculate iw1 and iw2 is " << iw1 << "," << iw2 << std::endl;

            // iw1 and iw2 never have overlap
            if (orb1_orb2_R[iw1][iw2].empty())
                continue;

            // e^i( ik_R*Rn - dk*tau1 )
            ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[ iw2it(iw1) ].tau[ iw2ia(iw1) ];
			ModuleBase::Vector3<double> tau2 = GlobalC::ucell.atoms[ iw2it(iw2) ].tau[ iw2ia(iw2) ];
			ModuleBase::Vector3<double> dtau = tau2 -tau1;
			for(int iR = 0; iR < orb1_orb2_R[iw1][iw2].size(); iR++)
			{
				//*
                double kRn = (kv.kvec_c[ik_R] * orb1_orb2_R[iw1][iw2][iR] - dk * tau1) * ModuleBase::TWO_PI;
                std::complex<double> kRn_phase(cos(kRn),sin(kRn));
				std::complex<double> orb_overlap( psi_psi[iw1][iw2][iR],(-dk * GlobalC::ucell.tpiba * psi_r_psi[iw1][iw2][iR]) );
				result = result + conj( lcao_wfc_global[ik_L][iband_L][iw1] ) * lcao_wfc_global[ik_R][iband_R][iw2] * kRn_phase * orb_overlap;
                //*/

                /*
                // test by jingan
                // R_tem is the vector of the center of the orbitals of iw1 and iw2
                ModuleBase::Vector3<double> R_tem = dtau + orb1_orb2_R[iw1][iw2][iR];
                double kRn = ( kv.kvec_c[ik_R] * orb1_orb2_R[iw1][iw2][iR] - dk * tau1 - 0.5 * dk * R_tem ) *
                ModuleBase::TWO_PI; std::complex<double>  kRn_phase(cos(kRn),sin(kRn)); double psi_r_psi_overlap = -dk *
                GlobalC::ucell.tpiba * psi_r_psi[iw1][iw2][iR] + 0.5 * dk * R_tem * ModuleBase::TWO_PI *
                psi_psi[iw1][iw2][iR]; std::complex<double> orb_overlap( psi_psi[iw1][iw2][iR], psi_r_psi_overlap );
                result = result + conj( lcao_wfc_global[ik_L][iband_L][iw1] ) * lcao_wfc_global[ik_R][iband_R][iw2] *
                kRn_phase * orb_overlap;
                // test by jingan
                */
            }
        }
    }

#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
	result = std::complex<double>(out_date_real,out_date_imag);
#endif	
	
	return result;
}

void unkOverlap_lcao::get_lcao_wfc_global_ik(const Grid_Technique& gt, std::complex<double> **ctot, std::complex<double> **cc)
{
	std::complex<double>* ctot_send = new std::complex<double>[GlobalV::NBANDS*GlobalV::NLOCAL];

#ifdef __MPI
	MPI_Status status;
#endif

	for (int i=0; i<GlobalV::DSIZE; i++)
	{
		if (GlobalV::DRANK==0)
		{
			if (i==0)
			{
				// get the wave functions from 'ctot',
				// save them in the matrix 'c'.
				for (int iw=0; iw<GlobalV::NLOCAL; iw++)
				{
					const int mu_local = gt.trace_lo[iw];
					if (mu_local >= 0)
					{
						for (int ib=0; ib<GlobalV::NBANDS; ib++)
						{
							//ctot[ib][iw] = cc[ib][mu_local];
							ctot_send[ib*GlobalV::NLOCAL+iw] = cc[ib][mu_local];
						}
					}
				}
			}
			else
			{
			#ifdef __MPI
				int tag;
				// receive lgd2
				int lgd2 = 0;
				tag = i * 3;
				MPI_Recv(&lgd2, 1, MPI_INT, i, tag, DIAG_WORLD, &status);
				if(lgd2==0)
				{

				}
				else
				{
					// receive trace_lo2
					tag = i * 3 + 1;
					int* trace_lo2 = new int[GlobalV::NLOCAL];
					MPI_Recv(trace_lo2, GlobalV::NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);
					// receive crecv
					std::complex<double>* crecv = new std::complex<double>[GlobalV::NBANDS*lgd2];
					ModuleBase::GlobalFunc::ZEROS(crecv, GlobalV::NBANDS*lgd2);
					tag = i * 3 + 2;
					MPI_Recv(crecv,GlobalV::NBANDS*lgd2,MPI_DOUBLE_COMPLEX,i,tag,DIAG_WORLD, &status);
				
					for (int ib=0; ib<GlobalV::NBANDS; ib++)
					{
						for (int iw=0; iw<GlobalV::NLOCAL; iw++)
						{
							const int mu_local = trace_lo2[iw];
							if (mu_local>=0)
							{
								//ctot[ib][iw] = crecv[mu_local*GlobalV::NBANDS+ib];
								ctot_send[ib*GlobalV::NLOCAL+iw] = crecv[mu_local*GlobalV::NBANDS+ib];
							}
						}
					}
				
					delete[] crecv;
					delete[] trace_lo2;
				}
			#endif
			}
		}// end GlobalV::DRANK=0
		else if ( i == GlobalV::DRANK)
		{
		#ifdef __MPI
			int tag;

			// send gt.lgd
			tag = GlobalV::DRANK * 3;
			MPI_Send(&gt.lgd, 1, MPI_INT, 0, tag, DIAG_WORLD);

			if(gt.lgd != 0)
			{
				// send trace_lo
				tag = GlobalV::DRANK * 3 + 1;
				MPI_Send(gt.trace_lo, GlobalV::NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

				// send cc
				std::complex<double>* csend = new std::complex<double>[GlobalV::NBANDS*gt.lgd];
				ModuleBase::GlobalFunc::ZEROS(csend, GlobalV::NBANDS*gt.lgd);

				for (int ib=0; ib<GlobalV::NBANDS; ib++)
				{
					for (int mu=0; mu<gt.lgd; mu++)
					{
						csend[mu*GlobalV::NBANDS+ib] = cc[ib][mu];
					}
				}
			
				tag = GlobalV::DRANK * 3 + 2;
				MPI_Send(csend, GlobalV::NBANDS*gt.lgd, MPI_DOUBLE_COMPLEX, 0, tag, DIAG_WORLD);

				delete[] csend;
			}
		#endif
		}// end i==GlobalV::DRANK
		#ifdef __MPI
		MPI_Barrier(DIAG_WORLD);
		#endif
	}

	#ifdef __MPI
	MPI_Bcast(ctot_send,GlobalV::NBANDS*GlobalV::NLOCAL,MPI_DOUBLE_COMPLEX,0,DIAG_WORLD);
	#endif

	for(int ib = 0; ib < GlobalV::NBANDS; ib++)
	{
		for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
		{
			ctot[ib][iw] = ctot_send[ib*GlobalV::NLOCAL+iw];
		}
	}

	delete[] ctot_send;
	return;
}

void unkOverlap_lcao::prepare_midmatrix_pblas(const int ik_L,
                                              const int ik_R,
                                              const ModuleBase::Vector3<double> dk,
                                              std::complex<double>*& midmatrix,
                                              const Parallel_Orbitals& pv,
                                              const K_Vectors& kv)
{
    // ModuleBase::Vector3<double> dk = kv.kvec_c[ik_R] - kv.kvec_c[ik_L];
    midmatrix = new std::complex<double>[pv.nloc];
    ModuleBase::GlobalFunc::ZEROS(midmatrix, pv.nloc);
    for (int iw_row = 0; iw_row < GlobalV::NLOCAL; iw_row++) // global
	{
		for (int iw_col = 0; iw_col < GlobalV::NLOCAL; iw_col++) // global
		{
			int ir = pv.trace_loc_row[ iw_row ]; // local
			int ic = pv.trace_loc_col[ iw_col ]; // local
			
			if(ir >= 0 && ic >= 0)
			{
				int index = ic*pv.nrow+ir;
				ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[ iw2it(iw_row) ].tau[ iw2ia(iw_row) ];		
				for(int iR = 0; iR < orb1_orb2_R[iw_row][iw_col].size(); iR++)
				{
                    double kRn
                        = (kv.kvec_c[ik_R] * orb1_orb2_R[iw_row][iw_col][iR] - dk * tau1) * ModuleBase::TWO_PI;
                    std::complex<double> kRn_phase(cos(kRn),sin(kRn));
					std::complex<double> orb_overlap( psi_psi[iw_row][iw_col][iR],(-dk * GlobalC::ucell.tpiba * psi_r_psi[iw_row][iw_col][iR]) );
					midmatrix[index] = midmatrix[index] + kRn_phase * orb_overlap;
				}
			}
		}
	}
}

std::complex<double> unkOverlap_lcao::det_berryphase(const int ik_L,
                                                     const int ik_R,
                                                     const ModuleBase::Vector3<double> dk,
                                                     const int occ_bands,
                                                     Local_Orbital_wfc& lowf,
                                                     const psi::Psi<std::complex<double>>* psi_in,
                                                     const K_Vectors& kv)
{
	const std::complex<double> minus = std::complex<double>(-1.0,0.0);
	std::complex<double> det = std::complex<double>(1.0,0.0);
	std::complex<double> *midmatrix = NULL;
	std::complex<double> *C_matrix = new std::complex<double>[lowf.ParaV->nloc];
	std::complex<double> *out_matrix = new std::complex<double>[lowf.ParaV->nloc];
	ModuleBase::GlobalFunc::ZEROS(C_matrix,lowf.ParaV->nloc);
	ModuleBase::GlobalFunc::ZEROS(out_matrix,lowf.ParaV->nloc);

    this->prepare_midmatrix_pblas(ik_L, ik_R, dk, midmatrix, *lowf.ParaV, kv);

    char transa = 'C';
	char transb = 'N';
	int occBands = occ_bands;
	int nlocal = GlobalV::NLOCAL;
	std::complex<double> alpha={1.0, 0.0}, beta={0.0, 0.0};
	int one = 1;
#ifdef __MPI
	pzgemm_(&transa,&transb,&occBands,&nlocal,&nlocal,&alpha,
			&psi_in[0](ik_L, 0, 0), &one, &one, lowf.ParaV->desc,
							  midmatrix,&one,&one,lowf.ParaV->desc,
													   &beta,
							   C_matrix,&one,&one,lowf.ParaV->desc);
							   
	pzgemm_(&transb,&transb,&occBands,&occBands,&nlocal,&alpha,
								 C_matrix,&one,&one,lowf.ParaV->desc,
			&psi_in[0](ik_R, 0, 0), &one, &one, lowf.ParaV->desc,
														 &beta,
							   out_matrix,&one,&one,lowf.ParaV->desc);	

	//int *ipiv = new int[ lowf.ParaV->nrow+lowf.ParaV->desc[4] ];
	int *ipiv = new int[ lowf.ParaV->nrow ];
	int info;
	pzgetrf_(&occBands,&occBands,out_matrix,&one,&one,lowf.ParaV->desc,ipiv,&info);

	for(int i = 0; i < occBands; i++) // global
	{	
		int ir = lowf.ParaV->trace_loc_row[ i ]; // local
		int ic = lowf.ParaV->trace_loc_col[ i ]; // local
		if(ir >= 0 && ic >= 0)
		{
			int index = ic*lowf.ParaV->nrow+ir;
			if(ipiv[ir] != (i+1))
			{
				det = minus * det * out_matrix[index];
			}
			else
			{
				det = det * out_matrix[index];
			}
		}
	}
	delete[] ipiv;
#endif
	delete[] midmatrix;
	delete[] C_matrix;
	delete[] out_matrix;

#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
	std::complex<double> result;
	MPI_Allreduce(&det , &result , 1, MPI_DOUBLE_COMPLEX , MPI_PROD , DIAG_WORLD);
	return result;
#endif
	
	return det;
}

void unkOverlap_lcao::test(const Grid_Technique& gt, std::complex<double>*** wfc_k_grid, const K_Vectors& kv)
{
	this->init(gt, wfc_k_grid, kv.nkstot);
	this->cal_R_number();
	this->cal_orb_overlap();

    for (int ik = 0; ik < this->kpoints_number; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
			{
				GlobalV::ofs_running << "the global lcao wfc : ik = " << ik << "  ib = " << ib << "  iw = " << iw << "  valuse = " << lcao_wfc_global[ik][ib][iw] << std::endl;
			}
    }

    /*
    for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
    {
        for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
        {
            if(!cal_tag[iw1][iw2]) continue;

            GlobalV::ofs_running << "the cal_tag is not 0: " << iw1 << "  " << iw2 << std::endl;
        }
    }
    */
    /*
    const int index_1 = 3;
    const int index_2 = 4;

    if(!orb1_orb2_R[index_1][index_2].empty())
    {
        for(int iR = 0; iR < orb1_orb2_R[index_1][index_2].size(); iR++)
            std::cout << "the R is " << orb1_orb2_R[index_1][index_2][iR].x << "," <<
    orb1_orb2_R[index_1][index_2][iR].y << "," << orb1_orb2_R[index_1][index_2][iR].z << " and overlap is " <<
    psi_psi[index_1][index_2][iR] << std::endl;
    }

    */
    /*
    ModuleBase::Vector3<double> dk = kv.kvec_c[0] - kv.kvec_c[0];
    GlobalV::ofs_running << "(" << 0 << "," << 0 << ") = " << std::abs(this->unkdotp_LCAO(0,0,0,0,dk)) << std::endl;
    */
    /*
    ModuleBase::Vector3<double> dk = kv.kvec_c[0] - kv.kvec_c[0];
    for(int ib = 0; ib < GlobalV::NBANDS; ib++)
        for(int ib2 = 0; ib2 < GlobalV::NBANDS; ib2++)
            GlobalV::ofs_running << "(" << ib2 << "," << ib << ") = " << std::abs(this->unkdotp_LCAO(0,0,ib2,ib,dk)) <<
    std::endl;
    */
    /*
    double result = 0;
    for(int iw = 0; iw < GlobalV::NLOCAL; iw++)
    {
        std::cout << "the wfc 11 is " << GlobalC::LOWF.wfc_k_grid[11][13][iw] << " and the 23 is " <<
    GlobalC::LOWF.wfc_k_grid[23][13][iw] << std::endl;
    }
    */
}

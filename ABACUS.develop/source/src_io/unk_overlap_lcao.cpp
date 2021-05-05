#include "unk_overlap_lcao.h"
#include "../src_lcao/LCAO_nnr.h"
#include "ctime"
#include "src_global/scalapack_connector.h"


unkOverlap_lcao::unkOverlap_lcao()
{
	allocate_flag = false;
	/*
	const int kpoints_number = kv.nkstot;
	lcao_wfc_global = new complex<double>**[kpoints_number];
	for(int ik = 0; ik < kpoints_number; ik++)
	{
		lcao_wfc_global[ik] = new complex<double>*[NBANDS];
		for(int ib = 0; ib < NBANDS; ib++)
		{
			lcao_wfc_global[ik][ib] = new complex<double>[NLOCAL];
			ZEROS(lcao_wfc_global[ik][ib], NLOCAL);
		}
	}
	
	cal_tag = new int*[NLOCAL];
	for(int iw = 0; iw < NLOCAL; iw++)
	{
		cal_tag[iw] = new int[NLOCAL];
		ZEROS(cal_tag[iw],NLOCAL);
	}
	*/
	//ofs_running << "this is unkOverlap_lcao()" << endl;
}

unkOverlap_lcao::~unkOverlap_lcao()
{
	if(allocate_flag)
	{
		for(int ik = 0; ik < kv.nkstot; ik++)
		{
			for(int ib = 0; ib < NBANDS; ib++)
			{
				delete lcao_wfc_global[ik][ib];
			}
			delete lcao_wfc_global[ik];
		}
		delete lcao_wfc_global;
	
		for(int iw = 0; iw < NLOCAL; iw++)
		{
			delete cal_tag[iw];
		}
		delete cal_tag;
	}
	
	//ofs_running << "this is ~unkOverlap_lcao()" << endl;
}


void unkOverlap_lcao::init()
{	
	//cout << "unkOverlap_lcao::init start" << endl;

	int Lmax_used, Lmax;

	MOT.allocate(
		ORB.get_ntype(),// number of atom types
		ORB.get_lmax(),// max L used to calculate overlap
		ORB.get_kmesh(), // kpoints, for integration in k space
		ORB.get_Rmax(),// max value of radial table
		ORB.get_dR(),// delta R, for making radial table
		ORB.get_dk()); // delta k, for integration in k space
		
	MOT.init_Table_Spherical_Bessel (2, 3, Lmax_used, Lmax, Exx_Abfs::Lmax);

	Ylm::set_coefficients ();

	MGT.init_Gaunt_CH( Lmax );
	MGT.init_Gaunt( Lmax );

	const int T = 0;  //任意选择的元素类型
	orb_r.set_orbital_info(
	ORB.Phi[T].PhiLN(0,0).getLabel(),  //atom label
	T,    //atom type
	1,    //angular momentum L
	1,    //number of orbitals of this L , just N
	ORB.Phi[T].PhiLN(0,0).getNr(),  //number of radial mesh
	ORB.Phi[T].PhiLN(0,0).getRab(), //the mesh interval in radial mesh 
	ORB.Phi[T].PhiLN(0,0).getRadial(),  // radial mesh value(a.u.)
	Numerical_Orbital_Lm::Psi_Type::Psi,
	ORB.Phi[T].PhiLN(0,0).getRadial(),  // radial wave function
	ORB.Phi[T].PhiLN(0,0).getNk(),
	ORB.Phi[T].PhiLN(0,0).getDk(),
	ORB.Phi[T].PhiLN(0,0).getDruniform(),
	false,
	true);
	
	// 数组初始化
	allocate_flag = true;
	const int kpoints_number = kv.nkstot;
	if(allocate_flag)
	{
		lcao_wfc_global = new complex<double>**[kpoints_number];
		for(int ik = 0; ik < kpoints_number; ik++)
		{
			lcao_wfc_global[ik] = new complex<double>*[NBANDS];
			for(int ib = 0; ib < NBANDS; ib++)
			{
				lcao_wfc_global[ik][ib] = new complex<double>[NLOCAL];
				ZEROS(lcao_wfc_global[ik][ib], NLOCAL);
			}
		}
	
		cal_tag = new int*[NLOCAL];
		for(int iw = 0; iw < NLOCAL; iw++)
		{
			cal_tag[iw] = new int[NLOCAL];
			ZEROS(cal_tag[iw],NLOCAL);
		}
	}
	
	
	//获取每个cpu核的原子轨道系数
	for(int ik = 0; ik < kpoints_number; ik++)
	{
		get_lcao_wfc_global_ik(lcao_wfc_global[ik],LOWF.WFC_K[ik]);
	}
	
	// 并行方案
	int nproc,myrank;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	const int total_term = NLOCAL * NLOCAL;
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
	int count = -1;
	for(int iw1 = 0; iw1 < NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < NLOCAL; iw2++)
		{
			count++;
			if(count >= start && count < (start + local_term))
			{
				cal_tag[iw1][iw2] = 1;
			}			
		}
	}
	
	
	
	for(int TA = 0; TA < ucell.ntype; TA++)
	{
		for (int TB = 0;  TB < ucell.ntype; TB++)
		{
			for (int LA=0; LA <= ORB.Phi[TA].getLmax() ; LA++)
			{
				for (int NA = 0; NA < ORB.Phi[TA].getNchi(LA); ++NA)
				{
					for (int LB = 0; LB <= ORB.Phi[TB].getLmax(); ++LB)
					{
						for (int NB = 0; NB < ORB.Phi[TB].getNchi(LB); ++NB)
						{
							center2_orb11[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb11(
									ORB.Phi[TA].PhiLN(LA,NA),								
									ORB.Phi[TB].PhiLN(LB,NB),
									MOT, MGT)));
						}
					}
				}
			}
		}
	}
	
	for(int TA = 0; TA < ucell.ntype; TA++)
	{
		for (int TB = 0;  TB < ucell.ntype; TB++)
		{
			for (int LA=0; LA <= ORB.Phi[TA].getLmax() ; LA++)
			{
				for (int NA = 0; NA < ORB.Phi[TA].getNchi(LA); ++NA)
				{
					for (int LB = 0; LB <= ORB.Phi[TB].getLmax(); ++LB)
					{
						for (int NB = 0; NB < ORB.Phi[TB].getNchi(LB); ++NB)
						{
							center2_orb21_r[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb21(
									ORB.Phi[TA].PhiLN(LA,NA),	
									orb_r,									
									ORB.Phi[TB].PhiLN(LB,NB),
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
	

	
	//cout << "unkOverlap_lcao::init end" << endl; 
	return;
	
}

int unkOverlap_lcao::iw2it(int iw)
{
    int ic, type;
    ic = 0;
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
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
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
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
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
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
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
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
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
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


//寻找近邻原子
void unkOverlap_lcao::cal_R_number()
{
	// 原子轨道1和原子轨道2之间存在overlap的数目,或者说是R的数目，为空时说明没有overlap
	orb1_orb2_R.resize(NLOCAL);
	for(int iw = 0; iw < NLOCAL; iw++)
	{
		orb1_orb2_R[iw].resize(NLOCAL);
	}
	
	Vector3<double> tau1, tau2, dtau;
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			GridD.Find_atom(tau1, T1, I1);
			
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];
				const double R_direct_x = (double)GridD.getBox(ad).x;
				const double R_direct_y = (double)GridD.getBox(ad).y;
				const double R_direct_z = (double)GridD.getBox(ad).z;
			
				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				if(distance < rcut - 1.0e-15)
				{
					// R_car 单位是 ucell.lat0
					Vector3<double> R_car = R_direct_x * ucell.a1 + 
											R_direct_y * ucell.a2 +
											R_direct_z * ucell.a3;
					
					for(int iw1 = 0; iw1 < atom1->nw; iw1++)
					{
						int orb_index_in_NLOCAL_1 = ucell.itiaiw2iwt( T1, I1, iw1 );
						for(int iw2 = 0; iw2 < atom2->nw; iw2++)
						{
							int orb_index_in_NLOCAL_2 = ucell.itiaiw2iwt( T2, I2, iw2 );
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
	//cout << "the cal_orb_overlap is start" << endl;
	psi_psi.resize(NLOCAL);
	psi_r_psi.resize(NLOCAL);
	for(int iw = 0; iw < NLOCAL; iw++)
	{
		psi_psi[iw].resize(NLOCAL);
		psi_r_psi[iw].resize(NLOCAL);
	}
	

	
	Vector3<double> origin_point(0.0,0.0,0.0); 
	
	
	for(int iw1 = 0; iw1 < NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < NLOCAL; iw2++)
		{
			//if ( !ParaO.in_this_processor(iw1,iw2) ) continue;
			
			// iw1 和 iw2 永远没有overlap
			if( orb1_orb2_R[iw1][iw2].empty() ) continue;
			
			int atomType1 = iw2it(iw1);  int ia1 = iw2ia(iw1);  int N1 = iw2iN(iw1);  int L1 = iw2iL(iw1);  int m1 = iw2im(iw1); 
			int atomType2 = iw2it(iw2);  int ia2 = iw2ia(iw2);  int N2 = iw2iN(iw2);  int L2 = iw2iL(iw2);  int m2 = iw2im(iw2);
			
			for(int iR = 0; iR < orb1_orb2_R[iw1][iw2].size(); iR++)
			{
				Vector3<double> r_distance = ( ucell.atoms[atomType2].tau[ia2] - ucell.atoms[atomType1].tau[ia1] + orb1_orb2_R[iw1][iw2][iR] ) * ucell.lat0;
				psi_psi[iw1][iw2].push_back(center2_orb11[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, m2 ));
				
				double overlap_x = -1 * sqrt(FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 1, m2 ); // m = 1
				double overlap_y = -1 * sqrt(FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 2, m2 ); // m = -1
				double overlap_z =      sqrt(FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 0, m2 ); // m =0
				Vector3<double> overlap( overlap_x,overlap_y,overlap_z );
				
				psi_r_psi[iw1][iw2].push_back(overlap);
			}

			
		}
	}
	
	//cout << "the cal_orb_overlap is end" << endl;
	
	return;
}

// dk 's unit is ucell.tpiba
complex<double> unkOverlap_lcao::unkdotp_LCAO(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const Vector3<double> dk)
{	
	//cout << "unkdotp_LCAO start" << endl;

	complex<double> result(0.0,0.0);
	
	for(int iw1 = 0; iw1 < NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < NLOCAL; iw2++)
		{
			//if ( !ParaO.in_this_processor(iw1,iw2) ) continue;
			if( !cal_tag[iw1][iw2] ) 
			{
				//ofs_running << "the no calculate iw1 and iw2 is " << iw1 << "," << iw2 << endl;
				continue;
			}
			
			//ofs_running << "the calculate iw1 and iw2 is " << iw1 << "," << iw2 << endl;
			
			// iw1 和 iw2 永远没有overlap
			if( orb1_orb2_R[iw1][iw2].empty() ) continue;
		
			
			// e^i( ik_R*Rn - dk*tau1 )
			Vector3<double> tau1 = ucell.atoms[ iw2it(iw1) ].tau[ iw2ia(iw1) ];
			Vector3<double> tau2 = ucell.atoms[ iw2it(iw2) ].tau[ iw2ia(iw2) ];
			Vector3<double> dtau = tau2 -tau1;
			for(int iR = 0; iR < orb1_orb2_R[iw1][iw2].size(); iR++)
			{
				//*
				double kRn = ( kv.kvec_c[ik_R] * orb1_orb2_R[iw1][iw2][iR] - dk * tau1 ) * TWO_PI;
				complex<double> kRn_phase(cos(kRn),sin(kRn));
				complex<double> orb_overlap( psi_psi[iw1][iw2][iR],(-dk * ucell.tpiba * psi_r_psi[iw1][iw2][iR]) );
				result = result + conj( lcao_wfc_global[ik_L][iband_L][iw1] ) * lcao_wfc_global[ik_R][iband_R][iw2] * kRn_phase * orb_overlap;
				//*/ 
				
				/*
				// test by jingan
				// R_tem 是 iw1 和 iw2 的轨道中心的矢量
				Vector3<double> R_tem = dtau + orb1_orb2_R[iw1][iw2][iR];
				double kRn = ( kv.kvec_c[ik_R] * orb1_orb2_R[iw1][iw2][iR] - dk * tau1 - 0.5 * dk * R_tem ) * TWO_PI;
				complex<double>  kRn_phase(cos(kRn),sin(kRn));
				double psi_r_psi_overlap = -dk * ucell.tpiba * psi_r_psi[iw1][iw2][iR] + 0.5 * dk * R_tem * TWO_PI * psi_psi[iw1][iw2][iR];
				complex<double> orb_overlap( psi_psi[iw1][iw2][iR], psi_r_psi_overlap );
				result = result + conj( lcao_wfc_global[ik_L][iband_L][iw1] ) * lcao_wfc_global[ik_R][iband_R][iw2] * kRn_phase * orb_overlap;
				// test by jingan
				*/
			}
	
		}
	}
	

#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the NPOOL = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
	result = complex<double>(out_date_real,out_date_imag);
#endif	
	
	return result;
}

void unkOverlap_lcao::get_lcao_wfc_global_ik(complex<double> **ctot, complex<double> **cc)
{
	complex<double>* ctot_send = new complex<double>[NBANDS*NLOCAL];

	MPI_Status status;

	for (int i=0; i<DSIZE; i++)
	{
		if (DRANK==0)
		{
			if (i==0)
			{
				// get the wave functions from 'ctot',
				// save them in the matrix 'c'.
				for (int iw=0; iw<NLOCAL; iw++)
				{
					const int mu_local = GridT.trace_lo[iw];
					if (mu_local >= 0)
					{
						for (int ib=0; ib<NBANDS; ib++)
						{
							//ctot[ib][iw] = cc[ib][mu_local];
							ctot_send[ib*NLOCAL+iw] = cc[ib][mu_local];
						}
					}
				}
			}
			else
			{
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
					int* trace_lo2 = new int[NLOCAL];
					MPI_Recv(trace_lo2, NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

					// receive crecv
					complex<double>* crecv = new complex<double>[NBANDS*lgd2];
					ZEROS(crecv, NBANDS*lgd2);
					tag = i * 3 + 2;
					MPI_Recv(crecv,NBANDS*lgd2,mpicomplex,i,tag,DIAG_WORLD, &status);
				
					for (int ib=0; ib<NBANDS; ib++)
					{
						for (int iw=0; iw<NLOCAL; iw++)
						{
							const int mu_local = trace_lo2[iw];
							if (mu_local>=0)
							{
								//ctot[ib][iw] = crecv[mu_local*NBANDS+ib];
								ctot_send[ib*NLOCAL+iw] = crecv[mu_local*NBANDS+ib];
							}
						}
					}
				
					delete[] crecv;
					delete[] trace_lo2;
				}
			}
		}// end DRANK=0
		else if ( i == DRANK)
		{
			int tag;

			// send GridT.lgd
			tag = DRANK * 3;
			MPI_Send(&GridT.lgd, 1, MPI_INT, 0, tag, DIAG_WORLD);

			if(GridT.lgd != 0)
			{
				// send trace_lo
				tag = DRANK * 3 + 1;
				MPI_Send(GridT.trace_lo, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

				// send cc
				complex<double>* csend = new complex<double>[NBANDS*GridT.lgd];
				ZEROS(csend, NBANDS*GridT.lgd);

				for (int ib=0; ib<NBANDS; ib++)
				{
					for (int mu=0; mu<GridT.lgd; mu++)
					{
						csend[mu*NBANDS+ib] = cc[ib][mu];
					}
				}
			
				tag = DRANK * 3 + 2;
				MPI_Send(csend, NBANDS*GridT.lgd, mpicomplex, 0, tag, DIAG_WORLD);

			

				delete[] csend;

			}
		}// end i==DRANK
		MPI_Barrier(DIAG_WORLD);
	}

	MPI_Bcast(ctot_send,NBANDS*NLOCAL,mpicomplex,0,DIAG_WORLD);

	for(int ib = 0; ib < NBANDS; ib++)
	{
		for(int iw = 0; iw < NLOCAL; iw++)
		{
			ctot[ib][iw] = ctot_send[ib*NLOCAL+iw];
		}
	}

	delete[] ctot_send;

	return;
}

void unkOverlap_lcao::prepare_midmatrix_pblas(const int ik_L, const int ik_R, const Vector3<double> dk, complex<double> *&midmatrix)
{
	//Vector3<double> dk = kv.kvec_c[ik_R] - kv.kvec_c[ik_L];
	midmatrix = new complex<double>[ParaO.nloc];
	ZEROS(midmatrix,ParaO.nloc);
	for (int iw_row = 0; iw_row < NLOCAL; iw_row++) // global
	{
		for (int iw_col = 0; iw_col < NLOCAL; iw_col++) // global
		{
			int ir = ParaO.trace_loc_row[ iw_row ]; // local
			int ic = ParaO.trace_loc_col[ iw_col ]; // local
			
			if(ir >= 0 && ic >= 0)
			{
				int index = ic*ParaO.nrow+ir;
				Vector3<double> tau1 = ucell.atoms[ iw2it(iw_row) ].tau[ iw2ia(iw_row) ];		
				for(int iR = 0; iR < orb1_orb2_R[iw_row][iw_col].size(); iR++)
				{
					double kRn = ( kv.kvec_c[ik_R] * orb1_orb2_R[iw_row][iw_col][iR] - dk * tau1 ) * TWO_PI;
					complex<double> kRn_phase(cos(kRn),sin(kRn));
					complex<double> orb_overlap( psi_psi[iw_row][iw_col][iR],(-dk * ucell.tpiba * psi_r_psi[iw_row][iw_col][iR]) );
					midmatrix[index] = midmatrix[index] + kRn_phase * orb_overlap;
				}
			}
		}
	}
	
}

complex<double> unkOverlap_lcao::det_berryphase(const int ik_L, const int ik_R, const Vector3<double> dk, const int occ_bands)
{
	const complex<double> minus = complex<double>(-1.0,0.0);
	complex<double> det = complex<double>(1.0,0.0);
	complex<double> *midmatrix = NULL;
	complex<double> *C_matrix = new complex<double>[ParaO.nloc];
	complex<double> *out_matrix = new complex<double>[ParaO.nloc];
	ZEROS(C_matrix,ParaO.nloc);
	ZEROS(out_matrix,ParaO.nloc);
	
	this->prepare_midmatrix_pblas(ik_L,ik_R,dk,midmatrix);
	
	//LOC.wfc_dm_2d.wfc_k
	char transa = 'C';
	char transb = 'N';
	int occBands = occ_bands;
	int nlocal = NLOCAL;
	double alpha=1.0, beta=0.0;
	int one = 1;
	pzgemm_(&transa,&transb,&occBands,&nlocal,&nlocal,&alpha,
			LOC.wfc_dm_2d.wfc_k[ik_L].c,&one,&one,ParaO.desc,
							  midmatrix,&one,&one,ParaO.desc,
													   &beta,
							   C_matrix,&one,&one,ParaO.desc);
							   
	pzgemm_(&transb,&transb,&occBands,&occBands,&nlocal,&alpha,
								 C_matrix,&one,&one,ParaO.desc,
			  LOC.wfc_dm_2d.wfc_k[ik_R].c,&one,&one,ParaO.desc,
														 &beta,
							   out_matrix,&one,&one,ParaO.desc);
	

	//int *ipiv = new int[ ParaO.nrow+ParaO.desc[4] ];
	int *ipiv = new int[ ParaO.nrow ];
	int info;
	pzgetrf_(&occBands,&occBands,out_matrix,&one,&one,ParaO.desc,ipiv,&info);
	
	for(int i = 0; i < occBands; i++) // global
	{	
		int ir = ParaO.trace_loc_row[ i ]; // local
		int ic = ParaO.trace_loc_col[ i ]; // local
		if(ir >= 0 && ic >= 0)
		{
			int index = ic*ParaO.nrow+ir;
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
	
	delete[] midmatrix;
	delete[] C_matrix;
	delete[] out_matrix;
	delete[] ipiv;
	
#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the NPOOL = 1.
	complex<double> result;
	MPI_Allreduce(&det , &result , 1, MPI_DOUBLE_COMPLEX , MPI_PROD , DIAG_WORLD);
	return result;
#endif
	
	return det;
}

void unkOverlap_lcao::test()
{
	this->init();
	this->cal_R_number();
	this->cal_orb_overlap();

	for(int ik = 0; ik < kv.nkstot; ik++)
	{
		for(int ib = 0; ib < NBANDS; ib++)
			for(int iw = 0; iw < NLOCAL; iw++)
			{
				ofs_running << "the global lcao wfc : ik = " << ik << "  ib = " << ib << "  iw = " << iw << "  valuse = " << lcao_wfc_global[ik][ib][iw] << endl;
			}
	}

	/*
	for(int iw1 = 0; iw1 < NLOCAL; iw1++)
	{
		for(int iw2 = 0; iw2 < NLOCAL; iw2++)
		{
			if(!cal_tag[iw1][iw2]) continue;
			
			ofs_running << "the cal_tag is not 0: " << iw1 << "  " << iw2 << endl;
		}
	}
	*/
	/*
	const int index_1 = 3;
	const int index_2 = 4;
	
	if(!orb1_orb2_R[index_1][index_2].empty())
	{
		for(int iR = 0; iR < orb1_orb2_R[index_1][index_2].size(); iR++)
			cout << "the R is " << orb1_orb2_R[index_1][index_2][iR].x << "," << orb1_orb2_R[index_1][index_2][iR].y << "," << orb1_orb2_R[index_1][index_2][iR].z << " and overlap is " << psi_psi[index_1][index_2][iR] << endl;
	}
	
	*/
	/*
	Vector3<double> dk = kv.kvec_c[0] - kv.kvec_c[0];
	ofs_running << "(" << 0 << "," << 0 << ") = " << abs(this->unkdotp_LCAO(0,0,0,0,dk)) << endl;
	*/
	/*
	Vector3<double> dk = kv.kvec_c[0] - kv.kvec_c[0];
	for(int ib = 0; ib < NBANDS; ib++)
		for(int ib2 = 0; ib2 < NBANDS; ib2++)
			ofs_running << "(" << ib2 << "," << ib << ") = " << abs(this->unkdotp_LCAO(0,0,ib2,ib,dk)) << endl;
	*/	
	/*
	double result = 0;
	for(int iw = 0; iw < NLOCAL; iw++)
	{
		cout << "the wfc 11 is " << LOWF.WFC_K[11][13][iw] << " and the 23 is " << LOWF.WFC_K[23][13][iw] << endl;
	}
	*/
}









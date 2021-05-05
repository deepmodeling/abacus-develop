#include "cal_r_overlap_R.h"

cal_r_overlap_R::cal_r_overlap_R(){}

cal_r_overlap_R::~cal_r_overlap_R()
{
	if(this->allocate_psi_r_psi)
	{
		for(int ix = 0; ix < R_x_num; ix++)
		{			
			for(int iy = 0; iy < R_y_num; iy++)
			{				
				for(int iz = 0; iz < R_z_num; iz++)
				{				
					delete[] psi_r_psi[ix][iy][iz];				
				}
				delete[] psi_r_psi[ix][iy];
			}
			delete[] psi_r_psi[ix];
		}
		delete[] psi_r_psi;
	}
}

void cal_r_overlap_R::init()
{
	TITLE("cal_r_overlap_R","init");

	this->R_x_num = GridD.getCellX();
    this->R_y_num = GridD.getCellY();
    this->R_z_num = GridD.getCellZ();
	this->R_minX = (int)GridD.getD_minX();
	this->R_minY = (int)GridD.getD_minY();
	this->R_minZ = (int)GridD.getD_minZ();
	
	// allocate for psi_r_psi	
	psi_r_psi = new Vector3<double> ***[R_x_num];
	for(int ix = 0; ix < R_x_num; ix++)
	{
		psi_r_psi[ix] = new Vector3<double> **[R_y_num];
		for(int iy = 0; iy < R_y_num; iy++)
		{
			psi_r_psi[ix][iy] = new Vector3<double> *[R_z_num];
			for(int iz = 0; iz < R_z_num; iz++)
			{				
				psi_r_psi[ix][iy][iz] = new Vector3<double> [ParaO.nloc];				
			}
		}
	}
	
	this->allocate_psi_r_psi = true;
	
	int Lmax_used=0;
	int Lmax=0;
	
	MOT.allocate(
		ORB.get_ntype(),// number of atom types
		ORB.get_lmax(),// max L used to calculate overlap
		ORB.get_kmesh(), // kpoints, for integration in k space
		ORB.get_Rmax(),// max value of radial table
		ORB.get_dR(),// delta R, for making radial table
		ORB.get_dk()); // delta k, for integration in k space
		
	MOT.init_Table_Spherical_Bessel (2, 3, Lmax_used, Lmax);

	Ylm::set_coefficients();

	MGT.init_Gaunt_CH( Lmax );
	MGT.init_Gaunt( Lmax );	

	int T = 0;  // atom type
	int mat_Nr = ORB.Phi[0].PhiLN(0,0).getNr();
	for(int it = 0; it < ucell.ntype; it++)
	{
		int count_Nr = ORB.Phi[it].PhiLN(0,0).getNr();
		if(count_Nr > mat_Nr) 
		{
			mat_Nr = count_Nr;
			T = it;
		}
	}

//	int new_kmesh = ORB.Phi[T].PhiLN(0,0).getNk() * 4;
//  if(new_kmesh%2 == 0) new_kmesh++;

	int new_kmesh = ORB.Phi[T].PhiLN(0,0).getNk();
//	cout << "new_kmesh = " << new_kmesh << endl;


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
	new_kmesh,
	ORB.Phi[T].PhiLN(0,0).getDk(),
	ORB.Phi[T].PhiLN(0,0).getDruniform(),
	false,
	true);

/*
	orbital_phi.resize(ucell.ntype);
	for(int it = 0; it < ucell.ntype; it++)
	{
		orbital_phi[it].resize(ORB.Phi[it].getLmax()+1);
		for(int iL = 0; iL <= ORB.Phi[it].getLmax(); iL++)
		{	
			orbital_phi[it][iL].resize(ORB.Phi[it].getNchi(iL));
			for(int iN = 0; iN < ORB.Phi[it].getNchi(iL); iN++)
			{
				orbital_phi[it][iL][iN].set_orbital_info(
					ORB.Phi[it].PhiLN(iL,iN).getLabel(),  //atom label
					ORB.Phi[it].PhiLN(iL,iN).getType(),    //atom type
					ORB.Phi[it].PhiLN(iL,iN).getL(),    //angular momentum L
					ORB.Phi[it].PhiLN(iL,iN).getChi(),    //number of orbitals of this L , just N
					ORB.Phi[it].PhiLN(iL,iN).getNr(),  //number of radial mesh
					ORB.Phi[it].PhiLN(iL,iN).getRab(), //the mesh interval in radial mesh 
					ORB.Phi[it].PhiLN(iL,iN).getRadial(),  // radial mesh value(a.u.)
					Numerical_Orbital_Lm::Psi_Type::Psi,
					ORB.Phi[it].PhiLN(iL,iN).getPsi(),  // radial wave function
					new_kmesh,
					ORB.Phi[it].PhiLN(iL,iN).getDk(),
					ORB.Phi[it].PhiLN(iL,iN).getDruniform(),
					false,
					true);
					
				//cout << "getDk:   " << ORB.Phi[it].PhiLN(iL,iN).getDk() << endl;
				//cout << "getDruniform:   " << ORB.Phi[it].PhiLN(iL,iN).getDruniform() << endl; 
				//cout << "getNk:   " << ORB.Phi[it].PhiLN(iL,iN).getNk() << endl; 
			} 
		}
	}

*/

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
							//center2_orb11[TA][TB][LA][NA][LB].insert( 
							//	make_pair(NB, Center2_Orb::Orb11(
							//		orbital_phi[TA][LA][NA],									
							//		orbital_phi[TB][LB][NB],
							//		MOT, MGT)));

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
							//center2_orb21_r[TA][TB][LA][NA][LB].insert( 
							//	make_pair(NB, Center2_Orb::Orb21(
							//		orbital_phi[TA][LA][NA],	
							//		orb_r,									
							//		orbital_phi[TB][LB][NB],
							//		MOT, MGT)));

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
	{
		for( auto &co2 : co1.second )
		{
			for( auto &co3 : co2.second )
			{
				for( auto &co4 : co3.second )
				{
					for( auto &co5 : co4.second )
					{
						for( auto &co6 : co5.second )
						{
							co6.second.init_radial_table();
						}
					}
				}
			}
		}
	}
	
	for( auto &co1 : center2_orb21_r )
	{
		for( auto &co2 : co1.second )
		{
			for( auto &co3 : co2.second )
			{
				for( auto &co4 : co3.second )
				{
					for( auto &co5 : co4.second )
					{
						for( auto &co6 : co5.second )
						{
							co6.second.init_radial_table();
						}
					}
				}
			}
		}
	}
						
	return;
}


void cal_r_overlap_R::out_r_overlap_R(const int nspin)
{	
	TITLE("cal_r_overlap_R","out_r_overlap_R");
	timer::tick("cal_r_overlap_R","out_r_overlap_R");

	Vector3<double> tau1, tau2, dtau;
	Vector3<double> origin_point(0.0,0.0,0.0);

    int R_x;
    int R_y;
    int R_z;

	double factor = sqrt(FOUR_PI/3.0);
	
	for(int ix = 0; ix < R_x_num; ix++)
    {
        int dRx = ix + R_minX;
        for(int iy = 0; iy < R_y_num; iy++)
        {
            int dRy = iy + R_minY;
            for(int iz = 0; iz < R_z_num; iz++)
            {
                int dRz = iz + R_minZ;	

				Vector3<double> R_car = Vector3<double>(dRx,dRy,dRz) * ucell.latvec;
				
				int ir,ic;
				for(int iw1 = 0; iw1 < NLOCAL; iw1++)
				{
					ir = ParaO.trace_loc_row[iw1];	
					if(ir >= 0)
					{
						for(int iw2 = 0; iw2 < NLOCAL; iw2++)
						{							
							ic = ParaO.trace_loc_col[iw2];
							if(ic >= 0)
							{
								int icc = ir + ic * ParaO.nrow;
								
								int orb_index_row = iw1 / NPOL;
								int orb_index_col = iw2 / NPOL;
								
								// soc中非对角项为零，两个对角项相同
								int new_index = iw1 - NPOL*orb_index_row 
									+ (iw2 - NPOL*orb_index_col)*NPOL;
								
								if(new_index == 0 || new_index == 3)
								{
									int it1 = iw2it(orb_index_row);  
									int ia1 = iw2ia(orb_index_row);  
									int N1 = iw2iN(orb_index_row);  
									int L1 = iw2iL(orb_index_row);  
									int m1 = iw2im(orb_index_row); 

									int it2 = iw2it(orb_index_col);  
									int ia2 = iw2ia(orb_index_col);  
									int N2 = iw2iN(orb_index_col);  
									int L2 = iw2iL(orb_index_col);  
									int m2 = iw2im(orb_index_col);

									Vector3<double> r_distance = ( ucell.atoms[it2].tau[ia2] 
									- ucell.atoms[it1].tau[ia1] + R_car ) * ucell.lat0;	

double overlap_o = center2_orb11[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, m2 );
									double overlap_x = -1 * factor * 
center2_orb21_r[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 1, m2 ); // m = 1
									double overlap_y = -1 * factor * 
center2_orb21_r[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 2, m2 ); // m = -1
									double overlap_z =      factor * 
center2_orb21_r[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 0, m2 ); // m = 0	

									psi_r_psi[ix][iy][iz][icc] = Vector3<double>( overlap_x,overlap_y,overlap_z ) 
+ ucell.atoms[it1].tau[ia1] * ucell.lat0 * overlap_o;

								}
								else
								{
									psi_r_psi[ix][iy][iz][icc] = Vector3<double>(0.0,0.0,0.0);
								}
							}
						}
					}
				}
			}
		}
	}
	
	// out r_overlap_R file
	ofstream out_r;
	stringstream ssh;
	ssh << global_out_dir << "data-rR-tr_SPIN" << nspin;
	if(DRANK == 0)
	{
		out_r.open(ssh.str().c_str());
		out_r << "Matrix Dimension of vector r(R): " << NLOCAL <<endl;
	}
	
	for(int ix = 0; ix < R_x_num; ix++)
    {
        int dRx = ix + R_minX;
        for(int iy = 0; iy < R_y_num; iy++)
        {
            int dRy = iy + R_minY;
            for(int iz = 0; iz < R_z_num; iz++)
            {
                int dRz = iz + R_minZ;				
				
				int ir,ic;
				for(int i = 0; i < NLOCAL; i++)
				{
					double *liner_x, *liner_y, *liner_z;
					liner_x = new double[NLOCAL];
					liner_y = new double[NLOCAL];
					liner_z = new double[NLOCAL];
					ZEROS(liner_x,NLOCAL);
					ZEROS(liner_y,NLOCAL);
					ZEROS(liner_z,NLOCAL);
					
					ir = ParaO.trace_loc_row[i];
					
					if(ir >= 0)
					{
						for(int j = 0; j < NLOCAL; j++)
						{
							ic = ParaO.trace_loc_col[j];
							if(ic >= 0)
							{
								int iic = ir + ic * ParaO.nrow;
								liner_x[j] = psi_r_psi[ix][iy][iz][iic].x;
								liner_y[j] = psi_r_psi[ix][iy][iz][iic].y;
								liner_z[j] = psi_r_psi[ix][iy][iz][iic].z;
					
							}
						}
					}
					
					Parallel_Reduce::reduce_double_all(liner_x,NLOCAL);
					Parallel_Reduce::reduce_double_all(liner_y,NLOCAL);
					Parallel_Reduce::reduce_double_all(liner_z,NLOCAL);
					
					if(DRANK == 0)
					{
						for(int j = 0; j < NLOCAL; j++)
						{
							if(i==0 && j==0)
							{
								out_r << dRx << " " << dRy << " " << dRz  
								<< "    //R vector(R2 - R1,unit: lattice vector)" <<endl;
							}
							
							out_r << setw(20) << setprecision(9) << setiosflags(ios::scientific) << liner_x[j] << " "
							      << setw(20) << setprecision(9) << setiosflags(ios::scientific) << liner_y[j] << " "
								  << setw(20) << setprecision(9) << setiosflags(ios::scientific) << liner_z[j] << " "
								  << endl;
						}
						
					}
					
					delete[] liner_x;
					delete[] liner_y;
					delete[] liner_z;
				
				}
				
				
			}
		}
	}
	
	if(DRANK == 0) out_r.close();


	timer::tick("cal_r_overlap_R","out_r_overlap_R");
	
	return;
}

int cal_r_overlap_R::iw2it(int iw)
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

// low efficiency ? -- mohan added 2021-02-14
int cal_r_overlap_R::iw2ia(int iw)
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

int cal_r_overlap_R::iw2iL(int iw)
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

int cal_r_overlap_R::iw2iN(int iw)
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

int cal_r_overlap_R::iw2im(int iw)
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

#include "cal_r_overlap_R.h"
#include "../src_lcao/global_fp.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

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

void cal_r_overlap_R::init(const Parallel_Orbitals &pv)
{
    ModuleBase::TITLE("cal_r_overlap_R", "init");
    this->ParaV = &pv;

	this->R_x_num = GlobalC::GridD.getCellX();
    this->R_y_num = GlobalC::GridD.getCellY();
    this->R_z_num = GlobalC::GridD.getCellZ();
	this->R_minX = (int)GlobalC::GridD.getD_minX();
	this->R_minY = (int)GlobalC::GridD.getD_minY();
	this->R_minZ = (int)GlobalC::GridD.getD_minZ();
	
	// allocate for psi_r_psi	
	psi_r_psi = new ModuleBase::Vector3<double> ***[R_x_num];
	for(int ix = 0; ix < R_x_num; ix++)
	{
		psi_r_psi[ix] = new ModuleBase::Vector3<double> **[R_y_num];
		for(int iy = 0; iy < R_y_num; iy++)
		{
			psi_r_psi[ix][iy] = new ModuleBase::Vector3<double> *[R_z_num];
			for(int iz = 0; iz < R_z_num; iz++)
			{				
				psi_r_psi[ix][iy][iz] = new ModuleBase::Vector3<double> [pv.nloc];				
			}
		}
	}
	
	this->allocate_psi_r_psi = true;
	
	int Lmax_used=0;
	int Lmax=0;
	
	MOT.allocate(
		GlobalC::ORB.get_ntype(),// number of atom types
		GlobalC::ORB.get_lmax(),// max L used to calculate overlap
		GlobalC::ORB.get_kmesh(), // kpoints, for integration in k space
		GlobalC::ORB.get_Rmax(),// max value of radial table
		GlobalC::ORB.get_dR(),// delta R, for making radial table
		GlobalC::ORB.get_dk()); // delta k, for integration in k space
		
#ifdef __MPI 
	MOT.init_Table_Spherical_Bessel (2, 3, Lmax_used, Lmax, Exx_Abfs::Lmax,GlobalC::ORB, GlobalC::ucell.infoNL.Beta);
#endif

	ModuleBase::Ylm::set_coefficients();

	MGT.init_Gaunt_CH( Lmax );
	MGT.init_Gaunt( Lmax );	

	int T = 0;  // atom type
	int mat_Nr = GlobalC::ORB.Phi[0].PhiLN(0,0).getNr();
	for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
		int count_Nr = GlobalC::ORB.Phi[it].PhiLN(0,0).getNr();
		if(count_Nr > mat_Nr) 
		{
			mat_Nr = count_Nr;
			T = it;
		}
	}

//	int new_kmesh = GlobalC::ORB.Phi[T].PhiLN(0,0).getNk() * 4;
//  if(new_kmesh%2 == 0) new_kmesh++;

	int new_kmesh = GlobalC::ORB.Phi[T].PhiLN(0,0).getNk();
//	std::cout << "new_kmesh = " << new_kmesh << std::endl;


	orb_r.set_orbital_info(
	GlobalC::ORB.Phi[T].PhiLN(0,0).getLabel(),  //atom label
	T,    //atom type
	1,    //angular momentum L
	1,    //number of orbitals of this L , just N
	GlobalC::ORB.Phi[T].PhiLN(0,0).getNr(),  //number of radial mesh
	GlobalC::ORB.Phi[T].PhiLN(0,0).getRab(), //the mesh interval in radial mesh 
	GlobalC::ORB.Phi[T].PhiLN(0,0).getRadial(),  // radial mesh value(a.u.)
	Numerical_Orbital_Lm::Psi_Type::Psi,
	GlobalC::ORB.Phi[T].PhiLN(0,0).getRadial(),  // radial wave function
	new_kmesh,
	GlobalC::ORB.Phi[T].PhiLN(0,0).getDk(),
	GlobalC::ORB.Phi[T].PhiLN(0,0).getDruniform(),
	false,
	true, GlobalV::CAL_FORCE);

/*
	orbital_phi.resize(GlobalC::ucell.ntype);
	for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
		orbital_phi[it].resize(GlobalC::ORB.Phi[it].getLmax()+1);
		for(int iL = 0; iL <= GlobalC::ORB.Phi[it].getLmax(); iL++)
		{	
			orbital_phi[it][iL].resize(GlobalC::ORB.Phi[it].getNchi(iL));
			for(int iN = 0; iN < GlobalC::ORB.Phi[it].getNchi(iL); iN++)
			{
				orbital_phi[it][iL][iN].set_orbital_info(
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getLabel(),  //atom label
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getType(),    //atom type
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getL(),    //angular momentum L
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getChi(),    //number of orbitals of this L , just N
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getNr(),  //number of radial mesh
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getRab(), //the mesh interval in radial mesh 
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getRadial(),  // radial mesh value(a.u.)
					Numerical_Orbital_Lm::Psi_Type::Psi,
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getPsi(),  // radial wave function
					new_kmesh,
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getDk(),
					GlobalC::ORB.Phi[it].PhiLN(iL,iN).getDruniform(),
					false,
					true);
					
				//std::cout << "getDk:   " << GlobalC::ORB.Phi[it].PhiLN(iL,iN).getDk() << std::endl;
				//std::cout << "getDruniform:   " << GlobalC::ORB.Phi[it].PhiLN(iL,iN).getDruniform() << std::endl; 
				//std::cout << "getNk:   " << GlobalC::ORB.Phi[it].PhiLN(iL,iN).getNk() << std::endl; 
			} 
		}
	}

*/

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
							//center2_orb11[TA][TB][LA][NA][LB].insert( 
							//	make_pair(NB, Center2_Orb::Orb11(
							//		orbital_phi[TA][LA][NA],									
							//		orbital_phi[TB][LB][NB],
							//		MOT, MGT)));

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
							//center2_orb21_r[TA][TB][LA][NA][LB].insert( 
							//	make_pair(NB, Center2_Orb::Orb21(
							//		orbital_phi[TA][LA][NA],	
							//		orb_r,									
							//		orbital_phi[TB][LB][NB],
							//		MOT, MGT)));

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


void cal_r_overlap_R::out_r_overlap_R()
{	
	ModuleBase::TITLE("cal_r_overlap_R","out_r_overlap_R");
	ModuleBase::timer::tick("cal_r_overlap_R","out_r_overlap_R");

	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> origin_point(0.0,0.0,0.0);

    int R_x;
    int R_y;
    int R_z;

	double factor = sqrt(ModuleBase::FOUR_PI/3.0);
	
	for(int ix = 0; ix < R_x_num; ix++)
    {
        int dRx = ix + R_minX;
        for(int iy = 0; iy < R_y_num; iy++)
        {
            int dRy = iy + R_minY;
            for(int iz = 0; iz < R_z_num; iz++)
            {
                int dRz = iz + R_minZ;	

				ModuleBase::Vector3<double> R_car = ModuleBase::Vector3<double>(dRx,dRy,dRz) * GlobalC::ucell.latvec;
				
				int ir,ic;
				for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
				{
					ir = this->ParaV->trace_loc_row[iw1];	
					if(ir >= 0)
					{
						for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
						{							
							ic = this->ParaV->trace_loc_col[iw2];
							if(ic >= 0)
							{
								int icc = ir + ic * this->ParaV->nrow;
								
								int orb_index_row = iw1 / GlobalV::NPOL;
								int orb_index_col = iw2 / GlobalV::NPOL;
								
								// soc中非对角项为零，两个对角项相同
								int new_index = iw1 - GlobalV::NPOL*orb_index_row 
									+ (iw2 - GlobalV::NPOL*orb_index_col)*GlobalV::NPOL;
								
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

									ModuleBase::Vector3<double> r_distance = ( GlobalC::ucell.atoms[it2].tau[ia2] 
									- GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;	

double overlap_o = center2_orb11[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, m2 );
									double overlap_x = -1 * factor * 
center2_orb21_r[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 1, m2 ); // m = 1
									double overlap_y = -1 * factor * 
center2_orb21_r[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 2, m2 ); // m = -1
									double overlap_z =      factor * 
center2_orb21_r[it1][it2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 0, m2 ); // m = 0	

									psi_r_psi[ix][iy][iz][icc] = ModuleBase::Vector3<double>( overlap_x,overlap_y,overlap_z ) 
+ GlobalC::ucell.atoms[it1].tau[ia1] * GlobalC::ucell.lat0 * overlap_o;

								}
								else
								{
									psi_r_psi[ix][iy][iz][icc] = ModuleBase::Vector3<double>(0.0,0.0,0.0);
								}
							}
						}
					}
				}
			}
		}
	}
	
	// out r_overlap_R file
	std::ofstream out_r;
	std::stringstream ssh;
	ssh << GlobalV::global_out_dir << "data-rR-tr.dat";
	if(GlobalV::DRANK == 0)
	{
		out_r.open(ssh.str().c_str());
		out_r << "Matrix Dimension of r(R): " << GlobalV::NLOCAL <<std::endl;
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
				for(int i = 0; i < GlobalV::NLOCAL; i++)
				{
					double *liner_x, *liner_y, *liner_z;
					liner_x = new double[GlobalV::NLOCAL];
					liner_y = new double[GlobalV::NLOCAL];
					liner_z = new double[GlobalV::NLOCAL];
					ModuleBase::GlobalFunc::ZEROS(liner_x,GlobalV::NLOCAL);
					ModuleBase::GlobalFunc::ZEROS(liner_y,GlobalV::NLOCAL);
					ModuleBase::GlobalFunc::ZEROS(liner_z,GlobalV::NLOCAL);
					
					ir = this->ParaV->trace_loc_row[i];
					
					if(ir >= 0)
					{
						for(int j = 0; j < GlobalV::NLOCAL; j++)
						{
							ic = this->ParaV->trace_loc_col[j];
							if(ic >= 0)
							{
								int iic = ir + ic * this->ParaV->nrow;
								liner_x[j] = psi_r_psi[ix][iy][iz][iic].x;
								liner_y[j] = psi_r_psi[ix][iy][iz][iic].y;
								liner_z[j] = psi_r_psi[ix][iy][iz][iic].z;
					
							}
						}
					}
					
					Parallel_Reduce::reduce_double_all(liner_x,GlobalV::NLOCAL);
					Parallel_Reduce::reduce_double_all(liner_y,GlobalV::NLOCAL);
					Parallel_Reduce::reduce_double_all(liner_z,GlobalV::NLOCAL);
					
					if(GlobalV::DRANK == 0)
					{
						for(int j = 0; j < GlobalV::NLOCAL; j++)
						{
							if(i==0 && j==0)
							{
								out_r << dRx << " " << dRy << " " << dRz  
								<< "    //R (lattice vector)" <<std::endl;
							}
							
							out_r << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) << liner_x[j] << " "
							      << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) << liner_y[j] << " "
								  << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) << liner_z[j] << " "
								  << std::endl;
						}
						
					}
					
					delete[] liner_x;
					delete[] liner_y;
					delete[] liner_z;
				
				}
				
				
			}
		}
	}
	
	if(GlobalV::DRANK == 0) out_r.close();


	ModuleBase::timer::tick("cal_r_overlap_R","out_r_overlap_R");
	
	return;
}

int cal_r_overlap_R::iw2it(int iw)
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

// low efficiency ? -- mohan added 2021-02-14
int cal_r_overlap_R::iw2ia(int iw)
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

int cal_r_overlap_R::iw2iL(int iw)
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

int cal_r_overlap_R::iw2iN(int iw)
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

int cal_r_overlap_R::iw2im(int iw)
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

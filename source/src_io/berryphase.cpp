#include "berryphase.h"

bool berryphase::berry_phase_flag=false;

berryphase::berryphase()
{
	GDIR = INPUT.gdir;
}

berryphase::~berryphase()
{
	//ofs_running << "this is ~berryphase()" << endl;
}

void berryphase::get_occupation_bands()
{
	double occupied_bands = static_cast<double>(CHR.nelec/DEGSPIN);	
	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0;
	}
	
	occ_nbands = (int) occupied_bands;
	if(occ_nbands > NBANDS) 
	{
		WARNING_QUIT("berryphase::get_occupation_bands","not enough bands for berryphase, increase band numbers.");
	}
	//ofs_running << "the berryphase's occ_nbands is " << occ_nbands << endl;
}

void berryphase::lcao_init()
{
	#ifdef __LCAO
	TITLE("berryphase","lcao_init");
	lcao_method.init();
	lcao_method.cal_R_number();
	lcao_method.cal_orb_overlap();
	#endif
	return;
}

// this routine 依赖于 kpoint的mp生成方式
void berryphase::set_kpoints(const int direction)
{
	TITLE("berryphase","set_kpoints");

	const int mp_x = kv.nmp[0]; // x,y,z方向k点的数目
	const int mp_y = kv.nmp[1];
	const int mp_z = kv.nmp[2];
	const int num_k = int(kv.nkstot/2);	

	if( direction == 1 ) // 计算x方向
	{
		const int num_string = mp_y * mp_z;
		
		if( NSPIN == 1 || NSPIN == 4 )
		{
			total_string = num_string;
			k_index.resize(total_string);
		}
		else if( NSPIN == 2 )
		{
			total_string = 2 * num_string;
			k_index.resize(total_string);
		}
		
		
		for(int istring = 0; istring < total_string; istring++)
		{
			k_index[istring].resize(mp_x+1); // 加 1 代表着每条string是从k=0到k=G	
		}
		
		int string_index = -1;
		for(int iz = 0; iz < mp_z; iz++)
		{
			for(int iy = 0; iy < mp_y; iy++)
			{
				string_index++;
				for(int ix = 0; ix < mp_x; ix++)
				{
					k_index[string_index][ix] = ix + iy * mp_x + iz * mp_x * mp_y;
					if( ix == (mp_x-1) ) k_index[string_index][ix+1] = k_index[string_index][0];
				}
			}
		}
		
		if( NSPIN == 2 )
		{
			for(int istring = num_string; istring < total_string; istring++)
			{
				for(int count = 0; count < mp_x+1; count++)
				{
					k_index[istring][count] = k_index[istring-num_string][count] + num_k;
				}
			}
		}
		
		nppstr = mp_x+1;
		
	}
	else if( direction == 2 ) // 计算y方向
	{
		const int num_string = mp_x * mp_z;
		
		if( NSPIN == 1 || NSPIN == 4 )
		{
			total_string = num_string;
			k_index.resize(total_string);
		}
		else if( NSPIN == 2 )
		{
			total_string = 2 * num_string;
			k_index.resize(total_string);
		}
	
		for(int istring = 0; istring < total_string; istring++)
		{
			k_index[istring].resize(mp_y+1); // 加 1 代表着每条string是从k=0到k=G	
		}
		
		int string_index = -1;
		for(int iz = 0; iz < mp_z; iz++)
		{
			for(int ix = 0; ix < mp_x; ix++)
			{
				string_index++;
				for(int iy = 0; iy < mp_y; iy++)
				{
					k_index[string_index][iy] = ix + iy * mp_x + iz * mp_x * mp_y;
					if( iy == (mp_y-1) ) k_index[string_index][iy+1] = k_index[string_index][0];
				}
			}
		}
		
		if( NSPIN == 2 )
		{
			for(int istring = num_string; istring < total_string; istring++)
			{
				for(int count = 0; count < mp_y+1; count++)
				{
					k_index[istring][count] = k_index[istring-num_string][count] + num_k;
				}
			}
		}
		
		nppstr = mp_y+1;
		
	}
	else if( direction == 3 ) // 计算z方向
	{
		const int num_string = mp_x * mp_y;
		
		if( NSPIN == 1 || NSPIN == 4 )
		{
			total_string = num_string;
			k_index.resize(total_string);
		}
		else if( NSPIN == 2 )
		{
			total_string = 2 * num_string;
			k_index.resize(total_string);
		}
		
		for(int istring = 0; istring < total_string; istring++)
		{
			k_index[istring].resize(mp_z+1); // 加 1 代表着每条string是从k=0到k=G
		}
		
		int string_index = -1;
		for(int iy = 0; iy < mp_y; iy++)
		{
			for(int ix = 0; ix < mp_x; ix++)
			{
				string_index++;
				for(int iz = 0; iz < mp_z; iz++)
				{
					k_index[string_index][iz] = ix + iy * mp_x + iz * mp_x * mp_y;
					if( iz == (mp_z-1) ) k_index[string_index][iz+1] = k_index[string_index][0];
				}
			}
		}
		
		if( NSPIN == 2 )
		{
			for(int istring = num_string; istring < total_string; istring++)
			{
				for(int count = 0; count < mp_z+1; count++)
				{
					k_index[istring][count] = k_index[istring-num_string][count] + num_k;
				}
			}
		}
		
		nppstr = mp_z+1;
		
	}

	// test by jingan
	/*
	ofs_running << "direction is " << direction << endl;
	ofs_running << "nppstr = " << nppstr << endl;
	ofs_running << "total string is " << total_string << endl;
	for(int istring = 0; istring < total_string; istring++)
	{
		ofs_running << " the string is " << istring << endl;
		for(int count = 0; count < nppstr; count++)
		{
			ofs_running << "(" << kv.kvec_c[ k_index[istring][count] ].x << ","
							   << kv.kvec_c[ k_index[istring][count] ].y << ","
							   << kv.kvec_c[ k_index[istring][count] ].z << ")" << endl;
		}
		
	}
	*/
	// test by jingan
	

}

double berryphase::stringPhase(int index_str, int nbands)
{
	complex<double> zeta(1.0, 0.0);
	ComplexMatrix mat(nbands,nbands);
	int ik_1;
	int ik_2;
	Vector3<double> G(0.0,0.0,0.0);
	Vector3<double> dk = kv.kvec_c[ k_index[index_str][1] ] - kv.kvec_c[ k_index[index_str][0] ];
	//ofs_running << "the string index is " << index_str << endl;
	
	for(int k_start = 0; k_start < (nppstr-1); k_start++)
	{
		ik_1 = k_index[index_str][k_start];
		ik_2 = k_index[index_str][k_start+1];
		
		if(BASIS_TYPE=="pw")
		{
			for (int mb = 0; mb < nbands; mb++)
			{

				for (int nb = 0; nb < nbands; nb++)
				{
					
					
					if(NSPIN!=4)
					{
						if ( k_start == (nppstr-2) )
						{
							if(direction == 1) 
							{
								Vector3<double> tem_G(1.0,0.0,0.0);
								G = tem_G;
							}
							if(direction == 2)
							{
								Vector3<double> tem_G(0.0,1.0,0.0);
								G = tem_G;
							}
							if(direction == 3)
							{
								Vector3<double> tem_G(0.0,0.0,1.0);
								G = tem_G;
							}
							
							mat(nb,mb) = pw_method.unkdotp_G0(ik_1, ik_2, nb, mb, wf.evc, G);
						}
						else 
						{
							mat(nb, mb) = pw_method.unkdotp_G(ik_1, ik_2, nb, mb, wf.evc);
						}
					}
					else
					{
						if ( k_start == (nppstr-2) )
						{
							if(direction == 1) 
							{
								Vector3<double> tem_G(1.0,0.0,0.0);
								G = tem_G;
							}
							if(direction == 2)
							{
								Vector3<double> tem_G(0.0,1.0,0.0);
								G = tem_G;
							}
							if(direction == 3)
							{
								Vector3<double> tem_G(0.0,0.0,1.0);
								G = tem_G;
							}
							
							mat(nb,mb) = pw_method.unkdotp_soc_G0(ik_1, ik_2, nb, mb, wf.evc, G);							
						}
						else  mat(nb, mb) = pw_method.unkdotp_soc_G(ik_1, ik_2, nb, mb, wf.evc);
					}
					
				} // nb
				
			} // mb
			
			complex<double> det(1.0,0.0);
			int info = 0;
			int *ipiv = new int[nbands];
			LapackConnector::zgetrf(nbands, nbands, mat, nbands, ipiv, &info);				
			for (int ib = 0; ib < nbands; ib++)
			{
				if (ipiv[ib] != (ib+1)) det = -det * mat(ib,ib);
				else det = det * mat(ib,ib);
			}
			
			zeta = zeta*det;
			
			delete[] ipiv;
		}
		#ifdef __LCAO
		else if(BASIS_TYPE=="lcao")
		{
			if(NSPIN!=4)
			{
				//complex<double> my_det = lcao_method.det_berryphase(ik_1,ik_2,dk,nbands);
				zeta = zeta * lcao_method.det_berryphase(ik_1,ik_2,dk,nbands);
				// test by jingan
				//ofs_running << "methon 1: det = " << my_det << endl;
				// test by jingan
			}
			else
			{
				
			}
			
			// test by jingan
			/*
			for (int mb = 0; mb < nbands; mb++)
			{

				for (int nb = 0; nb < nbands; nb++)
				{
					
					mat(nb, mb) = lcao_method.unkdotp_LCAO(ik_1,ik_2,nb,mb,dk);
				}
			}
			
			complex<double> det(1.0,0.0);
			int info = 0;
			int *ipiv = new int[nbands];
			LapackConnector::zgetrf(nbands, nbands, mat, nbands, ipiv, &info);				
			for (int ib = 0; ib < nbands; ib++)
			{
				if (ipiv[ib] != (ib+1)) det = -det * mat(ib,ib);
				else det = det * mat(ib,ib);
			}
			
			zeta = zeta*det;
			
			ofs_running << "methon 2: det = " << det << endl;
			
			delete[] ipiv;
			*/
			// test by jingan
		}
		#endif


	}
	
	
	return log(zeta).imag();
}

void berryphase::Berry_Phase(int nbands, double &pdl_elec_tot, int &mod_elec_tot)
{		
	complex<double> cave = 0.0;
	double *phik = new double[total_string];
	double phik_ave = 0.0;
	complex<double> *cphik = new complex<double>[total_string];
	double *wistring = new double[total_string];
	
	
	// electron polarization
	
	// get weight of every string 
	for(int istring = 0; istring < total_string; istring++)
	{
		wistring[istring] = 1.0 / total_string;
		if(NSPIN == 2) wistring[istring] = wistring[istring] * 2;
	}
	
	for(int istring = 0; istring < total_string; istring++)
	{
		phik[istring] = stringPhase(istring,nbands);
		// 将相位转换成复数形式
		cphik[istring] = complex<double>(cos(phik[istring]),sin(phik[istring]));	
		cave = cave + complex<double>(wistring[istring],0.0) * cphik[istring];
		
	}
	
	double theta0 = atan2(cave.imag(),cave.real());
	double dtheta = 0.0;
	for(int istring = 0; istring < total_string; istring++)
	{
		cphik[istring] = cphik[istring] / cave;
		dtheta = atan2(cphik[istring].imag(),cphik[istring].real());
		phik[istring] = (theta0 + dtheta) / (2 * PI);
		phik_ave = phik_ave + wistring[istring] * phik[istring];
		// test by jingan
		//ofs_running << "phik[" << istring << "] = " << phik[istring] << endl;
		// test by jingan
	}
	
	if(NSPIN == 1)
	{
		pdl_elec_tot = 2 * phik_ave;
	}
	else if( NSPIN == 2 || NSPIN == 4 )
	{
		pdl_elec_tot = phik_ave;
	}
	
	if(NSPIN == 1)  // remap to [-1,1]
	{
		pdl_elec_tot = pdl_elec_tot - 2.0 * round(pdl_elec_tot/2.0);
		mod_elec_tot = 2;
	}
	else if( NSPIN == 2 || NSPIN == 4 )  // remap to [-0.5,0.5]
	{
		pdl_elec_tot = pdl_elec_tot - 1.0 * round(pdl_elec_tot/1.0);
		mod_elec_tot = 1;
	}
	
	
	
	delete[] phik;
	delete[] cphik;
	delete[] wistring;
	
	
	//ofs_running << "Berry_Phase end " << endl;

}


void berryphase::Macroscopic_polarization()
{	
	get_occupation_bands();
	
	if( BASIS_TYPE == "lcao" ) this->lcao_init();
	
	
	ofs_running << "\n\n\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " | POLARIZATION CALCULATION:                                          |" << endl;
	ofs_running << " |                  Modern Theory of Polarization                     |" << endl;
	ofs_running << " | calculate the Macroscopic polarization of a crystalline insulator  |" << endl;
	ofs_running << " | by using Berry Phase method.                                       |" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	ofs_running << "\n\n\n\n";
	
	
	// ion polarization	
	double polarization_ion[3]; // 指的是晶格矢量R1，R2，R3方向
	ZEROS(polarization_ion,3);
	// 倒格矢
	Vector3<double> rcell_1(ucell.G.e11,ucell.G.e12,ucell.G.e13);
	Vector3<double> rcell_2(ucell.G.e21,ucell.G.e22,ucell.G.e23);
	Vector3<double> rcell_3(ucell.G.e31,ucell.G.e32,ucell.G.e33);
	int *mod_ion = new int[ucell.nat];
	double *pdl_ion_R1 = new double[ucell.nat];
	double *pdl_ion_R2 = new double[ucell.nat];
	double *pdl_ion_R3 = new double[ucell.nat];  
	ZEROS(mod_ion,ucell.nat);
	ZEROS(pdl_ion_R1,ucell.nat);
	ZEROS(pdl_ion_R2,ucell.nat);
	ZEROS(pdl_ion_R3,ucell.nat);
	
	bool lodd = false;
	int atom_index = 0;
	for(int it = 0; it < ucell.ntype; it++)
	{
		for(int ia = 0; ia < ucell.atoms[it].na; ia++)
		{
			if(ucell.atoms[it].zv % 2 == 1)
			{
				mod_ion[atom_index] = 1;
				lodd = true;
				atom_index++;
			}
			else
			{
				mod_ion[atom_index] = 2;
				atom_index++;
			}
		}
	}
	
	atom_index = 0;
	for(int it = 0; it < ucell.ntype; it++)
	{
		for(int ia = 0; ia < ucell.atoms[it].na; ia++)
		{
			pdl_ion_R1[atom_index] = ucell.atoms[it].zv * (ucell.atoms[it].tau[ia] * rcell_1);
			pdl_ion_R2[atom_index] = ucell.atoms[it].zv * (ucell.atoms[it].tau[ia] * rcell_2);
			pdl_ion_R3[atom_index] = ucell.atoms[it].zv * (ucell.atoms[it].tau[ia] * rcell_3);
			atom_index++;
		}
	}
	

	
	for(int i = 0; i < ucell.nat; i++)
	{		
		if(mod_ion[i] == 1) 
		{
			pdl_ion_R1[i] = pdl_ion_R1[i] - 1.0 * round(pdl_ion_R1[i] / 1.0);
			pdl_ion_R2[i] = pdl_ion_R2[i] - 1.0 * round(pdl_ion_R2[i] / 1.0);
			pdl_ion_R3[i] = pdl_ion_R3[i] - 1.0 * round(pdl_ion_R3[i] / 1.0);
		}
		else if(mod_ion[i] == 2) 
		{
			pdl_ion_R1[i] = pdl_ion_R1[i] - 2.0 * round(pdl_ion_R1[i] / 2.0);
			pdl_ion_R2[i] = pdl_ion_R2[i] - 2.0 * round(pdl_ion_R2[i] / 2.0);
			pdl_ion_R3[i] = pdl_ion_R3[i] - 2.0 * round(pdl_ion_R3[i] / 2.0);
		}
		
		polarization_ion[0] = polarization_ion[0] + pdl_ion_R1[i];
		polarization_ion[1] = polarization_ion[1] + pdl_ion_R2[i];
		polarization_ion[2] = polarization_ion[2] + pdl_ion_R3[i];
		
	}
	if(lodd)
	{
		polarization_ion[0] = polarization_ion[0] - 1.0 * round(polarization_ion[0] / 1.0);
		polarization_ion[1] = polarization_ion[1] - 1.0 * round(polarization_ion[1] / 1.0);
		polarization_ion[2] = polarization_ion[2] - 1.0 * round(polarization_ion[2] / 1.0);
	}
	else
	{
		polarization_ion[0] = polarization_ion[0] - 2.0 * round(polarization_ion[0] / 2.0);
		polarization_ion[1] = polarization_ion[1] - 2.0 * round(polarization_ion[1] / 2.0);
		polarization_ion[2] = polarization_ion[2] - 2.0 * round(polarization_ion[2] / 2.0);
	}
	
	delete[] mod_ion;
	delete[] pdl_ion_R1;
	delete[] pdl_ion_R2;
	delete[] pdl_ion_R3;
	
	// ion polarization	end
	
	// calculate Macroscopic polarization modulus because berry phase
	int modulus;
	if( (!lodd) && (NSPIN==1) ) modulus = 2;
	else modulus = 1;
	
	// test by jingan
	//ofs_running << "ion polarization end" << endl;
	// test by jingan


	// berry phase calculate begin	
	switch(GDIR)
	{
		case 1:
		{
			direction = 1;
			set_kpoints(direction);
			double pdl_elec_tot = 0.0;
			int mod_elec_tot = 0;
			Berry_Phase(occ_nbands, pdl_elec_tot, mod_elec_tot);
		
			const double rmod = ucell.a1.norm() * ucell.lat0;
			const double unit1 = rmod;
			const double unit2 = rmod / ucell.omega;
			const double unit3 = ( rmod / ucell.omega ) * ( 1.60097e-19/pow(5.29177e-11,2) );
			
			ofs_running << " VALUES OF POLARIZATION" << endl;
			ofs_running << endl;
			ofs_running << "  The Ionic Phase: " << setw(10) << fixed << setprecision(5) << polarization_ion[0] << endl;
			ofs_running << " Electronic Phase: " << setw(10) << fixed << setprecision(5) << pdl_elec_tot << endl;
			//ofs_running << " the electronic part polarization is P(ele) = " << rmod * pdl_elec_tot << "   (e/Omega).bohr   in R1 direction" << endl;
			
			// calculate total polarization,add electron part and ions part
			double total_polarization = pdl_elec_tot + polarization_ion[0] ;
			
			Vector3<double> polarization_xyz = ucell.a1;
			polarization_xyz.normalize();
			polarization_xyz = total_polarization * polarization_xyz;

			ofs_running << "\n" << "The calculated polarization direction is in R1 direction" << endl;
			ofs_running << "\n" << " P = " << outFormat(unit1*total_polarization, unit1*modulus, unit1*polarization_xyz) << "(e/Omega).bohr" << endl;
			ofs_running << "\n" << " P = " << outFormat(unit2*total_polarization, unit2*modulus, unit2*polarization_xyz) << "e/bohr^2" << endl;
			ofs_running << "\n" << " P = " << outFormat(unit3*total_polarization, unit3*modulus, unit3*polarization_xyz) << "C/m^2" << endl;			   
			ofs_running << endl;
		
			break;
		
		}
		case 2:
		{
			direction = 2;
			set_kpoints(direction);
			double pdl_elec_tot = 0.0;
			int mod_elec_tot = 0;
			Berry_Phase(occ_nbands, pdl_elec_tot, mod_elec_tot);
		
			const double rmod = ucell.a2.norm() * ucell.lat0;
			const double unit1 = rmod;
			const double unit2 = rmod / ucell.omega;
			const double unit3 = ( rmod / ucell.omega ) * ( 1.60097e-19/pow(5.29177e-11,2) );
			
			ofs_running << " VALUES OF POLARIZATION" << endl;
			ofs_running << endl;
			ofs_running << "  The Ionic Phase: " << setw(10) << fixed << setprecision(5) << polarization_ion[1] << endl;
			ofs_running << " Electronic Phase: " << setw(10) << fixed << setprecision(5) << pdl_elec_tot << endl;
			//ofs_running << " the electronic part polarization is P(ele) = " << rmod * pdl_elec_tot << "   (e/Omega).bohr   in R2 direction" << endl;
		
			// calculate total polarization,add electron part and ions part
			double total_polarization = pdl_elec_tot + polarization_ion[1] ;
			
			Vector3<double> polarization_xyz = ucell.a2;
			polarization_xyz.normalize();
			polarization_xyz = total_polarization * polarization_xyz;
		
			ofs_running << "\n" << "The calculated polarization direction is in R2 direction" << endl;
			ofs_running << "\n"  << " P = " << outFormat(unit1*total_polarization, unit1*modulus, unit1*polarization_xyz) << "(e/Omega).bohr" << endl;
			ofs_running << "\n"  << " P = " << outFormat(unit2*total_polarization, unit2*modulus, unit2*polarization_xyz) << "e/bohr^2" << endl;
			ofs_running << "\n"  << " P = " << outFormat(unit3*total_polarization, unit3*modulus, unit3*polarization_xyz) << "C/m^2" << endl;
			ofs_running << endl;
			
			break;
		
		}
		case 3:
		{
			direction = 3;
			set_kpoints(direction);
			double pdl_elec_tot = 0.0;
			int mod_elec_tot = 0;
			Berry_Phase(occ_nbands, pdl_elec_tot, mod_elec_tot);
		
			const double rmod = ucell.a3.norm() * ucell.lat0;
			const double unit1 = rmod;
			const double unit2 = rmod / ucell.omega;
			const double unit3 = ( rmod / ucell.omega ) * ( 1.60097e-19/pow(5.29177e-11,2) );
			
			ofs_running << " VALUES OF POLARIZATION" << endl;
			ofs_running << endl;
			ofs_running << "  The Ionic Phase: " << setw(10) << fixed << setprecision(5) << polarization_ion[2] << endl;
			ofs_running << " Electronic Phase: " << setw(10) << fixed << setprecision(5) << pdl_elec_tot << endl;
			//ofs_running << " the electronic part polarization is P(ele) = " << rmod * pdl_elec_tot << "   (e/Omega).bohr   in R3 direction" << endl;
		
			// calculate total polarization,add electron part and ions part
			double total_polarization = pdl_elec_tot + polarization_ion[2] ;
			
			Vector3<double> polarization_xyz = ucell.a3;
			polarization_xyz.normalize();
			polarization_xyz = total_polarization * polarization_xyz;

			ofs_running << "\n" << "The calculated polarization direction is in R3 direction" << endl;
			ofs_running << "\n" << " P = " << outFormat(unit1*total_polarization, unit1*modulus, unit1*polarization_xyz) << "(e/Omega).bohr" << endl;
			ofs_running << "\n" << " P = " << outFormat(unit2*total_polarization, unit2*modulus, unit2*polarization_xyz) << "e/bohr^2" << endl;
			ofs_running << "\n" << " P = " << outFormat(unit3*total_polarization, unit3*modulus, unit3*polarization_xyz) << "C/m^2" << endl;
			ofs_running << endl;
		
			break;
		}
	
	}

	//ofs_running << "the Macroscopic_polarization is over" << endl;
	
	return;
}

string berryphase::outFormat(const double polarization, const double modulus, const Vector3<double> project)
{
	stringstream outStr;
	outStr << setw(12) << fixed << setprecision(7) << polarization << "  (mod " ;
	outStr << setw(12) << fixed << setprecision(7) << modulus << ")  (";
	outStr << setw(12) << fixed << setprecision(7) << project.x << ",";
	outStr << setw(12) << fixed << setprecision(7) << project.y << ",";
	outStr << setw(12) << fixed << setprecision(7) << project.z << ") ";
	
	return outStr.str();
}




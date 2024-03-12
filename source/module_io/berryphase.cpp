#include "berryphase.h"

#include "module_cell/klist.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

bool berryphase::berry_phase_flag=false;

berryphase::berryphase()
{
	GDIR = INPUT.gdir;
}

#ifdef __LCAO
berryphase::berryphase(Local_Orbital_wfc &lowf_in) :
    lowf(&lowf_in)
{
	GDIR = INPUT.gdir;
}
#endif

berryphase::~berryphase()
{
	//GlobalV::ofs_running << "this is ~berryphase()" << std::endl;
}

void berryphase::get_occupation_bands()
{
	double occupied_bands = static_cast<double>(GlobalV::nelec/ModuleBase::DEGSPIN);	
	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0;
	}
	
	occ_nbands = (int) occupied_bands;
	if(occ_nbands > GlobalV::NBANDS) 
	{
		ModuleBase::WARNING_QUIT("berryphase::get_occupation_bands","not enough bands for berryphase, increase band numbers.");
	}
	//GlobalV::ofs_running << "the berryphase's occ_nbands is " << occ_nbands << std::endl;
}

#ifdef __LCAO
void berryphase::lcao_init(const K_Vectors& kv)
{
	ModuleBase::TITLE("berryphase","lcao_init");
	lcao_method.init(*this->lowf->gridt, this->lowf->wfc_k_grid, kv.nkstot);
	lcao_method.cal_R_number();
	lcao_method.cal_orb_overlap();
	return;
}
#endif

// this routine depends on kpoint's mp method
void berryphase::set_kpoints(const K_Vectors& kv, const int direction)
{
	ModuleBase::TITLE("berryphase","set_kpoints");

	const int mp_x = kv.nmp[0]; // no. of kpoints along x
	const int mp_y = kv.nmp[1]; // no. of kpoints along y
	const int mp_z = kv.nmp[2]; // no. of kpoints along z
	const int num_k = int(kv.nkstot/2);	

	if( direction == 1 ) // x direction calculation
	{
		const int num_string = mp_y * mp_z;
		
		if( GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4 )
		{
			total_string = num_string;
			k_index.resize(total_string);
		}
		else if( GlobalV::NSPIN == 2 )
		{
			total_string = 2 * num_string;
			k_index.resize(total_string);
		}
		
		
		for(int istring = 0; istring < total_string; istring++)
		{
			k_index[istring].resize(mp_x+1); //adding 1 means every string from k=0 to k=G
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
		
		if( GlobalV::NSPIN == 2 )
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
		
		if( GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4 )
		{
			total_string = num_string;
			k_index.resize(total_string);
		}
		else if( GlobalV::NSPIN == 2 )
		{
			total_string = 2 * num_string;
			k_index.resize(total_string);
		}
	
		for(int istring = 0; istring < total_string; istring++)
		{
			k_index[istring].resize(mp_y+1); // adding 1 means every string from k=0 to k=G
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
		
		if( GlobalV::NSPIN == 2 )
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
		
		if( GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4 )
		{
			total_string = num_string;
			k_index.resize(total_string);
		}
		else if( GlobalV::NSPIN == 2 )
		{
			total_string = 2 * num_string;
			k_index.resize(total_string);
		}
		
		for(int istring = 0; istring < total_string; istring++)
		{
			k_index[istring].resize(mp_z+1); // adding 1 means every string from k=0 to k=G
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
		
		if( GlobalV::NSPIN == 2 )
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
	GlobalV::ofs_running << "direction is " << direction << std::endl;
	GlobalV::ofs_running << "nppstr = " << nppstr << std::endl;
	GlobalV::ofs_running << "total std::string is " << total_string << std::endl;
	for(int istring = 0; istring < total_string; istring++)
	{
		GlobalV::ofs_running << " the std::string is " << istring << std::endl;
		for(int count = 0; count < nppstr; count++)
		{
			GlobalV::ofs_running << "(" << kv.kvec_c[ k_index[istring][count] ].x << ","
							   << kv.kvec_c[ k_index[istring][count] ].y << ","
							   << kv.kvec_c[ k_index[istring][count] ].z << ")" << std::endl;
		}
		
	}
	*/
	// test by jingan
	

}

#include "../module_base/complexmatrix.h"
double berryphase::stringPhase(int index_str,
                               int nbands,
                               const int npwx,
                               const psi::Psi<std::complex<double>>* psi_in,
                               const ModulePW::PW_Basis* rhopw,
                               const ModulePW::PW_Basis_K* wfcpw,
                               const K_Vectors& kv)
{
	std::complex<double> zeta(1.0, 0.0);
	ModuleBase::ComplexMatrix mat(nbands,nbands);
	int ik_1;
	int ik_2;
	ModuleBase::Vector3<double> G(0.0,0.0,0.0);
	ModuleBase::Vector3<double> dk = kv.kvec_c[ k_index[index_str][1] ] - kv.kvec_c[ k_index[index_str][0] ];
	//GlobalV::ofs_running << "the std::string index is " << index_str << std::endl;
	
	for(int k_start = 0; k_start < (nppstr-1); k_start++)
	{
		ik_1 = k_index[index_str][k_start];
		ik_2 = k_index[index_str][k_start+1];
		
		if(GlobalV::BASIS_TYPE=="pw")
		{
			for (int mb = 0; mb < nbands; mb++)
			{

				for (int nb = 0; nb < nbands; nb++)
				{
					
					
					if(GlobalV::NSPIN!=4)
					{
						if ( k_start == (nppstr-2) )
						{
							if(direction == 1) 
							{
								ModuleBase::Vector3<double> tem_G(1.0,0.0,0.0);
								G = tem_G;
							}
							if(direction == 2)
							{
								ModuleBase::Vector3<double> tem_G(0.0,1.0,0.0);
								G = tem_G;
							}
							if(direction == 3)
							{
								ModuleBase::Vector3<double> tem_G(0.0,0.0,1.0);
								G = tem_G;
							}

                            mat(nb, mb) = pw_method.unkdotp_G0(rhopw, wfcpw, ik_1, ik_2, nb, mb, psi_in, G);
                        }
						else 
						{
                            mat(nb, mb) = pw_method.unkdotp_G(wfcpw, ik_1, ik_2, nb, mb, psi_in);
                        }
					}
					else
					{
						if ( k_start == (nppstr-2) )
						{
							if(direction == 1) 
							{
								ModuleBase::Vector3<double> tem_G(1.0,0.0,0.0);
								G = tem_G;
							}
							if(direction == 2)
							{
								ModuleBase::Vector3<double> tem_G(0.0,1.0,0.0);
								G = tem_G;
							}
							if(direction == 3)
							{
								ModuleBase::Vector3<double> tem_G(0.0,0.0,1.0);
								G = tem_G;
							}

                            mat(nb, mb) = pw_method.unkdotp_soc_G0(rhopw, wfcpw, ik_1, ik_2, nb, mb, psi_in, G);
                        }
                        else
                            mat(nb, mb) = pw_method.unkdotp_soc_G(wfcpw, ik_1, ik_2, nb, mb, npwx, psi_in);
                    }
					
				} // nb
				
			} // mb
			
			std::complex<double> det(1.0,0.0);
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
		else if(GlobalV::BASIS_TYPE=="lcao")
		{
			if(GlobalV::NSPIN!=4)
			{
				//std::complex<double> my_det = lcao_method.det_berryphase(ik_1,ik_2,dk,nbands);
				zeta = zeta * lcao_method.det_berryphase(ik_1,ik_2,dk,nbands, *this->lowf, psi_in, kv);
				// test by jingan
				//GlobalV::ofs_running << "methon 1: det = " << my_det << std::endl;
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
					
					mat(nb, mb) = lcao_method.unkdotp_LCAO(ik_1,ik_2,nb,mb,dk,kv);
				}
			}
			
			std::complex<double> det(1.0,0.0);
			int info = 0;
			int *ipiv = new int[nbands];
			LapackConnector::zgetrf(nbands, nbands, mat, nbands, ipiv, &info);				
			for (int ib = 0; ib < nbands; ib++)
			{
				if (ipiv[ib] != (ib+1)) det = -det * mat(ib,ib);
				else det = det * mat(ib,ib);
			}
			
			zeta = zeta*det;
			
			GlobalV::ofs_running << "methon 2: det = " << det << std::endl;
			
			delete[] ipiv;
			*/
			// test by jingan
		}
		#endif


	}
	
	
	return log(zeta).imag();
}

void berryphase::Berry_Phase(int nbands,
                             double& pdl_elec_tot,
                             int& mod_elec_tot,
                             const int npwx,
                             const psi::Psi<std::complex<double>>* psi_in,
                             const ModulePW::PW_Basis* rhopw,
                             const ModulePW::PW_Basis_K* wfcpw,
                             const K_Vectors& kv)
{		
	std::complex<double> cave = 0.0;
	double *phik = new double[total_string];
	double phik_ave = 0.0;
	std::complex<double> *cphik = new std::complex<double>[total_string];
	double *wistring = new double[total_string];
	
	
	// electron polarization
	
	// get weight of every std::string 
	for(int istring = 0; istring < total_string; istring++)
	{
		wistring[istring] = 1.0 / total_string;
		if(GlobalV::NSPIN == 2) wistring[istring] = wistring[istring] * 2;
	}
	
	for(int istring = 0; istring < total_string; istring++)
	{
        phik[istring] = stringPhase(istring, nbands, npwx, psi_in, rhopw, wfcpw, kv);
        // transfer phase to complex number
		cphik[istring] = std::complex<double>(cos(phik[istring]),sin(phik[istring]));	
		cave = cave + std::complex<double>(wistring[istring],0.0) * cphik[istring];
		
	}
	
	double theta0 = atan2(cave.imag(),cave.real());
	double dtheta = 0.0;
	for(int istring = 0; istring < total_string; istring++)
	{
		cphik[istring] = cphik[istring] / cave;
		dtheta = atan2(cphik[istring].imag(),cphik[istring].real());
		phik[istring] = (theta0 + dtheta) / (2 * ModuleBase::PI);
		phik_ave = phik_ave + wistring[istring] * phik[istring];
		// test by jingan
		//GlobalV::ofs_running << "phik[" << istring << "] = " << phik[istring] << std::endl;
		// test by jingan
	}
	
	if(GlobalV::NSPIN == 1)
	{
		pdl_elec_tot = 2 * phik_ave;
	}
	else if( GlobalV::NSPIN == 2 || GlobalV::NSPIN == 4 )
	{
		pdl_elec_tot = phik_ave;
	}
	
	if(GlobalV::NSPIN == 1)  // remap to [-1,1]
	{
		pdl_elec_tot = pdl_elec_tot - 2.0 * round(pdl_elec_tot/2.0);
		mod_elec_tot = 2;
	}
	else if( GlobalV::NSPIN == 2 || GlobalV::NSPIN == 4 )  // remap to [-0.5,0.5]
	{
		pdl_elec_tot = pdl_elec_tot - 1.0 * round(pdl_elec_tot/1.0);
		mod_elec_tot = 1;
	}
	
	
	
	delete[] phik;
	delete[] cphik;
	delete[] wistring;
	
	
	//GlobalV::ofs_running << "Berry_Phase end " << std::endl;

}

void berryphase::Macroscopic_polarization(const int npwx,
                                          const psi::Psi<std::complex<double>>* psi_in,
                                          const ModulePW::PW_Basis* rhopw,
                                          const ModulePW::PW_Basis_K* wfcpw,
                                          const K_Vectors& kv)
{	
	get_occupation_bands();
#ifdef __LCAO	
	if( GlobalV::BASIS_TYPE == "lcao" ) this->lcao_init(kv);
#endif
	
	GlobalV::ofs_running << "\n\n\n\n";
	GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " | POLARIZATION CALCULATION:                                          |" << std::endl;
	GlobalV::ofs_running << " |                  Modern Theory of Polarization                     |" << std::endl;
	GlobalV::ofs_running << " | calculate the Macroscopic polarization of a crystalline insulator  |" << std::endl;
	GlobalV::ofs_running << " | by using Berry Phase method.                                       |" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	GlobalV::ofs_running << "\n\n\n\n";
	
	
	// ion polarization	
	double polarization_ion[3]; // means three lattice vector directions R1，R2，R3
	ModuleBase::GlobalFunc::ZEROS(polarization_ion,3);
	// reciprocal lattice
	ModuleBase::Vector3<double> rcell_1(GlobalC::ucell.G.e11,GlobalC::ucell.G.e12,GlobalC::ucell.G.e13);
	ModuleBase::Vector3<double> rcell_2(GlobalC::ucell.G.e21,GlobalC::ucell.G.e22,GlobalC::ucell.G.e23);
	ModuleBase::Vector3<double> rcell_3(GlobalC::ucell.G.e31,GlobalC::ucell.G.e32,GlobalC::ucell.G.e33);
	int *mod_ion = new int[GlobalC::ucell.nat];
	double *pdl_ion_R1 = new double[GlobalC::ucell.nat];
	double *pdl_ion_R2 = new double[GlobalC::ucell.nat];
	double *pdl_ion_R3 = new double[GlobalC::ucell.nat];  
	ModuleBase::GlobalFunc::ZEROS(mod_ion,GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(pdl_ion_R1,GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(pdl_ion_R2,GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(pdl_ion_R3,GlobalC::ucell.nat);
	
	bool lodd = false;
	int atom_index = 0;
	for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
		for(int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
		{
			if(GlobalC::ucell.atoms[it].ncpp.zv % 2 == 1)
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
	for(int it = 0; it < GlobalC::ucell.ntype; it++)
	{
		for(int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
		{
			pdl_ion_R1[atom_index] = GlobalC::ucell.atoms[it].ncpp.zv * (GlobalC::ucell.atoms[it].tau[ia] * rcell_1);
			pdl_ion_R2[atom_index] = GlobalC::ucell.atoms[it].ncpp.zv * (GlobalC::ucell.atoms[it].tau[ia] * rcell_2);
			pdl_ion_R3[atom_index] = GlobalC::ucell.atoms[it].ncpp.zv * (GlobalC::ucell.atoms[it].tau[ia] * rcell_3);
			atom_index++;
		}
	}
	

	
	for(int i = 0; i < GlobalC::ucell.nat; i++)
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
	if( (!lodd) && (GlobalV::NSPIN==1) ) modulus = 2;
	else modulus = 1;
	
	// test by jingan
	//GlobalV::ofs_running << "ion polarization end" << std::endl;
	// test by jingan


	// berry phase calculate begin	
	switch(GDIR)
	{
		case 1:
		{
			direction = 1;
			set_kpoints(kv, direction);
			double pdl_elec_tot = 0.0;
			int mod_elec_tot = 0;
            Berry_Phase(occ_nbands, pdl_elec_tot, mod_elec_tot, npwx, psi_in, rhopw, wfcpw, kv);

            const double rmod = GlobalC::ucell.a1.norm() * GlobalC::ucell.lat0;
			const double unit1 = rmod;
			const double unit2 = rmod / GlobalC::ucell.omega;
			const double unit3 = ( rmod / GlobalC::ucell.omega ) * ( 1.60097e-19/pow(5.29177e-11,2) );
			
			GlobalV::ofs_running << " VALUES OF POLARIZATION" << std::endl;
			GlobalV::ofs_running << std::endl;
			GlobalV::ofs_running << "  The Ionic Phase: " << std::setw(10) << std::fixed << std::setprecision(5) << polarization_ion[0] << std::endl;
			GlobalV::ofs_running << " Electronic Phase: " << std::setw(10) << std::fixed << std::setprecision(5) << pdl_elec_tot << std::endl;
			//GlobalV::ofs_running << " the electronic part polarization is P(ele) = " << rmod * pdl_elec_tot << "   (e/Omega).bohr   in R1 direction" << std::endl;
			
			// calculate total polarization,add electron part and ions part
			double total_polarization = pdl_elec_tot + polarization_ion[0] ;
			
			ModuleBase::Vector3<double> polarization_xyz = GlobalC::ucell.a1;
			polarization_xyz.normalize();
			polarization_xyz = total_polarization * polarization_xyz;

			GlobalV::ofs_running << "\n" << "The calculated polarization direction is in R1 direction" << std::endl;
			GlobalV::ofs_running << "\n" << " P = " << outFormat(unit1*total_polarization, unit1*modulus, unit1*polarization_xyz) << "(e/Omega).bohr" << std::endl;
			GlobalV::ofs_running << "\n" << " P = " << outFormat(unit2*total_polarization, unit2*modulus, unit2*polarization_xyz) << "e/bohr^2" << std::endl;
			GlobalV::ofs_running << "\n" << " P = " << outFormat(unit3*total_polarization, unit3*modulus, unit3*polarization_xyz) << "C/m^2" << std::endl;			   
			GlobalV::ofs_running << std::endl;
		
			break;
		
		}
		case 2:
		{
			direction = 2;
			set_kpoints(kv, direction);
			double pdl_elec_tot = 0.0;
			int mod_elec_tot = 0;
            Berry_Phase(occ_nbands, pdl_elec_tot, mod_elec_tot, npwx, psi_in, rhopw, wfcpw, kv);

            const double rmod = GlobalC::ucell.a2.norm() * GlobalC::ucell.lat0;
			const double unit1 = rmod;
			const double unit2 = rmod / GlobalC::ucell.omega;
			const double unit3 = ( rmod / GlobalC::ucell.omega ) * ( 1.60097e-19/pow(5.29177e-11,2) );
			
			GlobalV::ofs_running << " VALUES OF POLARIZATION" << std::endl;
			GlobalV::ofs_running << std::endl;
			GlobalV::ofs_running << "  The Ionic Phase: " << std::setw(10) << std::fixed << std::setprecision(5) << polarization_ion[1] << std::endl;
			GlobalV::ofs_running << " Electronic Phase: " << std::setw(10) << std::fixed << std::setprecision(5) << pdl_elec_tot << std::endl;
			//GlobalV::ofs_running << " the electronic part polarization is P(ele) = " << rmod * pdl_elec_tot << "   (e/Omega).bohr   in R2 direction" << std::endl;
		
			// calculate total polarization,add electron part and ions part
			double total_polarization = pdl_elec_tot + polarization_ion[1] ;
			
			ModuleBase::Vector3<double> polarization_xyz = GlobalC::ucell.a2;
			polarization_xyz.normalize();
			polarization_xyz = total_polarization * polarization_xyz;
		
			GlobalV::ofs_running << "\n" << "The calculated polarization direction is in R2 direction" << std::endl;
			GlobalV::ofs_running << "\n"  << " P = " << outFormat(unit1*total_polarization, unit1*modulus, unit1*polarization_xyz) << "(e/Omega).bohr" << std::endl;
			GlobalV::ofs_running << "\n"  << " P = " << outFormat(unit2*total_polarization, unit2*modulus, unit2*polarization_xyz) << "e/bohr^2" << std::endl;
			GlobalV::ofs_running << "\n"  << " P = " << outFormat(unit3*total_polarization, unit3*modulus, unit3*polarization_xyz) << "C/m^2" << std::endl;
			GlobalV::ofs_running << std::endl;
			
			break;
		
		}
		case 3:
		{
			direction = 3;
			set_kpoints(kv, direction);
			double pdl_elec_tot = 0.0;
			int mod_elec_tot = 0;
            Berry_Phase(occ_nbands, pdl_elec_tot, mod_elec_tot, npwx, psi_in, rhopw, wfcpw, kv);

            const double rmod = GlobalC::ucell.a3.norm() * GlobalC::ucell.lat0;
			const double unit1 = rmod;
			const double unit2 = rmod / GlobalC::ucell.omega;
			const double unit3 = ( rmod / GlobalC::ucell.omega ) * ( 1.60097e-19/pow(5.29177e-11,2) );
			
			GlobalV::ofs_running << " VALUES OF POLARIZATION" << std::endl;
			GlobalV::ofs_running << std::endl;
			GlobalV::ofs_running << "  The Ionic Phase: " << std::setw(10) << std::fixed << std::setprecision(5) << polarization_ion[2] << std::endl;
			GlobalV::ofs_running << " Electronic Phase: " << std::setw(10) << std::fixed << std::setprecision(5) << pdl_elec_tot << std::endl;
			//GlobalV::ofs_running << " the electronic part polarization is P(ele) = " << rmod * pdl_elec_tot << "   (e/Omega).bohr   in R3 direction" << std::endl;
		
			// calculate total polarization,add electron part and ions part
			double total_polarization = pdl_elec_tot + polarization_ion[2] ;
			
			ModuleBase::Vector3<double> polarization_xyz = GlobalC::ucell.a3;
			polarization_xyz.normalize();
			polarization_xyz = total_polarization * polarization_xyz;

			GlobalV::ofs_running << "\n" << "The calculated polarization direction is in R3 direction" << std::endl;
			GlobalV::ofs_running << "\n" << " P = " << outFormat(unit1*total_polarization, unit1*modulus, unit1*polarization_xyz) << "(e/Omega).bohr" << std::endl;
			GlobalV::ofs_running << "\n" << " P = " << outFormat(unit2*total_polarization, unit2*modulus, unit2*polarization_xyz) << "e/bohr^2" << std::endl;
			GlobalV::ofs_running << "\n" << " P = " << outFormat(unit3*total_polarization, unit3*modulus, unit3*polarization_xyz) << "C/m^2" << std::endl;
			GlobalV::ofs_running << std::endl;
		
			break;
		}
	
	}

	//GlobalV::ofs_running << "the Macroscopic_polarization is over" << std::endl;
	
	return;
}

std::string berryphase::outFormat(const double polarization, const double modulus, const ModuleBase::Vector3<double> project)
{
	std::stringstream outStr;
	outStr << std::setw(12) << std::fixed << std::setprecision(7) << polarization << "  (mod " ;
	outStr << std::setw(12) << std::fixed << std::setprecision(7) << modulus << ")  (";
	outStr << std::setw(12) << std::fixed << std::setprecision(7) << project.x << ",";
	outStr << std::setw(12) << std::fixed << std::setprecision(7) << project.y << ",";
	outStr << std::setw(12) << std::fixed << std::setprecision(7) << project.z << ") ";
	
	return outStr.str();
}




#include "dos.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#ifdef __LCAO
void Dos::calculate_Mulliken(const std::string &fa, Gint_Gamma &gg)
{
	ModuleBase::TITLE("Dos","calculate_Mulliken");
	std::ofstream ofs;
	
	if(GlobalV::MY_RANK==0)
	{
		ofs.open(fa.c_str());
		ofs << std::setiosflags(ios::left);
	}

	GlobalV::ofs_running << "\n CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;

	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		double** mulliken = new double* [GlobalV::NSPIN];
		for(int is=0; is<GlobalV::NSPIN; ++is)
		{
			mulliken[is] = new double[GlobalV::NLOCAL];
			ModuleBase::GlobalFunc::ZEROS(mulliken[is], GlobalV::NLOCAL);
		}
		
		gg.cal_mulliken( mulliken );	

		if(GlobalV::MY_RANK==0)
		{
			// normalize the mulliken charge.
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
				{
					mulliken[is][iw] *= GlobalC::ucell.omega/GlobalC::pw.ncxyz;
					if( abs(mulliken[is][iw]) < 1.0e-10 ) mulliken[is][iw] = 0.0; 
				}
			}

			// calculate the total charge of the system.
			double sch = 0.0;
			ofs << std::setprecision(8);
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				double sss = 0.0;
				for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
				{
					sch += mulliken[is][iw];
					sss += mulliken[is][iw];
				}
				ofs << sss << " (Total charge all spin " << is+1 << ")" << std::endl;
			}
			ofs << sch << " (Total charge of the system)" << std::endl;
			 
			 // output information for each atom.
			int iw_all=0;
			for(int it=0; it<GlobalC::ucell.ntype; ++it)
			{
				Atom* atom = &GlobalC::ucell.atoms[it];
				ofs << std::setw(5) << "TYPE" << std::setw(8) << "ATOM" << std::setw(5) << "SPIN";
				ofs << std::setprecision(3);
				for(int l=0; l<= atom->nwl; ++l)
				{
					for(int m=0; m<2*l+1; ++m)
					{
						if(l==0) ofs << std::setw(12) << "s";
						else if(l==1){ std::stringstream ss;ss << "p" << m+1;ofs << std::setw(12) << ss.str(); }
						else if(l==2){ std::stringstream ss;ss << "d" << m+1;ofs << std::setw(12) << ss.str(); }
						else if(l==3){ std::stringstream ss;ss << "f" << m+1;ofs << std::setw(12) << ss.str(); }
						else if(l==4){ std::stringstream ss;ss << "g" << m+1;ofs << std::setw(12) << ss.str(); }
					}
				} 

				ofs << std::setw(12) << "sum/zv";
				ofs << std::endl;

				double scht = 0.0;
				for(int ia=0; ia<atom->na; ++ia)
				{
					for(int is=0; is<GlobalV::NSPIN; is++)
					{
						int iw_alllll = iw_all;
						
						double sum = 0.0;
						ofs << std::setw(5) << atom->label 
							<< std::setw(8) << ia+1 << std::setw(5) << is+1;
						for(int l=0; l<=atom->nwl; ++l)
						{
							// sum up the multi-zeta charge.
							double *mmm = new double[2*l+1];
							ModuleBase::GlobalFunc::ZEROS(mmm, 2*l+1);
							for(int n=0; n<atom->l_nchi[l]; ++n)
							{
								for(int m=0; m<2*l+1; ++m)
								{
									mmm[m] += mulliken[is][iw_alllll];
									++iw_alllll;
								}
							}

							for(int m=0; m<2*l+1; ++m)
							{
								ofs << std::setw(12) << mmm[m];
								sum += mmm[m];
							}
							delete[] mmm;
						}
				
						ofs << sum << "/" << atom->zv/GlobalV::NSPIN;
						ofs << std::endl;

						scht += sum;
					}

					iw_all += atom->nw;
				}

				ofs << std::setprecision(8);
				ofs << scht << " (Total charge of atom species " << atom->label << ")" << std::endl;
			}
		}

		for(int is=0; is<GlobalV::NSPIN; ++is)
		{
			delete[] mulliken[is];
		}	
		delete[] mulliken;
	}
	else
	{
		ModuleBase::WARNING_QUIT("Mulliken Charge","Not implement yet.");	
	}	
	

	if(GlobalV::MY_RANK==0) ofs.close();

	return;
}
#endif

bool Dos::calculate_dos
(
	const int &is,
	const std::vector<int> &isk,
	const std::string &fa, //file address for DOS
	const std::string &fa1, //file address for DOS_smearing
	const double &de_ev, // delta energy in ev
	const double &emax_ev,
	const double &emin_ev,// minimal energy in ev.
	const int &nks,//number of k points
	const int &nkstot,
	const std::vector<double> &wk,//weight of k points
	const ModuleBase::matrix &wg,//weight of (kpoint,bands)
	const int &nbands,// number of bands
	double** ekb//store energy for each k point and each band
)
{
	ModuleBase::TITLE("Dos","calculae_dos");
	std::ofstream ofs;
	std::ofstream ofs1;
	if(GlobalV::MY_RANK==0)
	{
		ofs.open(fa.c_str());//make the file clear!!
		ofs1.open(fa1.c_str());//make the file clear!!
	}
	std::vector<double> dos;
	std::vector<double> ene;
	std::vector<double> dos_smearing; //dos_smearing

#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(de_ev <= 0)
	{
		ModuleBase::WARNING("DOS::calculate_dos","de <= 0 ");
		return 0; 
	}
	else if(emax_ev < emin_ev)
	{
		ModuleBase::WARNING("calculate_dos","emax_ev < emin_ev");
		return 0;
	}

	// mohan fixed bug 2010-1-18
	const int npoints = static_cast<int>(std::floor ( ( emax_ev - emin_ev ) / de_ev )) ;
	dos.clear();
	ene.clear();
	if(npoints <= 0)
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npoints",npoints);
		ModuleBase::WARNING("calculate_dos","npoints < 0");
		return 0;
	}
	if(GlobalV::MY_RANK==0)
	{
		ofs << npoints << std::endl;
		ofs << GlobalC::kv.nkstot << std::endl;
	}

	GlobalV::ofs_running << "\n OUTPUT DOS FILE IN: " << fa << std::endl;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"min state energy (eV)",emin_ev);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max state energy (eV)",emax_ev);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"delta energy interval (eV)",  de_ev);
	
	double *e_mod = new double[npoints];
	ModuleBase::GlobalFunc::ZEROS(e_mod, npoints);

	double sum   = 0.0;
	double e_new = emin_ev;
	double e_old = 0.0;
	while( e_new < emax_ev)
	{
//		GlobalV::ofs_running << " enew=" << e_new << std::endl;
		double count = 0.0;
		e_old = e_new ;
		e_new += de_ev;
		for(int ik=0;ik<nks;ik++)
		{
			if(is == GlobalC::kv.isk[ik])
			{
				for(int ib = 0; ib < nbands; ib++)
				{
					//  compare et and e_old(e_new) in ev unit.
					if( ekb[ik][ib]*ModuleBase::Ry_to_eV >= e_old && ekb[ik][ib]*ModuleBase::Ry_to_eV < e_new)
					{
						// because count is 'double' type,so
						// we can't write count++ or ++count
						count += wk[ik]*nkstot; //mohanix bug 2012-04-23
//						GlobalV::ofs_running << " count = " << count << " wk = " << wk[ik] << " nks = " << nks << std::endl;
					}		
				}
			}
		}
#ifdef __MPI
		Parallel_Reduce::reduce_double_allpool(count);
#endif
		count = count / static_cast<double>(nkstot);
		sum += count;
		if(GlobalV::MY_RANK==0)
		{
			ofs << e_new << " " << count << std::endl;
			dos.push_back(count);
			ene.push_back(e_new);
		}

	}

	//now use Gaussian smearing to smooth the dos and write to DOS_is_smearing
	if(GlobalV::MY_RANK==0)
	{
		dos_smearing.resize(dos.size()-1);

		//double b = INPUT.b_coef;
		double b = sqrt(2.0)*GlobalC::en.bcoeff;
		for(int i=0;i<dos.size()-1;i++)
		{
			double Gauss=0.0;

			for(int j=0;j<dos.size()-1;j++)
			{
				double de = ene[j] - ene[i];
				double de2 = de * de;
				//----------------------------------------------------------
				// EXPLAIN : if en
				//----------------------------------------------------------
				Gauss = exp(-de2/b/b)/sqrt(3.1415926)/b;
				dos_smearing[j] += dos[i]*Gauss;
			}
		}

		//----------------------------------------------------------
		// EXPLAIN : output DOS_smearing.dat
		//----------------------------------------------------------
		double sum2=0.0;
		for(int i=0;i<dos.size()-1;i++)
		{
			sum2 += dos_smearing[i];
			ofs1 <<std::setw(20)<<ene[i]
				<<std::setw(20)<<dos_smearing[i]
				<<std::setw(20)<<sum2<<"\n";
		}
	}


	if(GlobalV::MY_RANK==0)
	{
		ofs.close();
		ofs1.close();
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of bands",nbands);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"sum up the states", sum);
	delete[] e_mod;

	return 1;
}


void Dos::nscf_fermi_surface(const std::string &out_band_dir,
	const int &nks,
	const int &nband,
	double **ekb)
{
#ifdef __MPI

	int start = 1;
	int end = GlobalV::NBANDS;

	assert(GlobalC::wf.allocate_ekb);

	std::ofstream ofs;
	if(GlobalV::MY_RANK==0)
	{
		ofs.open(out_band_dir.c_str());//make the file clear!!
		ofs << std::setprecision(6);
		ofs.close();	
	}

	for(int ik=0; ik<GlobalC::kv.nkstot; ik++)
	{
		if ( GlobalV::MY_POOL == GlobalC::Pkpoints.whichpool[ik] )
		{
			if( GlobalV::RANK_IN_POOL == 0)
			{
				std::ofstream ofs(out_band_dir.c_str(),ios::app);
				ofs << std::setprecision(8);

				if(ik==0)
				{
					ofs << " BEGIN_INFO" << std::endl;
					ofs << "   #" << std::endl;
					ofs << "   # this is a Band-XCRYSDEN-Structure-File" << std::endl;
					ofs << "   # aimed at Visualization of Fermi Surface" << std::endl;
					ofs << "   #" << std::endl;
					ofs << "   # Case: " << GlobalC::ucell.latName << std::endl;
					ofs << "   #" << std::endl;	
					ofs << " Fermi Energy: " << GlobalC::en.ef << std::endl;
					ofs << " END_INFO" << std::endl;
					ofs << " BEGIN_BLOCK_BANDGRID_3D" << std::endl;
					ofs << " band_energies" << std::endl;
					ofs << " BANDGRID_3D_BANDS" << std::endl;
					ofs << " " << end-start+1 << std::endl;
					ofs << " NKX NKY NKZ" << std::endl;
					ofs << " 0 0 0" << std::endl;
					ofs << " " << GlobalC::ucell.G.e11 << " " << GlobalC::ucell.G.e12 << " " << GlobalC::ucell.G.e13 << std::endl; 
					ofs << " " << GlobalC::ucell.G.e21 << " " << GlobalC::ucell.G.e22 << " " << GlobalC::ucell.G.e23 << std::endl; 
					ofs << " " << GlobalC::ucell.G.e31 << " " << GlobalC::ucell.G.e32 << " " << GlobalC::ucell.G.e33 << std::endl; 
				}

				const int ik_now = ik - GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL];
				ofs << "ik= " << ik << std::endl;
				ofs << GlobalC::kv.kvec_c[ik_now].x << " " << GlobalC::kv.kvec_c[ik_now].y << " " << GlobalC::kv.kvec_c[ik_now].z << std::endl;  

				for(int ib = 0; ib < nband; ib++)
				{
					ofs << " " << ekb[ik_now][ib] * ModuleBase::Ry_to_eV;
				}
				ofs << std::endl;

				// the last k point
				if(ik==GlobalC::kv.nkstot-1)
				{
					ofs << " END_BANDGRID_3D" << std::endl;
					ofs << " END_BLOCK_BANDGRID_3D" << std::endl;
				}
				ofs.close();

			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

#else


#endif

}


void Dos::nscf_band(
	const int &is,
	const std::string &out_band_dir, 
	const int &nks, 
	const int &nband,
	const double &fermie,
	double** ekb)
{
	ModuleBase::TITLE("Dos","nscf_band");

#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
		std::ofstream ofs(out_band_dir.c_str());//make the file clear!!
		ofs.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<double> klength;
	klength.resize(nks);
	klength[0] = 0.0;
	for(int ik=0; ik<nks; ik++)
	{
		if (ik>0)
		{
			auto delta=GlobalC::kv.kvec_c[ik]-GlobalC::kv.kvec_c[ik-1];
			klength[ik] = klength[ik-1] + delta.norm();
		}
		if ( GlobalV::MY_POOL == GlobalC::Pkpoints.whichpool[ik] )
		{
			const int ik_now = ik - GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL];
			if( GlobalC::kv.isk[ik_now+is*nks] == is )
			{ 
				if ( GlobalV::RANK_IN_POOL == 0)
				{
					std::ofstream ofs(out_band_dir.c_str(),ios::app);
					ofs << std::setprecision(8);
					//start from 1
					ofs << ik+1;
					ofs << " " << klength[ik] << " ";
					for(int ib = 0; ib < nband; ib++)
					{
						ofs << " " << (ekb[ik_now+is*nks][ib]-fermie) * ModuleBase::Ry_to_eV;
					}
					ofs << std::endl;
					ofs.close();	
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// old version
	/*
	for(int ip=0;ip<GlobalV::KPAR;ip++)
	{
		if(GlobalV::MY_POOL == ip && GlobalV::RANK_IN_POOL == 0)
		{
			std::ofstream ofs(out_band_dir.c_str(),ios::app);
			for(int ik=0;ik<nks;ik++)
			{
				ofs<<std::setw(12)<<ik;
				for(int ib = 0; ib < nband; ib++)
				{
					ofs <<std::setw(12)<< ekb[ik][ib] * ModuleBase::Ry_to_eV;
				}
				ofs<<std::endl;
			}
			ofs.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	*/
#else
//	std::cout<<"\n nband = "<<nband<<std::endl;
//	std::cout<<out_band_dir<<std::endl;

	std::ofstream ofs(out_band_dir.c_str());
	for(int ik=0;ik<nks;ik++)
	{
		if( GlobalC::kv.isk[ik] == is)
		{
			ofs<<std::setw(12)<<ik + 1;
			for(int ibnd = 0; ibnd < nband; ibnd++)
			{
				ofs <<std::setw(15) << (ekb[ik][ibnd]-fermie) * ModuleBase::Ry_to_eV;
			}
			ofs<<std::endl;
		}
	}
	ofs.close();
#endif
	return;
}

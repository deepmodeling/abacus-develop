#include "dos.h"
#include "../src_pw/global.h"
#ifdef __LCAO
void Dos::calculate_Mulliken(const string &fa)
{
	TITLE("Dos","calculate_Mulliken");
	ofstream ofs;
	
	if(MY_RANK==0)
	{
		ofs.open(fa.c_str());
		ofs << setiosflags(ios::left);
	}

	ofs_running << "\n CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << endl;

	if(GAMMA_ONLY_LOCAL)
	{
		double** mulliken = new double* [NSPIN];
		for(int is=0; is<NSPIN; ++is)
		{
			mulliken[is] = new double[NLOCAL];
			ZEROS(mulliken[is], NLOCAL);
		}
		
		UHM.GG.cal_mulliken( mulliken );	

		if(MY_RANK==0)
		{
			// normalize the mulliken charge.
			for(int is=0; is<NSPIN; ++is)
			{
				for(int iw=0; iw<NLOCAL; ++iw)
				{
					mulliken[is][iw] *= ucell.omega/pw.ncxyz;
					if( abs(mulliken[is][iw]) < 1.0e-10 ) mulliken[is][iw] = 0.0; 
				}
			}

			// calculate the total charge of the system.
			double sch = 0.0;
			ofs << setprecision(8);
			for(int is=0; is<NSPIN; ++is)
			{
				double sss = 0.0;
				for(int iw=0; iw<NLOCAL; ++iw)
				{
					sch += mulliken[is][iw];
					sss += mulliken[is][iw];
				}
				ofs << sss << " (Total charge all spin " << is+1 << ")" << endl;
			}
			ofs << sch << " (Total charge of the system)" << endl;
			 
			 // output information for each atom.
			int iw_all=0;
			for(int it=0; it<ucell.ntype; ++it)
			{
				Atom* atom = &ucell.atoms[it];
				ofs << setw(5) << "TYPE" << setw(8) << "ATOM" << setw(5) << "SPIN";
				ofs << setprecision(3);
				for(int l=0; l<= atom->nwl; ++l)
				{
					for(int m=0; m<2*l+1; ++m)
					{
						if(l==0) ofs << setw(12) << "s";
						else if(l==1){ stringstream ss;ss << "p" << m+1;ofs << setw(12) << ss.str(); }
						else if(l==2){ stringstream ss;ss << "d" << m+1;ofs << setw(12) << ss.str(); }
						else if(l==3){ stringstream ss;ss << "f" << m+1;ofs << setw(12) << ss.str(); }
						else if(l==4){ stringstream ss;ss << "g" << m+1;ofs << setw(12) << ss.str(); }
					}
				} 

				ofs << setw(12) << "sum/zv";
				ofs << endl;

				double scht = 0.0;
				for(int ia=0; ia<atom->na; ++ia)
				{
					for(int is=0; is<NSPIN; is++)
					{
						int iw_alllll = iw_all;
						
						double sum = 0.0;
						ofs << setw(5) << atom->label 
							<< setw(8) << ia+1 << setw(5) << is+1;
						for(int l=0; l<=atom->nwl; ++l)
						{
							// sum up the multi-zeta charge.
							double *mmm = new double[2*l+1];
							ZEROS(mmm, 2*l+1);
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
								ofs << setw(12) << mmm[m];
								sum += mmm[m];
							}
							delete[] mmm;
						}
				
						ofs << sum << "/" << atom->zv/NSPIN;
						ofs << endl;

						scht += sum;
					}

					iw_all += atom->nw;
				}

				ofs << setprecision(8);
				ofs << scht << " (Total charge of atom species " << atom->label << ")" << endl;
			}
		}

		for(int is=0; is<NSPIN; ++is)
		{
			delete[] mulliken[is];
		}	
		delete[] mulliken;
	}
	else
	{
		WARNING_QUIT("Mulliken Charge","Not implement yet.");	
	}	
	

	if(MY_RANK==0) ofs.close();

	return;
}
#endif

bool Dos::calculate_dos
(
	const int &is,
	const int *isk,
	const string &fa, //file address
	const double &de_ev, // delta energy in ev
	const double &emax_ev,
	const double &emin_ev,// minimal energy in ev.
	const int &nks,//number of k points
	const int &nkstot,
	const double *wk,//weight of k points
	const matrix &wg,//weight of (kpoint,bands)
	const int &nbands,// number of bands
	double** ekb//store energy for each k point and each band
)
{
	TITLE("Dos","calculae_dos");
	ofstream ofs;
	if(MY_RANK==0)
	{
		ofs.open(fa.c_str());//make the file clear!!
	}

#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(de_ev <= 0)
	{
		WARNING("DOS::calculate_dos","de <= 0 ");
		return 0; 
	}
	else if(emax_ev < emin_ev)
	{
		WARNING("calculate_dos","emax_ev < emin_ev");
		return 0;
	}

	// mohan fixed bug 2010-1-18
	const int npoints = static_cast<int>(std::floor ( ( emax_ev - emin_ev ) / de_ev )) ;
	if(npoints <= 0)
	{
		OUT(ofs_running,"npoints",npoints);
		WARNING("calculate_dos","npoints < 0");
		return 0;
	}
	if(MY_RANK==0)
	{
		ofs << npoints << endl;
		ofs << kv.nkstot << endl;
	}

	ofs_running << "\n OUTPUT DOS FILE IN: " << fa << endl;
	OUT(ofs_running,"min state energy (eV)",emin_ev);
	OUT(ofs_running,"max state energy (eV)",emax_ev);
	OUT(ofs_running,"delta energy interval (eV)",  de_ev);
	
	double *e_mod = new double[npoints];
	ZEROS(e_mod, npoints);

	double sum   = 0.0;
	double e_new = emin_ev;
	double e_old = 0.0;
	while( e_new < emax_ev)
	{
//		ofs_running << " enew=" << e_new << endl;
		double count = 0.0;
		e_old = e_new ;
		e_new += de_ev;
		for(int ik=0;ik<nks;ik++)
		{
			if(is == kv.isk[ik])
			{
				for(int ib = 0; ib < nbands; ib++)
				{
					//  compare et and e_old(e_new) in ev unit.
					if( ekb[ik][ib]*Ry_to_eV >= e_old && ekb[ik][ib]*Ry_to_eV < e_new)
					{
						// because count is 'double' type,so
						// we can't write count++ or ++count
						count += wk[ik]*nkstot; //mohanix bug 2012-04-23
//						ofs_running << " count = " << count << " wk = " << wk[ik] << " nks = " << nks << endl;
					}		
				}
			}
		}
#ifdef __MPI
		Parallel_Reduce::reduce_double_allpool(count);
#endif
		count = count / static_cast<double>(nkstot);
		sum += count;
		if(MY_RANK==0)
		{
			ofs << e_new << " " << count << endl;
		}

	}
	if(MY_RANK==0)
	{
		ofs.close();
	}
	OUT(ofs_running,"number of bands",nbands);
	OUT(ofs_running,"sum up the states", sum);
	delete[] e_mod;

	return 1;
}


void Dos::nscf_fermi_surface(const string &out_band_dir,
	const int &nks,
	const int &nband,
	double **ekb)
{
#ifdef __MPI

	int start = 1;
	int end = NBANDS;

	assert(wf.allocate_ekb);

	ofstream ofs;
	if(MY_RANK==0)
	{
		ofs.open(out_band_dir.c_str());//make the file clear!!
		ofs << setprecision(6);
		ofs.close();	
	}

	for(int ik=0; ik<kv.nkstot; ik++)
	{
		if ( MY_POOL == Pkpoints.whichpool[ik] )
		{
			if( RANK_IN_POOL == 0)
			{
				ofstream ofs(out_band_dir.c_str(),ios::app);
				ofs << setprecision(8);

				if(ik==0)
				{
					ofs << " BEGIN_INFO" << endl;
					ofs << "   #" << endl;
					ofs << "   # this is a Band-XCRYSDEN-Structure-File" << endl;
					ofs << "   # aimed at Visualization of Fermi Surface" << endl;
					ofs << "   #" << endl;
					ofs << "   # Case: " << ucell.latName << endl;
					ofs << "   #" << endl;	
					ofs << " Fermi Energy: " << en.ef << endl;
					ofs << " END_INFO" << endl;
					ofs << " BEGIN_BLOCK_BANDGRID_3D" << endl;
					ofs << " band_energies" << endl;
					ofs << " BANDGRID_3D_BANDS" << endl;
					ofs << " " << end-start+1 << endl;
					ofs << " NKX NKY NKZ" << endl;
					ofs << " 0 0 0" << endl;
					ofs << " " << ucell.G.e11 << " " << ucell.G.e12 << " " << ucell.G.e13 << endl; 
					ofs << " " << ucell.G.e21 << " " << ucell.G.e22 << " " << ucell.G.e23 << endl; 
					ofs << " " << ucell.G.e31 << " " << ucell.G.e32 << " " << ucell.G.e33 << endl; 
				}

				const int ik_now = ik - Pkpoints.startk_pool[MY_POOL];
				ofs << "ik= " << ik << endl;
				ofs << kv.kvec_c[ik_now].x << " " << kv.kvec_c[ik_now].y << " " << kv.kvec_c[ik_now].z << endl;  

				for(int ib = 0; ib < nband; ib++)
				{
					ofs << " " << ekb[ik_now][ib] * Ry_to_eV;
				}
				ofs << endl;

				// the last k point
				if(ik==kv.nkstot-1)
				{
					ofs << " END_BANDGRID_3D" << endl;
					ofs << " END_BLOCK_BANDGRID_3D" << endl;
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
	const string &out_band_dir, 
	const int &nks, 
	const int &nband,
	const double &fermie,
	double** ekb)
{
	TITLE("Dos","nscf_band");

#ifdef __MPI
	if(MY_RANK==0)
	{
		ofstream ofs(out_band_dir.c_str());//make the file clear!!
		ofs.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(int ik=0; ik<nks; ik++)
	{
		if ( MY_POOL == Pkpoints.whichpool[ik] )
		{
			const int ik_now = ik - Pkpoints.startk_pool[MY_POOL];
			if( kv.isk[ik_now+is*nks] == is )
			{ 
				if ( RANK_IN_POOL == 0)
				{
					ofstream ofs(out_band_dir.c_str(),ios::app);
					ofs << setprecision(8);
					//start from 1
					ofs << ik+1;
					for(int ib = 0; ib < nband; ib++)
					{
						ofs << " " << (ekb[ik_now+is*nks][ib]-fermie) * Ry_to_eV;
					}
					ofs << endl;
					ofs.close();	
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// old version
	/*
	for(int ip=0;ip<NPOOL;ip++)
	{
		if(MY_POOL == ip && RANK_IN_POOL == 0)
		{
			ofstream ofs(out_band_dir.c_str(),ios::app);
			for(int ik=0;ik<nks;ik++)
			{
				ofs<<setw(12)<<ik;
				for(int ib = 0; ib < nband; ib++)
				{
					ofs <<setw(12)<< ekb[ik][ib] * Ry_to_eV;
				}
				ofs<<endl;
			}
			ofs.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	*/
#else
//	cout<<"\n nband = "<<nband<<endl;
//	cout<<out_band_dir<<endl;

	ofstream ofs(out_band_dir.c_str());
	for(int ik=0;ik<nks;ik++)
	{
		if( kv.isk[ik] == is)
		{
			ofs<<setw(12)<<ik;
			for(int ibnd = 0; ibnd < nband; ibnd++)
			{
				ofs <<setw(15) << (ekb[ik][ibnd]-fermie) * Ry_to_eV;
			}
			ofs<<endl;
		}
	}
	ofs.close();
#endif
	return;
}

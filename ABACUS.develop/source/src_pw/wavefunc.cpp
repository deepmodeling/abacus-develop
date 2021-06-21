#include "wavefunc.h"
#include "global.h"
#include "../src_pw/wavefunc_in_pw.h"
#include "src_io/winput.h"
#include "src_io/chi0_hilbert.h"

wavefunc::wavefunc()
{
	allocate_ekb = false;
}

wavefunc::~wavefunc()
{ 
	if(test_deconstructor)
	{
		cout << " ~wavefunc()" << endl;
	}
	if(allocate_ekb)
	{
		// bug still remains, hard to find!
		// it might be somewhere out there,
		// may be in diagH_LAPACK.
		// I don't know why.......
		// mohan 2010-08-08
		//for(int ik=0; ik<kv.nks-1; ik++) delete[] ekb[ik];
		//delete[] ekb;
	}
}

void wavefunc::allocate_ekb_wg(const int nks)
{
    TITLE("wavefunc","init_local");
    this->npwx = this->setupIndGk(pw, nks);

	// band energies
	this->ekb = new double*[nks];
	for(int ik=0; ik<nks; ik++)
	{
		ekb[ik] = new double[NBANDS];
		ZEROS(ekb[ik],NBANDS);
	}
	this->allocate_ekb = true;

	// the weight of each k point and band
    this->wg.create(nks, NBANDS);
    Memory::record("wavefunc","ekb",nks*NBANDS,"double");
    Memory::record("wavefunc","wg",nks*NBANDS,"double");

    return;
}

void wavefunc::allocate(const int nks)
{	
	TITLE("wavefunc","allocate");

	this->npwx = this->setupIndGk(pw, nks);
	OUT(ofs_running,"npwx",npwx);

	assert(npwx > 0);
	assert(nks > 0);
	if( (CALCULATION!="scf-sto" && CALCULATION!="relax-sto" && CALCULATION!="md-sto") ) //qianrui add 
	assert(NBANDS > 0);

	// allocate for kinetic energy
	delete[] g2kin;
	this->g2kin = new double[npwx];
	ZEROS(g2kin, npwx);
	Memory::record("wavefunc","g2kin",npwx,"double");

	// if use spin orbital, do not double nks but double allocate evc and wanf2.
	int prefactor = 1;
	if(NSPIN==4) prefactor = NPOL;//added by zhengdy-soc
	
	this->ekb = new double*[nks];
	for(int ik=0; ik<nks; ik++)
	{
		ekb[ik] = new double[NBANDS];
		ZEROS(ekb[ik], NBANDS);
	}
	this->allocate_ekb = true;
	
	// the weight of each k point and band
	this->wg.create(nks, NBANDS);
	Memory::record("wavefunc","et",nks*NBANDS,"double");
	Memory::record("wavefunc","wg",nks*NBANDS,"double");

	delete[] evc;
	delete[] wanf2;

	const int nks2 = nks;

	if(CALCULATION=="nscf" && wf.mem_saver==1)
	{
		// mohan add 2010-09-07
		this->evc = new ComplexMatrix[1];
		this->wanf2 = new ComplexMatrix[1];
		
		// //added by zhengdy-soc
		evc[0].create(NBANDS, npwx * NPOL);

		if(BASIS_TYPE=="lcao_in_pw")
		{
			wanf2[0].create(NLOCAL, npwx * NPOL);
			cout << " Memory for wanf2 (MB): " << 
				Memory::record("wavefunc","wanf2",NLOCAL*(prefactor*npwx),"complexmatrix") << endl;
		}
		cout << " MEMORY FOR PSI (MB)  : " << 
			Memory::record("wavefunc","evc",NBANDS*(prefactor*npwx),"complexmatrix") << endl;
	}
	else
	{
		this->evc = new ComplexMatrix [nks2];
		this->wanf2 = new ComplexMatrix [nks2];

		for (int ik = 0; ik < nks2; ik++)
		{
			this->evc[ik].create(NBANDS, npwx * NPOL);//added by zhengdy-soc

			//Mohan add 2010-1-10
			if((BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2)
			{
				this->wanf2[ik].create(NLOCAL, npwx * NPOL);//added by zhengdy-soc
			}
		};

		cout << " MEMORY FOR PSI (MB)  : " << 
		Memory::record("wavefunc","evc",nks2*NBANDS*(prefactor*npwx),"complexmatrix") << endl;
	}

	//showMemStats();
	return;
}

//===================================================================
// This routine computes an estimate of the start_ wavefunctions
// from superposition of atomic wavefunctions or random wave functions.
//===================================================================
#include "occupy.h"
void wavefunc::wfcinit(void)
{
    TITLE("wavefunc","wfcinit");
    timer::tick("wavefunc","wfcinit",'C');

    this->wfcinit_k();

    en.demet = 0.0;

    //================================
    // Occupations are computed here
    //================================
	if(BASIS_TYPE=="pw")
	{
		// mohan fix bug 2011-02-25,
		// in nscf, occupations is not needed,
		if(CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax") //pengfei 2014-10-13
		{
    		Occupy::calculate_weights();
		}
	}
    if (test_wf>2)
    {
        out.printrm(ofs_running, " wg  ",  wg);
        this->check_psi(evc);
    }

    timer::tick("wavefunc","wfcinit",'C');
    return;
}

int wavefunc::get_starting_nw(void)const
{
    if (start_wfc == "file")
    {
		throw runtime_error("wavefunc::get_starting_nw. start_ wfc from file: not implemented yet! "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__)); 	// Peize Lin change 2019-05-01
        //WARNING_QUIT("wfcinit_k","\n start_ wfc from file: not implemented yet!");
        //**********************************************************************
        // ... read the wavefunction into memory (if it is not done in c_bands)
        //**********************************************************************
    }
    else if (start_wfc.substr(0,6) == "atomic")
    {
        if (ucell.natomwfc >= NBANDS)
        {
            if(test_wf)ofs_running << " Start wave functions are all pseudo atomic wave functions." << endl;
        }
        else
        {
            if(test_wf)ofs_running << " Start wave functions are atomic + "
            << NBANDS - ucell.natomwfc
            << " random wave functions." << endl;
        }
        return max(ucell.natomwfc,  NBANDS);
    }
    else if (start_wfc == "random")
    {
        if(test_wf)ofs_running << " Start wave functions are all random." << endl;
        return NBANDS;
    }
    else
    {
		throw runtime_error("wavefunc::get_starting_nw. Don't know what to do! Please Check source code! "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__)); 	// Peize Lin change 2019-05-01
        //WARNING_QUIT("get_starting_nw","Don't know what to do! Please Check source code!");
    }
}




void wavefunc::LCAO_in_pw_k(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","LCAO_in_pw_k");
	timer::tick("wavefunc","LCAO_in_pw_k",'G');

	assert(BASIS_TYPE=="lcao_in_pw");
	
	static bool ltable = false;

	if(!ltable)
	{
		this->table_local.create(ucell.ntype, ucell.nmax_total, NQX);
		Wavefunc_in_pw::make_table_q(ORB.orbital_file, this->table_local);
		ltable = true;
	}
	
	Wavefunc_in_pw::produce_local_basis_in_pw(ik, wvf, this->table_local);

	//-------------------------------------------------------------
	// (2) diago to get wf.ekb, then the weights can be calculated.
	//-------------------------------------------------------------
    hm.hpw.allocate(this->npwx, NPOL, ppcell.nkb, pw.nrxx);
	hm.hpw.init_k(ik);
	
	//hm.diagH_subspace(ik ,NLOCAL, NBANDS, wvf, wvf, ekb[ik]);
//	for(int ib=0; ib<NBANDS; ib++)
//	{
//		cout << " ib=" << ib << " e=" << ekb[ik][ib] << endl;
//	}

//	DONE(ofs_running,"CONSTRUCT_LOCAL_BASIS_IN_PW");

	timer::tick("wavefunc","LCAO_in_pw_k",'G');
	return;
}


void wavefunc::LCAO_in_pw_k_q(const int &ik, ComplexMatrix &wvf, Vector3<double> q)   // pengfei  2016-11-23
{
	TITLE("wavefunc","LCAO_in_pw_k_q");
	timer::tick("wavefunc","LCAO_in_pw_k_q",'G');
	//assert(LOCAL_BASIS==4); xiaohui modify 2013-09-01
	assert(BASIS_TYPE=="lcao_in_pw"); //xiaohui add 2013-09-01. Attention! How about "BASIS_TYPE=="lcao""???

	Wavefunc_in_pw::produce_local_basis_q_in_pw(ik, wvf, this->table_local, q);

	timer::tick("wavefunc","LCAO_in_pw_k_q",'G');
	return;
}


void wavefunc::diago_PAO_in_pw_k(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","diago_PAO_in_pw_k");

	hm.hpw.init_k(ik);
    this->diago_PAO_in_pw_k2(ik, wvf);

	return;
}

void wavefunc::diago_PAO_in_pw_k2(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","diago_PAO_in_pw_k2");
	// (6) Prepare for atmoic orbitals or random orbitals
	const int starting_nw = this->get_starting_nw();
	assert(starting_nw > 0);

	ComplexMatrix wfcatom(starting_nw, npwx * NPOL);//added by zhengdy-soc
	if(test_wf)OUT(ofs_running, "starting_nw", starting_nw);
	if(start_wfc.substr(0,6)=="atomic")
	{
		this->atomic_wfc(ik, this->npw, ucell.lmax_ppwf, wfcatom, ppcell.tab_at, NQX, DQ);
		if( start_wfc == "atomic+random" && starting_nw == ucell.natomwfc )//added by qianrui 2021-5-16
		{
			double rr, arg;
			for(int ib = 0 ; ib < starting_nw ; ++ib )
			{
				int startig = 0;
				for(int ip = 0 ; ip < NPOL; ++ip)
				{
					for(int ig = 0 ; ig < npw ; ++ig)
					{
						rr = rand()/double(RAND_MAX);
						arg = TWO_PI * rand()/double(RAND_MAX);
						wfcatom(ib,startig+ig) *= (1.0 + 0.05 * complex<double>(rr * cos(arg), rr * sin(arg)));
					}
					startig += npwx;
				}
			}
		}
		
		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		this->random(wfcatom, ucell.natomwfc, NBANDS, ik);
	}
	else if(start_wfc=="random")
	{
		this->random(wfcatom,0,NBANDS,ik);
	}

	// (7) Diago with cg method.
	double *etatom  = new double[starting_nw];
	ZEROS(etatom, starting_nw);
	//if(DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
	if(KS_SOLVER=="cg") //xiaohui add 2013-09-02
	{
		hm.diagH_subspace(ik ,starting_nw, NBANDS, wfcatom, wfcatom, etatom);
	}

	/*
	ofs_running << " " << "ik = " << ik << " Bands(eV)" << endl;
	for (int ib=0;ib<starting_nw;ib++)
	{
		ofs_running << " " << setw(15) << etatom[ib]*Ry_to_eV << endl;
	}
	*/

	assert(wvf.nr <= wfcatom.nr);
	for (int ib=0; ib<NBANDS; ib++)
	{
		for (int ig=0; ig<this->npwx; ig++)
		{
			wvf(ib, ig) = wfcatom(ib, ig);
			if(NPOL==2) wvf(ib,ig + this->npwx) = wfcatom(ib,ig + this->npwx);
		}
	}

	//added by zhengdy-soc
/*	for(int i = 0;i<starting_nw;i++)
	{
		ekb[ik][i] = etatom[i];
	}*/
	delete[] etatom;
}

void wavefunc::wfcinit_k(void)
{
	TITLE("wavefunc","wfcinit_k");

	if(mem_saver) 
	{
		return;
	}

	for(int ik=0; ik<kv.nks; ik++)
	{
		if(BASIS_TYPE=="lcao_in_pw")
		{
			// just get the numerical local basis wave functions
			// in plane wave basis
			this->LCAO_in_pw_k(ik, wf.wanf2[ik]);
		}
		else if (BASIS_TYPE=="pw")
		{
			// get the wave functions 
			// by first diagolize PAO
			// wave functions.
			this->diago_PAO_in_pw_k(ik, wf.evc[ik]);
		}
	}
	
	//---------------------------------------------------
	//  calculte the overlap <i,0 | e^{i(q+G)r} | j,R>
	//---------------------------------------------------
	
	if((!chi0_hilbert.epsilon) && chi0_hilbert.kmesh_interpolation )    // pengfei  2016-11-23
	{
		chi0_hilbert.Parallel_G();    // for parallel: make sure in each core, G(all_gcars(pw.gcars))  are the same
		
		// iw1->i, iw2->j, R store the positions of the neighbor unitcells that |i,0> and |j,R> have overlaps
		R = new Vector3<int>** [NLOCAL];
		for(int iw1=0; iw1<NLOCAL; iw1++)
		{
			R[iw1] = new Vector3<int>* [NLOCAL];
			for(int iw2=0; iw2<NLOCAL; iw2++)
			{
				R[iw1][iw2] = new Vector3<int>[chi0_hilbert.lcao_box[0]*chi0_hilbert.lcao_box[1]*chi0_hilbert.lcao_box[2]];
			}
		}
		
		Rmax = new int* [NLOCAL];
		for(int iw=0; iw<NLOCAL; iw++)
		{
			Rmax[iw] = new int[NLOCAL];
		}
		
		int NR; // The Max number of the overlaps for each iw1,iw2; 
		NR = get_R(chi0_hilbert.lcao_box[0],chi0_hilbert.lcao_box[1],chi0_hilbert.lcao_box[2]); 
		
		// store the overlap relationship to "nearest.dat"
		stringstream ss;
		ss << global_out_dir <<"nearest.dat";
		ofstream ofs(ss.str().c_str());
		ofs << NR << endl;
		cout <<"NR = "<<NR<<endl;    // Max
		for(int iw1=0; iw1<NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<NLOCAL; iw2++)
			{
				ofs<<iw1<<"   "<<iw2<<"    "<<Rmax[iw1][iw2]<<endl;   // iw1, iw2, and how many overlaps between iw1 and iw2
				for(int i=0; i<Rmax[iw1][iw2]; i++)
				{
					ofs<<R[iw1][iw2][i].x<<"  "<<R[iw1][iw2][i].y<<"  "<<R[iw1][iw2][i].z<<endl;   // positions
				}
			}
		}
		ofs.close();
		
		int NG = chi0_hilbert.dim;  // chi0's dimension
		
		complex<double> ***wanf2_q;    // <j,0 | k+G+q>
		
		wanf2_q = new complex<double> **[kv.nks];
		for(int ik=0; ik<kv.nks; ik++)
		{
			wanf2_q[ik] = new complex<double> *[NLOCAL];
			for(int iw=0; iw<NLOCAL; iw++)
			{
				wanf2_q[ik][iw] = new complex<double>[npwx];
			}
		}
		
		complex<double> overlap_aux[NLOCAL][NLOCAL][NG][NR];     // <i,0 | e^{i(q+G)r} | j,R> 
		complex<double> overlap[NLOCAL][NLOCAL][NG][NR];
		double overlap_aux_R[NLOCAL][NLOCAL][NG][NR];       //real part
		double overlap_R[NLOCAL][NLOCAL][NG][NR];
		double overlap_aux_I[NLOCAL][NLOCAL][NG][NR];       //imag part
		double overlap_I[NLOCAL][NLOCAL][NG][NR];
		
		ComplexMatrix Mat;
		Mat.create(NLOCAL,npwx);
		Vector3<double> qg; // q+G
		Vector3<double> gkqg, Rcar[NLOCAL][NLOCAL][NR];  // k+G+qg, Rcartesian
		
		for(int iw1=0;iw1<NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<NLOCAL; iw2++)
			{
				for(int i=0; i<NR; i++)
				{
					Rcar[iw1][iw2][i].x = 0.0; Rcar[iw1][iw2][i].y = 0.0;Rcar[iw1][iw2][i].z = 0.0;
				}
			}
		}

		for(int iw1=0; iw1<NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<NLOCAL; iw2++)
			{
				for(int i=0; i<Rmax[iw1][iw2]; i++)
				{
					Rcar[iw1][iw2][i].x = ucell.latvec.e11 * R[iw1][iw2][i].x 
					+ ucell.latvec.e21 * R[iw1][iw2][i].y + ucell.latvec.e31 * R[iw1][iw2][i].z;
					Rcar[iw1][iw2][i].y = ucell.latvec.e12 * R[iw1][iw2][i].x 
					+ ucell.latvec.e22 * R[iw1][iw2][i].y + ucell.latvec.e32 * R[iw1][iw2][i].z;
					Rcar[iw1][iw2][i].z = ucell.latvec.e13 * R[iw1][iw2][i].x 
					+ ucell.latvec.e23 * R[iw1][iw2][i].y + ucell.latvec.e33 * R[iw1][iw2][i].z;
				}
			}
		}
 
		
		double arg; 
		complex<double> phase;
		for(int iq=0; iq<chi0_hilbert.nq; iq++)
		{
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					for(int ig=0; ig<NG; ig++)
					{
						ZEROS( overlap_aux[iw1][iw2][ig], NR);
					}
				}
			}
		}
		
		// main
		for(int iq=0; iq<chi0_hilbert.nq; iq++)    // loop over iq 
		{
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					for(int ig=0; ig<NG; ig++)
					{
						ZEROS( overlap_aux[iw1][iw2][ig], NR);
					}
				}
			}
			
			for(int g=0; g<NG; g++)    // loop over ig
			{
				qg.x = chi0_hilbert.qcar[iq][0] + chi0_hilbert.all_gcar[g].x;   
				qg.y = chi0_hilbert.qcar[iq][1] + chi0_hilbert.all_gcar[g].y; 
				qg.z = chi0_hilbert.qcar[iq][2] + chi0_hilbert.all_gcar[g].z;
				cout <<"qg = "<<qg.x<<" "<<qg.y<<" "<<qg.z<<endl;
				for(int ik=0; ik<kv.nks; ik++)
				{
					this->LCAO_in_pw_k_q(ik, Mat, qg);
					for(int iw=0; iw<NLOCAL; iw++)
					{
						for(int ig=0; ig<kv.ngk[ik]; ig++)
						{
							wanf2_q[ik][iw][ig] = Mat(iw,ig);
						}										
					}
				}
				for(int iw1=0; iw1<NLOCAL; iw1++)    // loop over iw1
				{
					for(int iw2=0; iw2<NLOCAL; iw2++)  // loop over iw2
					{
						for(int ik=0;ik<kv.nks;ik++)         // loop over ik
						{
							for(int ig=0;ig<kv.ngk[ik];ig++)    // loop over ig
							{
								gkqg = pw.get_GPlusK_cartesian(ik, wf.igk(ik, ig)) + qg;
								for(int ir=0; ir<Rmax[iw1][iw2]; ir++)   // Rmax
								{
									arg = gkqg * Rcar[iw1][iw2][ir] * TWO_PI;
									phase = complex<double>( cos(arg),  -sin(arg) );
									overlap_aux[iw1][iw2][g][ir] += conj(wf.wanf2[ik](iw1,ig)) 
									* wanf2_q[ik][iw2][ig] * phase/static_cast<double>(kv.nks);		
									// Peize Lin add static_cast 2018-07-14
								}
							}
						}	
					}
				}
			}
			
			for(int g=0; g<NG; g++)
			{
				for(int iw1=0; iw1<NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<NLOCAL; iw2++)
					{
						for(int ir=0; ir<NR; ir++)
						{
							overlap_aux_R[iw1][iw2][g][ir] = overlap_aux[iw1][iw2][g][ir].real();
							overlap_aux_I[iw1][iw2][g][ir] = overlap_aux[iw1][iw2][g][ir].imag();
						}
					}
				}
			}
			
#ifdef __MPI
			MPI_Allreduce(overlap_aux_R,overlap_R,NLOCAL * NLOCAL * NG * NR,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
			MPI_Allreduce(overlap_aux_I,overlap_I,NLOCAL * NLOCAL * NG * NR,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif
			
			for(int g=0; g<NG; g++)
			{
				for(int iw1=0; iw1<NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<NLOCAL; iw2++)
					{
						for(int ir=0; ir<NR; ir++)
						{
							overlap[iw1][iw2][g][ir] = complex<double>( overlap_R[iw1][iw2][g][ir], 
							overlap_I[iw1][iw2][g][ir]);
						}
					}
				}
			}
						
			//------------------------------
			// store the overlap in q_(iq)
			//------------------------------
			stringstream ss1;
			ss1 << global_out_dir <<"q_"<<iq;
			ofstream ofs1(ss1.str().c_str());
			ofs1<<NG<<endl;
			for(int g=0; g<NG; g++)
			{
				for(int iw1=0; iw1<NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<NLOCAL; iw2++)
					{
						for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
						{
							ofs1<<overlap[iw1][iw2][g][ir]<<"  ";
						}
						ofs1<<endl;
					}
				}
				// mohan update 2021-02-24
				ofs1<<endl; ofs1<<endl;
			}
			ofs1.close();
		}
		
		for(int iw1=0; iw1<NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<NLOCAL; iw2++)
			{
				delete[] R[iw1][iw2];
			}
			delete[] R[iw1];
		}
		delete[] R;
		
		for(int iw=0; iw<NLOCAL; iw++)
		{
			delete[] Rmax[iw];
		}
		delete[] Rmax;
		
		for(int ik=0; ik<kv.nks; ik++)
		{
			for(int iw=0; iw<NLOCAL; iw++)
			{
				delete[] wanf2_q[ik][iw];
			}
			delete[] wanf2_q[ik];
		}
		delete[] wanf2_q;
	
	}
	
	return;
}

//--------------------------------------------
// get the nearest unitcell positions 
// that exist overlaps between two orbitals
// iw1 and iw2
//--------------------------------------------
int wavefunc::get_R(int ix, int iy, int iz)   // pengfei 2016-11-23
{
	int count;
	Vector3<double> r,r1,r2;

	for(int iw1=0; iw1<NLOCAL; iw1++)
	{
		for(int iw2=0; iw2<NLOCAL; iw2++)
		{
			int it1 = iw2it(iw1); int ia1 = iw2ia(iw1);
			int it2 = iw2it(iw2); int ia2 = iw2ia(iw2);
			//cout <<"iw1= "<<iw1<<" iw2= "<<iw2<<" it1= "<<it1<<" ia1= "<<ia1<<" it2= "<<it2<<" ia2= "<<ia2<<endl;
			count = 0;

			for(int nx=-int(ix/2);nx<=int(ix/2);nx++)
			{
				for(int ny=-int(iy/2);ny<=int(iy/2);ny++)
				{
					for(int nz=-int(iz/2);nz<=int(iz/2);nz++)
					{
						//cout <<"count = "<<count<<endl;
						//cout<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<endl;
						r1.x = ucell.atoms[it1].tau[ia1].x * ucell.lat0;
						r1.y = ucell.atoms[it1].tau[ia1].y * ucell.lat0;
						r1.z = ucell.atoms[it1].tau[ia1].z * ucell.lat0;
						r2.x = (ucell.atoms[it2].tau[ia2].x 
						+ ucell.latvec.e11 * nx + ucell.latvec.e21 * ny + ucell.latvec.e31 * nz) * ucell.lat0;
						r2.y = (ucell.atoms[it2].tau[ia2].y 
						+ ucell.latvec.e12 * nx + ucell.latvec.e22 * ny + ucell.latvec.e32 * nz) * ucell.lat0;
						r2.z = (ucell.atoms[it2].tau[ia2].z 
						+ ucell.latvec.e13 * nx + ucell.latvec.e23 * ny + ucell.latvec.e33 * nz) * ucell.lat0;
						r = r2 - r1;
						double distance = sqrt(r*r);

						if(distance < (ucell.atoms[it1].Rcut + ucell.atoms[it2].Rcut))
						{
							R[iw1][iw2][count].x = nx; 
							R[iw1][iw2][count].y = ny; 
							R[iw1][iw2][count].z = nz; 
							count++;
						}
					}
				}
			}
			Rmax[iw1][iw2] = count;
		}
	}

	int NR = 0;
	for(int iw1=0;iw1<NLOCAL;iw1++)
	{
		for(int iw2=0;iw2<NLOCAL;iw2++)
		{
			if(Rmax[iw1][iw2] > NR)
			{
				NR = Rmax[iw1][iw2];
			}
		}
	}

	return NR;
}


int wavefunc::iw2it( int iw)    // pengfei 2016-11-23
{
    int ic, type;
    ic =0;
    for(int it =0; it<ucell.ntype; it++)
	{
        for(int ia = 0; ia<ucell.atoms[it].na; ia++)
        {
            for(int L=0; L<ucell.atoms[it].nwl+1; L++)
			{
                for(int N=0; N<ucell.atoms[it].l_nchi[L]; N++)
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

int wavefunc::iw2ia( int iw)    // pengfei 2016-11-23
{
	int ic, na;
	ic =0;
	for(int it =0; it<ucell.ntype; it++)
	{
		for(int ia = 0; ia<ucell.atoms[it].na; ia++)
		{
			for(int L=0; L<ucell.atoms[it].nwl+1; L++)
				for(int N=0; N<ucell.atoms[it].l_nchi[L]; N++)
				{
					for(int i=0; i<(2*L+1); i++)
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
	return na;
}

//LiuXh add a new function here,
//20180515
void wavefunc::init_after_vc(const int nks)
{
    static bool done_once = false;
    if(done_once)
    {
        //return; //LiuXh add 2017-12-12
    }
    else
    {
        done_once = true;
    }

    TITLE("wavefunc","init");
    //this->npwx = this->setupIndGk(pw, nks);
    OUT(ofs_running,"npwx",npwx);

    assert(npwx > 0);
    assert(nks > 0);
    assert(NBANDS > 0);

    delete[] g2kin;
    this->g2kin = new double[npwx];   // [npw],kinetic energy
    ZEROS(g2kin, npwx);
    Memory::record("wavefunc","g2kin",npwx,"double");
    if(test_wf)OUT(ofs_running,"g2kin allocation","Done");

    int prefactor = 1;
    this->ekb = new double*[nks];
    for(int ik=0; ik<nks; ik++)
    {
        ekb[ik] = new double[NBANDS];
        ZEROS(ekb[ik], NBANDS);
    }
    this->allocate_ekb = true;

    this->wg.create(nks, NBANDS);       // the weight of each k point and band
    Memory::record("wavefunc","et",nks*NBANDS,"double");
    Memory::record("wavefunc","wg",nks*NBANDS,"double");
    if(test_wf)OUT(ofs_running, "et allocation","Done");
    if(test_wf)OUT(ofs_running, "wg allocation","Done");

    delete[] evc;
    delete[] wanf2;

    const int nks2 = nks * prefactor;

    if(CALCULATION=="nscf" && wf.mem_saver==1)
    {
        this->evc = new ComplexMatrix[1];
        this->wanf2 = new ComplexMatrix[1];

        evc[0].create(NBANDS*prefactor, npwx);
        if(BASIS_TYPE=="lcao_in_pw")
        {
            wanf2[0].create(NLOCAL*prefactor, npwx);
            cout << " Memory for wanf2 (MB): " <<
            Memory::record("wavefunc","wanf2",(NLOCAL*prefactor)*npwx,"complexmatrix") << endl;
        }
        cout << " MEMORY FOR PSI (MB)  : " <<
        Memory::record("wavefunc","evc",(NBANDS*prefactor)*npwx,"complexmatrix") << endl;
    }
    else
    {
        this->evc = new ComplexMatrix [nks2];
        this->wanf2 = new ComplexMatrix [nks2];

        for (int ik = 0; ik < nks2; ik++)
        {
            this->evc[ik].create(NBANDS*prefactor, npwx);

            if((BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2)
            {
                this->wanf2[ik].create(NLOCAL, npwx);
            }
        }

        cout << " MEMORY FOR PSI (MB)  : " <<
        Memory::record("wavefunc","evc",(nks*prefactor)*(NBANDS*prefactor)*npwx,"complexmatrix") << endl;
    }

    if(test_wf)
    {
        OUT(ofs_running,"evc allocation","Done");
        if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw")
        {
            OUT(ofs_running,"wanf2 allocation","Done");
        }
    }

    return;
}

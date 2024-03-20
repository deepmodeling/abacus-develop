#include "ReadData.h"
#include "../src_parallel/parallel_common.h"

ReadData::ReadData()
: test(0)
{
	weight = new double[1];
}

ReadData::~ReadData()
{
	if(TEST1) cout << "\n ~ReadData()" << endl;
	delete[] weight;
}

//==========================================================
// Readin Overlap Q(k) = < J(k,ie) | Bloch(k) >
// The number of local(ie) orbitals is : nwfc;
// The number of Bloch orbitals is : nbands;
// The number of ie is : ne
// The number of k points : nks;
//==========================================================
void ReadData::OverlapQandS(const string &name)
{
	TITLE(ofs_running, "ReadData", "OverlapQandS");
	
	ifstream ifs;
	
	double* weighttmp = new double[NKSTOT];
	double* carkx = new double[NKSTOT];
	double* carky = new double[NKSTOT];
	double* carkz = new double[NKSTOT];
	ZEROS(weighttmp, NKSTOT);
	
	if(MY_RANK==0)
	{
		ifs.open(name.c_str());

		if (!ifs)
		{
			cout << "\n Can't find file : "  << name;
			WARNING_QUIT("ReadData::OverlapQandS","Can't find file.");	
		}
		else
		{
			cout << " FILE : " << name << endl;
		}

		double lat0;
		ifs >> lat0;
		assert(lat0 == LAT0);

		double e11,e12,e13,e21,e22,e23,e31,e32,e33;

		ifs >> e11 >> e12 >> e13;
		ifs >> e21 >> e22 >> e23;
		ifs >> e31 >> e32 >> e33;

		assert( e11 == LATVEC.e11);
		assert( e12 == LATVEC.e12);
		assert( e13 == LATVEC.e13);
		assert( e21 == LATVEC.e21);
		assert( e22 == LATVEC.e22);
		assert( e23 == LATVEC.e23);
		assert( e31 == LATVEC.e31);
		assert( e32 == LATVEC.e32);
		assert( e33 == LATVEC.e33);

		int ntype;
		READ_VALUE(ifs, ntype); // 1
		assert(ntype==NTYPE);
	
		double x,y,z;
		for(int it=0; it<ntype; it++) 
		{
			string tmp_label;
			READ_VALUE(ifs, tmp_label); // 2
			assert(tmp_label == LABEL[it]);
			
			int na;
			READ_VALUE(ifs, na);// 3 
			assert(na == NA[it]);
			for(int ia=0; ia<NA[it]; ia++)
			{
				//mohan fix bug 2010-09-27
				// x,y,z is useless.
				// And CARPOSX, CASPOSY, CASPOSZ
				// only contains the first structure.
				ifs >> x >> y >> z;
			}
		}

		double tmp_ecutwfc;
		double tmp_ecutwfc_jlq;
		double tmp_rcut;

		READ_VALUE(ifs, tmp_ecutwfc);// 4
		READ_VALUE(ifs, tmp_ecutwfc_jlq);// 5
		READ_VALUE(ifs, tmp_rcut);// 6

		assert(tmp_ecutwfc == ECUT);
		assert(tmp_ecutwfc_jlq == ECUT_JLQ);
		assert(tmp_rcut == RCUT);

		// mohan add 2009-08-28
		// make sure all structures have 
		// same SMOOTH and SIGMA value!
		bool tmp_smooth;
		double tmp_sigma;
		double tmp_tolerence;
	
		READ_VALUE(ifs, tmp_smooth); // 7
		READ_VALUE(ifs, tmp_sigma); // 8
		READ_VALUE(ifs, tmp_tolerence); // 9

		assert(tmp_smooth == SMOOTH);
		assert(tmp_sigma == SIGMA);
		assert(tmp_tolerence == TOLERENCE);

		int lmaxall, nkstot, nbands, nwfcall, ne;

		READ_VALUE(ifs, lmaxall); // 10
		READ_VALUE(ifs, nkstot);// 11
		READ_VALUE(ifs, nbands);// 12
		READ_VALUE(ifs, nwfcall);// 13
		READ_VALUE(ifs, ne);// 14

		assert(lmaxall == LMAXALL);
		assert(nkstot == NKSTOT);
		assert(nbands == NBANDS);
		assert(nwfcall == NWFCALL);
		assert(ne == NE);

		//==========================================
		// if use IBZ k-points, weight is different
		// for each k point
		//==========================================
		// 0 means don't need to search from start.
		if( SCAN_BEGIN(ifs,"<WEIGHT_OF_KPOINTS>",0) )
		{
			ofs_running << " Kx/Ky/Kz/Kweight" ;
			for(int ik=0; ik<NKSTOT; ik++) 
			{
				if(ik%4==0) ofs_running << endl;
				ifs >> carkx[ik] 
				>> carky[ik] 
				>> carkz[ik]
				>> weighttmp[ik];
				ofs_running << " " << carkx[ik]
				<< " " << carky[ik]
				<< " " << carkz[ik] 
				<< " " << weighttmp[ik] << endl;
			}
			SCAN_END(ifs,"</WEIGHT_OF_KPOINTS>");
		}

		double sum = 0;
		for(int ik=0; ik<NKSTOT; ik++) sum+= weighttmp[ik];

		//ofs_running << "\n sum kpoint weight = " << sum;
		if( abs(sum-1.0) > 1.0e-5 )
		{
			ofs_running << "\n sum of k weight = " << sum << endl;
			WARNING_QUIT("ReadData::OverlapQandS","sum of k weight is wrong.");
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_double(weighttmp, NKSTOT);
	Parallel_Common::bcast_double(carkx, NKSTOT);
	Parallel_Common::bcast_double(carky, NKSTOT);
	Parallel_Common::bcast_double(carkz, NKSTOT);
#endif

	ofs_running << " NKS=" << NKS << endl;


	delete[] this->weight;
	this->weight = new double[NKS];
	ZEROS(weight, NKS);

	CARKX = new double[NKS];
	CARKY = new double[NKS];
	CARKZ = new double[NKS];
	ZEROS(CARKX, NKS);
	ZEROS(CARKY, NKS);
	ZEROS(CARKZ, NKS);


	for(int ik=0; ik<NKSTOT; ik++)
	{
#ifdef __MPI
		const int iknow = ik - Pkpoints.startk_pool[MY_POOL];
		const int pool = Pkpoints.whichpool[ik];
		if(MY_POOL == pool)
		{
			this->weight[iknow] = weighttmp[ik];
			CARKX[iknow] = carkx[ik];
			CARKY[iknow] = carky[ik];
			CARKZ[iknow] = carkz[ik];
			//ofs_running << "\n iknow = " << iknow << " ik = " << ik << " weight = " << weight[iknow];
		}
#else
		weight[ik] = weighttmp[ik];
		CARKX[ik] = carkx[ik];
		CARKY[ik] = carky[ik];
		CARKZ[ik] = carkz[ik];
#endif
	}


	delete[] weighttmp;
	delete[] carkx;
	delete[] carky;
	delete[] carkz;

	

	this->OverlapQ(ifs);



	this->OverlapSq1q2(ifs);



	return;
}

void ReadData::OverlapSq1q2(ifstream &ifs)
{
//	TITLE("ReadData", "OverlapSq1q2");
	timer::tick("ReadData", "OverlapSq1q2");

	USEPW = true;
	if(MY_RANK==0)
	{
		bool restart = false;
		bool quit = false;
		if( SCAN_BEGIN(ifs,"<OVERLAP_Sq>",restart,quit) )
		{
			USEPW = false;
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_bool(USEPW);
#endif
	
	if(USEPW)
	{
		ofs_running << " USE PLANE WAVE BASIS" << endl;
	}
	else
	{
		ofs_running << " USE Q and S matrix" << endl;
	}

	if(CALSPI==0) return;

	if(!USEPW)
	{
		this->Sq1q2 = new ComplexArray[ NKS ];
		for(int ik=0; ik<NKS; ik++)
		{
			this->Sq1q2[ik].create(NWFCALL, NWFCALL, NE, NE);
		}
	}


#ifdef __MPI
	if(!USEPW)
	{
		const int ndata = NWFCALL * NWFCALL * NE * NE;
		double* buffera = new double[ndata];
		double* bufferb = new double[ndata];
		for(int ik=0; ik<NKSTOT; ik++)
		{
			ZEROS(buffera, ndata);
			ZEROS(bufferb, ndata);
			if(MY_RANK==0)
			{
				for(int id=0; id<ndata; id++)
				{
					ifs >> buffera[id] >> bufferb[id];
				}
			}
			//ofs_running << "\n ik=" << ik;
			MPI_Status ierror;
			const int pool = Pkpoints.whichpool[ik];
			const int iknow = ik - Pkpoints.startk_pool[MY_POOL];
			const int startdata = 0; // this is differnt from Q.
			
			//ofs_running << " this k point belong to pool : " << pool;
			//ofs_running << " my pool : " << MY_POOL << endl;

			if(MY_RANK==0)
			{
				// always ready to send data to other processors.
				if(pool==0)
				{
					// send copies (nproc - 1).
					for(int ip=1; ip<Pkpoints.nproc_pool[pool]; ip++)
					{
						MPI_Send(buffera, ndata, MPI_DOUBLE, ip, ik, MPI_COMM_WORLD);
						MPI_Send(bufferb, ndata, MPI_DOUBLE, ip, ik, MPI_COMM_WORLD);
			//			ofs_running << "\n send data to processor " << ip;
					}

					// local copy
					for(int id=0; id<ndata; id++)
					{
						this->Sq1q2[iknow].ptr[id + startdata] = complex<double> (buffera[id], bufferb[id]);
					}
			//		ofs_running << "\n local data copy ";
				}
				else
				{
					// send copys: nproc
					for(int ip=0; ip<Pkpoints.nproc_pool[pool]; ip++)
					{
						const int ipall = ip + Pkpoints.startpro_pool[pool];
						MPI_Send(buffera, ndata, MPI_DOUBLE, ipall, ik, MPI_COMM_WORLD);
						MPI_Send(bufferb, ndata, MPI_DOUBLE, ipall, ik, MPI_COMM_WORLD);
			//			ofs_running << "\n send data to processor " << ipall;
					}
				}
			}//end MY_RANK==0
			else
			{
				if(pool == MY_POOL)
				{
			//		ofs_running << "\n begin receive, ndata : " << ndata << endl;
					MPI_Recv(buffera, ndata, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD, &ierror);
					MPI_Recv(bufferb, ndata, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD, &ierror);
					
					for(int id=0; id<ndata; id++)
					{
						this->Sq1q2[iknow].ptr[id + startdata] = complex<double> (buffera[id], bufferb[id]);
					}
			//		ofs_running << "\n receive data end" << endl;
				}
				else
				{
				//	ofs_running << "\n do nothing";
				}
			}		
			MPI_Barrier(MPI_COMM_WORLD);
		}
		delete[] buffera;
		delete[] bufferb;
	}// end begin
#else
	assert(NKS==NKSTOT);
	if(!USEPW)
	{
		double a,b;
		for (int ik = 0; ik < NKSTOT; ik++)
		{
			this->Sq1q2[ik].create(NWFCALL, NWFCALL, NE, NE);
	
			for (int i = 0; i < Sq1q2[ik].getSize(); i++)
			{
				ifs >> a >> b;
				Sq1q2[ik].ptr[i] = complex<double>(a,b);
			}
		}
	}
#endif

//	SCAN_END(ifs,"</OVERLAP_Sq>");

	timer::tick("ReadData", "OverlapSq1q2");
	return;
}

void ReadData::OverlapQ(ifstream &ifs)
{
	if(CALSPI==0) return;
//	TITLE("ReadData", "OverlapQ");
	timer::tick("ReadData", "OverlapQ");

	assert(NKS > 0);
	assert(NBANDS > 0);
	assert(NWFCALL > 0);
	assert(NE > 0);

	this->Qin.create(NKS, NBANDS, NWFCALL, NE);

	ofs_running << " DIMENSION OF Q=<Psi|Jlq>" << endl;
	OUT(ofs_running,"NKS",NKS);
	OUT(ofs_running,"NBANDS",NBANDS);
	OUT(ofs_running,"NWFCALL",NWFCALL);
	OUT(ofs_running,"NE",NE);

	bool begin = false;
	
	if(MY_RANK ==0 )
	{
		if( SCAN_BEGIN(ifs,"<OVERLAP_Q>",0) )
		{
			begin = true;
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_bool(begin);
#endif



	if(begin)
	{
		double a,b;
#ifdef __MPI
		const int ndata = NBANDS * NWFCALL * NE;
		double* buffera = new double[ndata];
		double* bufferb = new double[ndata];
		for(int ik=0; ik<NKSTOT; ik++)
		{
			ZEROS(buffera, ndata);
			ZEROS(bufferb, ndata);
			if(MY_RANK==0)
			{
				for(int id=0; id<ndata; id++)
				{
					ifs >> buffera[id] >> bufferb[id];
				}
			}

			//ofs_running << "\n ik=" << ik;
			
			MPI_Status ierror;
			const int pool = Pkpoints.whichpool[ik];
			const int iknow = ik - Pkpoints.startk_pool[MY_POOL];
			const int startdata = iknow * ndata;

			//ofs_running << "\n this k point belong to pool : " << pool;
			//ofs_running << "\n my pool : " << MY_POOL;

			if(MY_RANK==0)
			{
				// always ready to send data to other processors.
				if(pool==0)
				{
					// send copies (nproc - 1).
					for(int ip=1; ip<Pkpoints.nproc_pool[pool]; ip++)
					{
						MPI_Send(buffera, ndata, MPI_DOUBLE, ip, ik, MPI_COMM_WORLD);
						MPI_Send(bufferb, ndata, MPI_DOUBLE, ip, ik, MPI_COMM_WORLD);
			//			ofs_running << "\n send data to processor " << ip;
					}
				
					// local copy
					for(int id=0; id<ndata; id++)
					{
						this->Qin.ptr[id + startdata] = complex<double> (buffera[id], bufferb[id]);
					}
					//ofs_running << "\n local data copy ";
				}
				else
				{
					// send copys: nproc
					for(int ip=0; ip<Pkpoints.nproc_pool[pool]; ip++)
					{
						const int ipall = ip + Pkpoints.startpro_pool[pool];
						MPI_Send(buffera, ndata, MPI_DOUBLE, ipall, ik, MPI_COMM_WORLD);
						MPI_Send(bufferb, ndata, MPI_DOUBLE, ipall, ik, MPI_COMM_WORLD);
					//	ofs_running << "\n send data to processor " << ipall;
					}
				}	
			}// end MY_RANK==0
			else
			{
				if(pool == MY_POOL)
				{
					//ofs_running << "\n begin receive, ndata : " << ndata << endl;
					MPI_Recv(buffera, ndata, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD, &ierror);
					MPI_Recv(bufferb, ndata, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD, &ierror);
					for(int id=0; id<ndata; id++)
					{
						this->Qin.ptr[id + startdata] = complex<double> (buffera[id], bufferb[id]);
					}
					//ofs_running << "\n receive data end" << endl;
				}
				else
				{
					//ofs_running << "\n do nothing";
				}
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
		}// endik
		delete[] buffera;
		delete[] bufferb;
#else
		for(int i=0; i<Qin.getSize(); i++)
		{
			ifs >> a >> b;
			this->Qin.ptr[i] = complex<double>( a, b );
		}
#endif
	}



	//SCAN_END(ifs,"</OVERLAP_Q>");
	
	timer::tick("ReadData", "OverlapQ");

	return;
}

//==========================================================
// Readin inverse S matrix.
// Sinv is fail to use if C4 change.
//==========================================================
void ReadData::OverlapSinv(const string &name)
{
//	TITLE("ReadData", "OverlapSinv");
	ifstream ifs(name.c_str());

	if (!ifs)
	{
		cout << "\n Can't find file : "  << name;
		cout << "\n But maybe we don't need it.";
		return;
	}

	assert(NKS > 0);
	assert(NWFCALL > 0);

	this->Sinv.create(NKS, NWFCALL, NWFCALL);

	return;

	string tmp;
	int row, col;

	for (int i = 0; i < NKS; i++)
	{
		ifs >> tmp >> row >> col;
//		cout << "\n" << tmp << " row=" << row << " col=" << col << endl;
		assert(row == NWFCALL);
		assert(col == NWFCALL);

		for (int j = 0; j < row; j++)
		{
			for (int k = 0; k < col; k++)
			{
				double a,b;
				ifs >> a >> b;
				Sinv(i,j,k) = complex<double>(a,b);
			}
		}
	}
	ifs.close();
	return;
}

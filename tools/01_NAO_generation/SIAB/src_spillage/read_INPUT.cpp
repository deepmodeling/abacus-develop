#include "read_INPUT.h"
#include "../src_parallel/parallel_common.h"

Read_INPUT::Read_INPUT()
{
    QS_data = new ReadData[1];
    test = 0;
    qsfile = new string[1];
}

Read_INPUT::~Read_INPUT()
{
    if (TEST1) cout << "\n ~Read_INPUT()" << endl;
    delete[] QS_data;
    delete[] qsfile;
}

void Read_INPUT::init(void)
{
    if (test==1) TITLE("Read_INPUT","init");

	stringstream ss;
	ss << "running_" << MY_RANK+1 << ".txt";
	ofs_running.open( ss.str().c_str() );

    if (MY_RANK==0)
    {
        ifs.open("INPUT");
		ifs2.open("INPUTs");
       
	   	if( ifs && ifs2)
		{
			cout << " \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cout << " Confused with the coexsitance of INPUT and INPUTs files." << endl;
			cout << " Please keep only one." << endl;
			cout << " INPUT (used to generate local orbitals, need overlap data)." << endl;
			cout << " INPUTs (used to generat overlap data)." << endl;
			cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            WARNING_QUIT("Read_INPUT","Find both files: INPUT and INPUTs");
		}
		if( ifs && !ifs2)
		{
			cout << " \n ========================================================" << endl;
			cout << " Find file: INPUT (used to generate local orbitals)." << endl;
			cout << " Can't find file: INPUTs (used to generate overlap data)" << endl;
			cout << " Minimize the spillage now." << endl;
			cout << " ========================================================" << endl;
		}
		if( !ifs && ifs2)
		{
			cout << " \n ========================================================" << endl;
			cout << " Can't find file: INPUT (used to generate local orbitals)." << endl;
			cout << " Find file: INPUTs (used to generate overlap data)" << endl;
			cout << " Calculate the overlap now." << endl;
			cout << " ========================================================" << endl;
		}
		if (!ifs && !ifs2)
        {
            WARNING_QUIT("Read_INPUT","can't find either file: INPUT or INPUTs");
        }

		if( ifs )
		{
			// mohan add 2010-05-01
			if ( SCAN_BEGIN(ifs,"<BANDS>",1,0) )
			{
				READ_VALUE( ifs, BANDS_CONTROL);
				if (BANDS_CONTROL)
				{
					READ_VALUE(ifs, BANDS_START);
					assert(BANDS_START > 0);
					BANDS_START -= 1;
					// because the band index start from 0;
					READ_VALUE(ifs, BANDS_END);
					// BANDS_END need not -=1;
				}

				ofs_running << " ATTENTION! YOU SELECT THE BANDS: " << endl;
				ofs_running << " BANDS_START=" << BANDS_START+1 << " BANDS_END=" << BANDS_END << endl;
			}

			// mohan add 2009-08-26
			if ( SCAN_BEGIN(ifs,"<BLOCK_NE>",1,0) )
			{
				READ_VALUE( ifs, BLOCK_NE);
				assert(BLOCK_NE>=0);

				READ_VALUE( ifs, BLOCK_NE_MIN);
				assert(BLOCK_NE >= BLOCK_NE_MIN);
				// static_cast<int>( sqrt( 2.0 * ecut )* rcut/3.1415926535897932 );

				ofs_running << " ATTENTION! YOU SELECT BLOCK_NE=" << BLOCK_NE << endl;
				ofs_running << " ALSO BLOCK_NE_MIN=" << BLOCK_NE_MIN << endl;
			}

			// mohan add 2009-08-31
			if ( SCAN_BEGIN(ifs,"<BLOCK_NE_ENERGY>",1,0) )
			{
				double block_ne_energy;
				double rcut_in;
				READ_VALUE( ifs, block_ne_energy);
				READ_VALUE( ifs, rcut_in );
				BLOCK_NE = static_cast<int>( sqrt( 2.0 * block_ne_energy )* rcut_in/3.1415926535897932 );
				assert(BLOCK_NE >= 0);

				ofs_running << "\n Attential! Now, input BLOCK_NE parameter be replaced!" << endl;
				ofs_running << "\n New BLOCK_NE = " << BLOCK_NE
					<< " energy_cut = " << block_ne_energy
					<< " rcut_in = " << rcut_in << endl;
			}

			// mohan add 2009-08-28
			// 1: max spillage value.
			// 2: average spillage value.
			// 3: add kapa value.
			if ( SCAN_BEGIN(ifs,"<SCHEME_VALUE>",1,0) )
			{
				READ_VALUE( ifs, SCHEME_VALUE );
				assert(abs(SCHEME_VALUE)>0);
				assert(abs(SCHEME_VALUE)<4);
			}

			if ( SCAN_BEGIN(ifs,"<PARALLEL>",1,0) )
			{
				READ_VALUE( ifs, NKSTOT );
				assert(NKSTOT>0);
#ifdef __MPI
				READ_VALUE( ifs, KPAR );
				assert(KPAR>0);
#endif	
			}

		}// end ifs

		else if( ifs2 )
		{
			if ( SCAN_BEGIN(ifs2,"<FILE>",1,0) )
			{
				READ_VALUE( ifs2, WFC_FILE);
			}
			
			if ( SCAN_BEGIN(ifs2,"<SPHERICAL_BESSEL>",1,0) )
			{
				READ_VALUE( ifs2, SMOOTH);
				assert(SMOOTH == 0 || SMOOTH == 1);
				READ_VALUE( ifs2, SIGMA);
				assert(SIGMA >= 0.0);
				READ_VALUE( ifs2, ECUT_JLQ);
				assert(ECUT_JLQ > 0 && ECUT_JLQ < 1000);
				// mohan set the boundaries)
				READ_VALUE( ifs2, RCUT);
				assert(RCUT > 0 && RCUT < 15);
				READ_VALUE( ifs2, TOLERENCE);	
				assert(TOLERENCE> 0.0);
			}
			
            if ( SCAN_BEGIN(ifs2,"<PARALLEL>",1,0) )
            {
                READ_VALUE( ifs2, NKSTOT );
                assert(NKSTOT>0);
#ifdef __MPI
                READ_VALUE( ifs2, KPAR );
                assert(KPAR>0);
#endif
            }
		}// end ifs2
    } // end my rank0.

    //------------------------------------------------
	// second part
    //------------------------------------------------
	// 2.1
	if( ifs )
	{
		if (MY_RANK==0)
		{
			this->read_PW_QS_1();
		}
	}
	else if ( ifs2 )
	{
		this->read_WFC();
	}

	// 2.2
	// bcast NKSTOT.
	// bcast the global information.
	// need to first bcast KPAR.
	this->bcast();

	// set the parallel information 
	// about pool.
	Pkpoints.init();

	// init the k point information.
	Pkpoints.kinfo(NKSTOT);

	// calculate the local number of k points.
#ifdef __MPI
	NKS = Pkpoints.nks_pool[MY_POOL];
#else
	NKS = NKSTOT;
#endif

	OUT(ofs_running,"NKS",NKS);

	// 2.3
	if( ifs )
	{
		this->read_PW_QS_2();
	}
	else if( ifs2 )
	{
		this->read_WFC_2();
		USEPW = 1;	
	}

    return;
}

void Read_INPUT::bcast(void)
{
#ifdef __MPI
    Parallel_Common::bcast_bool( BANDS_CONTROL );
    Parallel_Common::bcast_int( BANDS_START );
    Parallel_Common::bcast_int( BANDS_END );

    Parallel_Common::bcast_int( BLOCK_NE );
    Parallel_Common::bcast_int( BLOCK_NE_MIN );

    Parallel_Common::bcast_int( SCHEME_VALUE );
    Parallel_Common::bcast_int( KPAR );

    Parallel_Common::bcast_int( CALSPI );
    Parallel_Common::bcast_bool( RESTART );
    Parallel_Common::bcast_bool( this->output );
    Parallel_Common::bcast_int( STRNUM );

	OUT(ofs_running,"STRNUM",STRNUM);

    if (MY_RANK!=0)
    {
        delete[] qsfile;
        qsfile = new string[STRNUM];
    }

    for (int i=0; i<STRNUM; i++)
    {
        Parallel_Common::bcast_string( qsfile, STRNUM );
    }

	// information about NTYPE...
	Parallel_Common::bcast_int( NTYPE );
	if(MY_RANK!=0)
	{
		LABEL = new string[NTYPE];
		NA = new int[NTYPE];
	}

	Parallel_Common::bcast_string( LABEL, NTYPE);
	Parallel_Common::bcast_int( NA, NTYPE);

	if(MY_RANK!=0)
	{
		CARPOSX = new double*[NTYPE];
		CARPOSY = new double*[NTYPE];
		CARPOSZ = new double*[NTYPE];
		for(int it=0; it<NTYPE; it++)
		{
			CARPOSX[it] = new double[NA[it]];
			CARPOSY[it] = new double[NA[it]];
			CARPOSZ[it] = new double[NA[it]];
		}
	}

	for(int it=0; it<NTYPE; it++)
	{
		Parallel_Common::bcast_double(CARPOSX[it], NA[it]);
		Parallel_Common::bcast_double(CARPOSY[it], NA[it]);
		Parallel_Common::bcast_double(CARPOSZ[it], NA[it]);
	}
	
	Parallel_Common::bcast_double( LAT0 );
	
	Parallel_Common::bcast_double( LATVEC.e11 );
	Parallel_Common::bcast_double( LATVEC.e12 );
	Parallel_Common::bcast_double( LATVEC.e13 );
	Parallel_Common::bcast_double( LATVEC.e21 );
	Parallel_Common::bcast_double( LATVEC.e22 );
	Parallel_Common::bcast_double( LATVEC.e23 );
	Parallel_Common::bcast_double( LATVEC.e31 );
	Parallel_Common::bcast_double( LATVEC.e32 );
	Parallel_Common::bcast_double( LATVEC.e33 );

	// about ecut...
	Parallel_Common::bcast_double( ECUT );
	Parallel_Common::bcast_double( ECUT_JLQ );
	Parallel_Common::bcast_double( RCUT );
	Parallel_Common::bcast_bool( SMOOTH );
	Parallel_Common::bcast_double( SIGMA );
	Parallel_Common::bcast_double( TOLERENCE );
	Parallel_Common::bcast_int( LMAXALL );
	Parallel_Common::bcast_int( NKSTOT );
	Parallel_Common::bcast_int( NBANDS );
	Parallel_Common::bcast_int( NWFCALL );
	Parallel_Common::bcast_int( NE );

#endif
    return;
}

void Read_INPUT::read_PW_QS_2(void)
{
    TITLE("Read_INPUT","read_PW_QS_2");
    // SV is one member of class 'SpillageValue'.
    this->SV.allocate( STRNUM );
    delete[] this->QS_data;
    this->QS_data = new ReadData[ STRNUM ];

    for (int i=0; i<STRNUM; i++)
    {
        // don't change the order !! Because first readin
        // nks, nbands, nwfc, ne
        this->QS_data[i].OverlapQandS( qsfile[i].c_str() );
        //this->QS_data[in].OverlapSinv( sinv_file.c_str() ); // not used
    }

    return;
}

void Read_INPUT::read_PW_QS_1(void)
{
    TITLE(ofs_running, "Read_INPUT","read_PW_QS_1");

    if ( SCAN_BEGIN(ifs,"<PW_QS>") )
    {
        // <1> calculate spillage or not.
        READ_VALUE( ifs, CALSPI );
        // <1.5> restart from file or not.
        // mohan add this function 2009-08-27
        READ_VALUE( ifs, RESTART);
        // <2> output file or not.
        READ_VALUE( ifs, this->output );
        // <3> number of structures.
        READ_VALUE( ifs, STRNUM );
		OUT(ofs_running,"STRNUM",STRNUM);
        if (STRNUM<=0)
        {
            WARNING_QUIT("Read_INPUT::read_PW_QS_1","STRNUM<=0");
        }

        delete[] qsfile;
        this->qsfile = new string[STRNUM];

        // <4> read in each file.
        for (int in=0; in< STRNUM; in++)
        {
            READ_VALUE(ifs, qsfile[in]);
            ifstream check(qsfile[in].c_str());
            if (!check)
            {
                ofs_running << qsfile[in] << endl;
                WARNING_QUIT("Read_INPUT","the file is not exist!");
            }
            else
            {
				OUT(ofs_running,"FILE",qsfile[in]);
            	static bool already_read = false;
            	if (!already_read)
				{
					check >> LAT0;
					OUT(ofs_running,"LAT0",LAT0);
					check >> LATVEC.e11 >> LATVEC.e12 >> LATVEC.e13;
					check >> LATVEC.e21 >> LATVEC.e22 >> LATVEC.e23;
					check >> LATVEC.e31 >> LATVEC.e32 >> LATVEC.e33;
					OUT(ofs_running,"a1",LATVEC.e11,LATVEC.e12,LATVEC.e13);
					OUT(ofs_running,"a2",LATVEC.e21,LATVEC.e22,LATVEC.e23);
					OUT(ofs_running,"a3",LATVEC.e31,LATVEC.e32,LATVEC.e33);
					READ_VALUE(check, NTYPE);
					OUT(ofs_running,"NTYPE",NTYPE);
					assert( NTYPE > 0 );
					NA = new int[NTYPE];
					LABEL = new string[NTYPE];
					CARPOSX = new double*[NTYPE];
					CARPOSY = new double*[NTYPE];
					CARPOSZ = new double*[NTYPE];
					for(int it=0; it<NTYPE; it++)
					{
						READ_VALUE(check, LABEL[it]);
						READ_VALUE(check, NA[it]);
						OUT(ofs_running,"LABEL",LABEL[it]);
						OUT(ofs_running,"NA",NA[it]);
						CARPOSX[it] = new double[NA[it]];
						CARPOSY[it] = new double[NA[it]];
						CARPOSZ[it] = new double[NA[it]];
						for(int ia=0; ia<NA[it]; ia++)
						{
							check >> CARPOSX[it][ia] >> CARPOSY[it][ia] >> CARPOSZ[it][ia];
							OUT(ofs_running,"POS",CARPOSX[it][ia],CARPOSY[it][ia],CARPOSZ[it][ia]);
						}
					}

					READ_VALUE(check, ECUT);
					READ_VALUE(check, ECUT_JLQ);
					READ_VALUE(check, RCUT);
					READ_VALUE(check, SMOOTH);
					READ_VALUE(check, SIGMA);

					READ_VALUE(check, TOLERENCE);
					READ_VALUE(check, LMAXALL);
					READ_VALUE(check, NKSTOT);
					READ_VALUE(check, NBANDS);
					READ_VALUE(check, NWFCALL);
					READ_VALUE(check, NE); 

					OUT(ofs_running,"ECUT(Ry)",ECUT);
					OUT(ofs_running,"ECUT_JLQ(Ry)",ECUT_JLQ);
					OUT(ofs_running,"RCUT(Bohr)",RCUT);
					OUT(ofs_running,"SMOOTH",SMOOTH);
					OUT(ofs_running,"SIGMA",SIGMA);
					OUT(ofs_running,"TOLERENCE",TOLERENCE);
					OUT(ofs_running,"LMAXALL",LMAXALL);
					OUT(ofs_running,"NKSTOT",NKSTOT);
					OUT(ofs_running,"NBANDS",NBANDS);
					OUT(ofs_running,"NWFCALL",NWFCALL);
					OUT(ofs_running,"NE",NE);

					if(SMOOTH) assert(SIGMA!=0.0);

                	already_read = true;
				}//only read once.
                check.close();
            }
        }
    }
    return;
}

void Read_INPUT::read_WFC(void)
{
    TITLE("Read_INPUT","read_WFC");

	ifstream ifswfc( WFC_FILE.c_str() );

	if(!ifswfc)
	{
		cout << " Can't find wfc file : " << WFC_FILE << endl;
		WARNING_QUIT("Read_INPUT::read_WFC","Can't find wavefunctions file.");
	}
	else
	{
		cout << " Find file: " << WFC_FILE << endl;
	}

	if(MY_RANK==0)
	{

		ifswfc >> LAT0;
		cout << "\n LAT0 = " << LAT0;
		ifswfc >> LATVEC.e11 >> LATVEC.e12 >> LATVEC.e13;
		ifswfc >> LATVEC.e21 >> LATVEC.e22 >> LATVEC.e23;
		ifswfc >> LATVEC.e31 >> LATVEC.e32 >> LATVEC.e33;
		cout << "\n a1 = " << LATVEC.e11 << " " << LATVEC.e12 << " " << LATVEC.e13; 
		cout << "\n a2 = " << LATVEC.e21 << " " << LATVEC.e22 << " " << LATVEC.e23; 
		cout << "\n a3 = " << LATVEC.e31 << " " << LATVEC.e32 << " " << LATVEC.e33; 
		READ_VALUE(ifswfc, NTYPE);
		cout << "\n NTYPE = " << NTYPE << endl;
		assert( NTYPE > 0 );
		NA = new int[NTYPE];
		LABEL = new string[NTYPE];
		CARPOSX = new double*[NTYPE];
		CARPOSY = new double*[NTYPE];
		CARPOSZ = new double*[NTYPE];
		for(int it=0; it<NTYPE; it++)
		{
			READ_VALUE(ifswfc, LABEL[it]);
			READ_VALUE(ifswfc, NA[it]);
			cout << "\n LABEL = " << LABEL[it];
			cout << "\n NA=" << NA[it];
			CARPOSX[it] = new double[NA[it]];
			CARPOSY[it] = new double[NA[it]];
			CARPOSZ[it] = new double[NA[it]];
			for(int ia=0; ia<NA[it]; ia++)
			{
				ifswfc >> CARPOSX[it][ia] >> CARPOSY[it][ia] >> CARPOSZ[it][ia];
				cout << "\n " << CARPOSX[it][ia] << " " << CARPOSY[it][ia] << " " << CARPOSZ[it][ia];
			}
		}

		READ_VALUE(ifswfc, ECUT);
		READ_VALUE(ifswfc, NKSTOT);
		READ_VALUE(ifswfc, NBANDS);

		cout << "\n ECUT = " << ECUT;
		cout << "\n ECUT_JLQ = " << ECUT_JLQ;
		cout << "\n RCUT = " << RCUT;
		cout << "\n SMOOTH = " << SMOOTH;
		cout << "\n SIGMA = " << SIGMA;
		cout << "\n TOLERENCE = " << TOLERENCE;
		cout << "\n NKSTOT = " << NKSTOT;
		cout << "\n NBANDs = " << NBANDS;
	}

	ifswfc.close();
	return;
}


void Read_INPUT::read_WFC_2(void)
{
	TITLE("Read_INPUT","read_WFC_2");
	QUIT();
}


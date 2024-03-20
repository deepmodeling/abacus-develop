#include "MultiZeta.h"
#include "SpillageStep.h"
#include "Out_Orbital.h"
#include "../src_parallel/parallel_common.h"
#include "tools.h"

MultiZeta::MultiZeta()
{
    test = 0; // output function name
    Level = new SpillageStep[1];
    lmax_type = new int[1];
}

MultiZeta::~MultiZeta()
{
    delete[] Level;
    delete[] lmax_type;
}

// be called in main.cpp program.
void MultiZeta::init(void)
{
    TITLE(ofs_running,"MultiZeta","init");
    //=============================
    // nlevel : optimization levels
    //=============================
    for (int is=0; is<nlevel; is++)
    {
        // the spillage value of each level is "smaller" than the previous level,
        // because we get the "left" Bloch wave functions each time.
        cout << "\n\n =================================================== ";
        cout << "\n                       level = " << is+1;
        cout << "\n =================================================== " << endl;

        ofs_running << "\n\n =================================================== ";
        ofs_running << "\n                       level = " << is+1;
        ofs_running << "\n =================================================== " << endl;

        this->ilevel = is+1;// 0 means level 1.

        // mohan add 2009-09-25
        //input.Coef.accumulate_num_zero();

        // if start temperature is different for each level
        // ( I don't think this is a very good idea for high level orbitals!!
        // mohan note:2009-8-20, it's better to make temperature from "high enough" to
        // "low enough" for each level)
        metro.reset_temperature( is );

        this->Level[is].ilevel = is;
        this->Level[is].set_nchi();
        this->Level[is].set_nwfc();
        this->Level[is].set_iw00();// to get Sq1q2
        this->Level[is].allocate_data( 0 ); // to get Qin


        // if is=0, means it's the first time, spillage0 = 1.0;
        if (is==0)
        {
            for (int istr=0; istr<STRNUM; istr++)
            {
                this->Level[is].data[istr].spillage0 = 1.0;
            }
        }
        else if (is>0)
        {
            for (int istr=0; istr<STRNUM; istr++)
            {
                this->Level[is].data[istr].spillage0 = input.SV.value_old[istr];
            }
            cout << "\n Orthogonal......." << endl;
            ofs_running << "\n Orthogonal......." << endl;
            //****************************************
            // (2): Orthogonal to get new Bloch space!
            // This subroutine is very important.
            //****************************************
            Orth.start( this->Level[is-1], this->Level[is] );
        }

		if(USEPW)
		{
			PW.allocate_psi1d(is);
			PW.allocate_psi3d(is);
		}

        for (int istr=0; istr<STRNUM; istr++)
        {
            this->Level[is].data[istr].initQS();
            this->Level[is].init_QS_matrix(istr);
            input.SV.value_old[istr] = this->Level[is].get_spillage(istr, 0, 0);
        }

        cout << "\n Initial Spillage Value at this Level: " << endl;
        // output spillage value
        input.SV.out();
        cout << endl;

        //*******************************************************************
        // main part: Use Metropolis algorithm, move different temperature,
        //*******************************************************************
        this->metro.move_various_temperature( this->Level[is]);

        // save spillage value for each level
        input.SV.save_level(is);

        this->saveC();

    }

    return;;
}



// read in paramters from BLOCK <OPTIMIZE>, called by MultiZeta::init();
void MultiZeta::set_levels(void)
{
    TITLE(ofs_running,"MultiZeta","set_levels");

	bool begin = false;

	if(MY_RANK==0)
	{
		if ( SCAN_BEGIN(input.ifs, "<OPTIMIZE>") ) begin = true;
	}

#ifdef __MPI
	Parallel_Common::bcast_bool(begin);
#endif

    // read in parameters in BLOCK : 'OPTIMIZE'
    if ( begin )
    {
        // (1) read in important parameter: nlevel.
        // nlevel means the number of orbital shells we need.
        
		if(MY_RANK==0)
		{
			READ_VALUE(input.ifs, this->nlevel);
        	assert(nlevel>-1);
        	assert(nlevel<100);
		}
	
#ifdef __MPI	
		Parallel_Common::bcast_int(nlevel);
#endif
			
        input.nlevel = this->nlevel; // mohan add 2009-06-09
//        cout << " nlevel = " << nlevel << endl;
		OUT(ofs_running,"nlevel",nlevel);

        // (2) init 'SpillageStep' class, each level has
        // different parameters.
        delete[] this->Level;
        this->Level = new SpillageStep[nlevel];

        string useless;

		if(MY_RANK==0)
		{
        	READ_VALUE(input.ifs, useless);
		}

#ifdef __MPI
		Parallel_Common::bcast_string(useless);
#endif

        // (3) read in information concerning each level.
        for (int i=0; i<nlevel; i++)
        {
            ofs_running << "\n >>>>> Level=" << i+1 ;

            this->Level[i].set_info(NTYPE);
            for (int it=0; it<NTYPE; it++)
            {
                // Read in information for each level.
                // (a) atom id (element periodic table).
                // (b) number of atoms for each atom type.
                // (c) lmax used in this Level.
        
				if(MY_RANK==0)
				{
		        	input.ifs >> Level[i].info[it].id //id for this type
                	>> Level[i].info[it].na //number of atoms.
                	>> Level[i].info[it].state// new or skip, mohan 2009-08-27
                	>> Level[i].info[it].lmax;//lmax in this level
				}

#ifdef __MPI
				Parallel_Common::bcast_int(Level[i].info[it].id);
				Parallel_Common::bcast_int(Level[i].info[it].na);
				Parallel_Common::bcast_string(Level[i].info[it].state);
				Parallel_Common::bcast_int(Level[i].info[it].lmax);
#endif

                if (!RESTART && Level[i].info[it].state=="skip")
                {
                    ofs_running << "\n Level = " << i;
                    ofs_running << "\n it = " << it;
                    WARNING_QUIT("MultiZeta::set_levels","skip is not available if RESTART = false.");
                }
                if ( Level[i].info[it].state!="new" && Level[i].info[it].state!="skip")
                {
                    ofs_running << "\n Level = " << i;
                    ofs_running << "\n it = " << it;
                    ofs_running << "\n state = " << Level[i].info[it].state << endl;
                    WARNING_QUIT("MultiZeta::set_levels","new or skip?");
                }

                ofs_running << "\n id=" << Level[i].info[it].id;
                ofs_running << " na=" << Level[i].info[it].na;
                ofs_running << " lmax=" << Level[i].info[it].lmax;

				if(Level[i].info[it].id==0)
				{
					ofs_running << " The id of elements should be numbers!" << endl;
					ofs_running << " check the id of elements." << endl;
					WARNING_QUIT("MultiZeta::set_levels","id=0");
				}

                // lmax=0(s), 1(p), 2(d), 3(f), 4(g)
                if (Level[i].info[it].lmax>=5 || Level[i].info[it].lmax<0)
                {
                    ofs_running << "\n lmax=" << Level[i].info[it].lmax << endl;
                    WARNING_QUIT("MultiZeta::set_levels"," lmax out of range.");
                }

                // 'info' stands for "Type_Information"
                // init: n, nbase;
                // allocate
                Level[i].info[it].init();

                for (int L=0; L<Level[i].info[it].lmax+1; L++)
                {

                    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    //	read in nchi for each L and full/average orbital.
                    int nz;
                    string fain;
                    string combine;

                    // read in characters such as : 1/a
                    input.ifs >> combine;

#ifdef __MPI
                  	Parallel_Common::bcast_string(combine);
#endif
				  
				    //ofs_running << "\n lenght = " << combine.length() << endl;

                    // firstly, we need to find the position of the label: '/'
                    const int index1 = combine.find_first_of("/");
                    //ofs_running << "\n index1 = " << index1 << endl;

                    // if index1==-1, means the combine didn't contain
                    // '/' label
                    if (index1 == -1)
                    {
                        // f stands for 'full'
                        fain = "f";
                        // this 'atoi' command change string to int type.
                        nz = std::atoi( combine.c_str() );
                    }
                    else // mean containing "/"
                    {
                        nz = std::atoi( combine.substr(0,index1+1).c_str() );
                        // substrate the front int and remain the character.
                        // the second parameter '1' means the position is 'index+2'
                        fain = combine.substr(index1+1, 1);
                    }

                    ofs_running << "\n L=" << L;
                    ofs_running << " Zeta=" << nz;
                    ofs_running << " fain(FULL/AVERAGE)=" << fain;
                    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                    // assert number of zeta < 20, you can change on your own,
                    // typically, we only need nz=2 or 3.
                    assert( nz < 20 );
                    Level[i].info[it].n[L] = nz;
                    Level[i].info[it].fa[L] = fain;

                    // mohan added 2009-05-22
                    // i stands for level
                    // nbase stands for the radial wave functions before
                    // this level for each L of each atom type.
                    // (hard to understand, right? Try to figure out!)
                    // this parameter is very useful.
                    // nbase provide a base to count radial wave functions.
                    if (i==0)
                    {
                        Level[i].info[it].nbase[L] = 0;
                    }
                    else if (i>0)
                    {
                        for (int j=0; j<i; j++)
                        {
                            if ( L < Level[j].info[it].lmax+1 )
                            {
                                Level[i].info[it].nbase[L] = Level[j].info[it].nbase[L] + Level[j].info[it].n[L];
                            }
                        }
                    }
                    ofs_running << " start_n(" << L << ")=" << Level[i].info[it].nbase[L]+1;
                }// end L

                Level[i].info[it].cal_nmax();

                ofs_running << endl;
            }// end Type
        }// end Level
    }// end <OPTIMIZE>
    return;
}


// be called by MultiZeta::init()
void MultiZeta::cal_lmax_nmax(void)
{
    TITLE(ofs_running,"MultZeta","cal_lmax_nmax");
    //====================================
    // (1) calculate LMAXUSED and NMAXUSED
    //====================================
    // initialize, LMAXUSED will save
    // the max L in all levels and
    // all atom types.

    // initialize, lmax_type will save
    // the max L for each atom type
    // in all levels.
    delete[] lmax_type;
    this->lmax_type = new int[NTYPE];
    for (int it=0; it<NTYPE; it++)
    {
        this->lmax_type[it] = -1;
    }

    // sum information from all levels.
	NMAXUSED = 0;
    for (int k=0; k<this->nlevel; k++)
    {
        for (int it=0; it< NTYPE; it++)
        {
            this->lmax_type[it] = std::max( this->lmax_type[it], this->Level[k].info[it].lmax );
            LMAXUSED = std::max( LMAXUSED, this->lmax_type[it] );
            NMAXUSED = std::max( NMAXUSED, this->Level[k].info[it].nmax );
        }
    }

    // for test.
    if (test==2)
    {
        for (int it=0; it<NTYPE; it++)
        {
            cout << "\n lmax_type(" << it+1 << ")=" << lmax_type[it];
        }
        cout << "\n LMAXUSED=" << LMAXUSED;
        cout << "\n NMAXUSED=" << NMAXUSED;
    }

    //=============================
    // (2) calculate l_nchi
    //=============================
    this->l_nchi.create(NTYPE, LMAXUSED+1);

    // after sum up information from all levels,
    // we will know clearly hwo many radial wave functions
    // for each L of each atom type.
    for (int it=0; it<NTYPE; it++)
    {
        for (int L=0; L<lmax_type[it]+1; L++)
        {
            for (int k=0; k<this->nlevel; k++)
            {
                if (L < this->Level[k].info[it].lmax+1)
                {
                    this->l_nchi(it, L) += this->Level[k].info[it].n[L];
                }
            }
            ofs_running << " Type=" << it+1 << " L=" << L << " Zeta=" << l_nchi(it,L) << endl;
        }
    }

    // check for safety.
    for (int is=0; is<STRNUM; is++)
    {
        if ( LMAXUSED > LMAXALL )
        {
            cout << "\n In structure " << is;
            cout << "\n (1) LMAXUSED in MultiZeta.cpp = " << LMAXUSED;
            cout << "\n (2) LMAXALL in QS_data = " << LMAXALL;
            WARNING_QUIT("SpillageStep::get_data","More lmax required in INPUT than in PW output file!");
        }
    }

    //=============================
    // (3) NCHIUSED
    //=============================
    // calculate the total number of radial wave functions.
    NCHIUSED=0;
    for (int it=0; it<NTYPE; it++)
    {
        for (int L=0; L<lmax_type[it]+1; L++)
        {
            for (int N=0; N<l_nchi(it,L); N++)
            {
                ++NCHIUSED;
            }
        }
    }
    ofs_running << " N Zeta(including all atom species) = " << NCHIUSED << endl;
    return;
}




void MultiZeta::saveC(void)
{
    Out_Orbital Orb;
    Orb.write();
    return;
}



#include "Coefficients.h"
#include "../src_parallel/parallel_common.h"

bool Coefficients::write_c4_flag;

Coefficients::Coefficients()
{
    test = 0;
    output_c4_name = new string[1];
    accumulate_num = 0;//mohan add 2009-09-25
}

Coefficients::~Coefficients()
{
    if (TEST1) cout << "\n ~Coefficients()" << endl;
    delete[] output_c4_name;
}

void Coefficients::init(
    const int &ntype_in,
    const int &lmax_in,
    const int &nmax_in,
    const double &ecut_in,
    const double &ecut_jlq_in, // mohan add 2009-07-23
    const double &rcut_in,
    const double &tolerence_in)
{
    TITLE(ofs_running,"Coefficients", "Constructor");

    ntype = ntype_in;
    lmax = lmax_in;
    nmax = nmax_in;
    ecut = ecut_in;
    ecut_jlq = ecut_jlq_in; // mohan add 2009-07-23
    rcut = rcut_in;
//	in Rydberg Unit
    enumber = static_cast<int>( sqrt( ecut )* rcut/PI );
    enumber_jlq = static_cast<int> ( sqrt( ecut_jlq) * rcut/PI ); // mohan add 2009-07-23
    tolerence =  tolerence_in;

	OUT(ofs_running,"ntype",ntype);
	OUT(ofs_running,"lmax",lmax);
	OUT(ofs_running,"nmax",nmax);
	OUT(ofs_running,"ecut(Ry)",ecut);
	OUT(ofs_running,"ecut_jlq(Ry)",ecut_jlq);
	OUT(ofs_running,"rcut(Bohr)",rcut);
	OUT(ofs_running,"enumber",enumber);
	OUT(ofs_running,"enumber_jlq",enumber_jlq);
	OUT(ofs_running,"tolerence",tolerence);
	
    assert(ntype > 0);
    assert(lmax >= 0);
    assert(nmax > 0);
    assert(enumber > 0);
    assert(enumber_jlq > 0);

    assert(ecut > 0.0);
    assert(ecut_jlq > 0.0);
    assert(rcut > 0.0);
    assert(tolerence > 0.0);

    this->allocate_C4();

    delete[] output_c4_name;
    this->output_c4_name = new string[1];

    return;
}

void Coefficients::allocate_C4(void)
{
    if (test==1)TITLE("Coefficients", "allocate_C4");

    this->C4.create(ntype, lmax + 1, nmax, enumber);
    this->C4_old.create(ntype, lmax + 1, nmax, enumber);
    this->C4_accumulate.create(ntype, lmax + 1, nmax, enumber); // mohan add 2009-09-25
//	this->acc_count.create(ntype, lmax+1, nmax, enumber);
    this->accept_number.create(ntype, lmax+1, nmax, enumber); // mohan add 2009-10-16
    this->accept_rate.create(ntype, lmax+1, nmax, enumber);// mohan add 2009-10-31

    // init accept_rate
    for (int i=0; i<accept_rate.getSize(); i++)
    {
        this->accept_rate.ptr[i] = 1.0;
    }

    for (int it = 0; it < ntype; it++)
    {
        for (int il = 0; il < lmax+1; il++)
        {
            for (int in = 0; in < nmax; in++)
            {
                ofs_running << "\n Init C4: T=" << it+1 << " L=" << il << " N=" << in;
                for (int ie = 0; ie < enumber; ie++)
                {
                    // mohan modify 2009-08-26
                    int ne_cut;
                    int ne_cut_min = -1; // must initalized
                    if (BLOCK_NE_MIN>0)
                    {
                        ne_cut_min = BLOCK_NE_MIN;
                    }

                    if (BLOCK_NE<0)
                    {
                        ne_cut = enumber_jlq;
                    }
                    else
                    {
                        ne_cut = BLOCK_NE;
                    }

                    if (ie<ne_cut && ie>ne_cut_min)
                    {
                        this->C4(it, il, in, ie) = Random::between0and1();
                        this->C4_old(it, il, in, ie) = this->C4(it, il, in, ie);
                    }
                    else
                    {
                        this->C4(it, il, in, ie) = 0.0;
                        this->C4_old(it, il, in, ie) = 0.0;
                    }
                    if (ie%4==0)ofs_running<<endl;
                    ofs_running << setw(20) << C4_old(it, il, in, ie);
                }
            }
        }
    }

    // mohan add 2009-08-27, update 2010-04-17
    // read in C4 from ORBITALR_RESULTS.txt file.
    // make clear in the <C4> BLOCK,
    // how many C4 you want to read in!
    // remain bug when ntype > 1.
    if (RESTART)
    {
        if (MY_RANK==0)
        {
            ifstream ifs( "ORBITAL_RESULTS.txt" );
            if (!ifs)
            {
                WARNING_QUIT("Coefficients::allocate_C4()","can't find the restart file: ORBITAL_RESULTS.txt");
            }
            else
            {
                cout << "\n Read in RESTART information from : ORBITAL_RESULTS.txt" << endl;
            }

            double ecut_in;
            double rcut_in;
            int ne_in;
            if ( SCAN_BEGIN(ifs,"<INPUTS>") )
            {
                READ_VALUE(ifs, ecut_in);
                READ_VALUE(ifs, rcut_in);
                READ_VALUE(ifs, ne_in);
                if (ne_in != enumber)
                {
                    cout << "\n ne_in = " << ne_in;
                    cout << "\n enumber = " << enumber;
                    WARNING_QUIT("Coefficients::Allocate_C4()","enumber != ne_in");
                }
            }

            if ( SCAN_BEGIN(ifs,"<Coefficient>") )
            {
                int nchi_in;
                READ_VALUE(ifs, nchi_in);
                assert(nchi_in > 0);

                string no_use1, no_use2, no_use3;
                int it2, il2, in2;
                bool hold_on = false;
                int count_nchi = 0;
                for (int it = 0; it < ntype; it++)
                {
                    for (int il = 0; il < lmax+1; il++)
                    {
                        for (int in = 0; in < nmax; in++)
                        {
                            if (!hold_on)
                            {
                                ifs >> no_use1 >> no_use2 >> no_use3;
                                assert(no_use1 == "Type");
                                assert(no_use2 == "L");
                                assert(no_use3 == "Zeta-Orbital");
                                ifs >> it2>> il2>> in2;
								it2 -= 1;
								in2 -= 1;
                            }
                            if (it2!=it || il2!=il || in2!=in)
                            {
                                hold_on = true;
                                continue;
                            }
                            else
                            {
                                hold_on = false;
                                cout << "\n\n Read in T=" << it+1 << " L=" << il << " Zeta-Orbital=" << in;
                                for (int ie=0; ie<ne_in; ie++)
                                {
                                    ifs >> C4(it, il, in, ie);
                                    C4_old(it, il, in, ie) = C4(it, il, in, ie);
                                    if (ie%4==0) cout << endl;
                                    cout 
										<< setiosflags(ios::fixed) 
										<< setprecision(20) 
										<< setiosflags(ios::showpoint)
										<< setw(25) << C4_old(it, il, in, ie);
                                }
                                ++count_nchi;
                                // if count_nchi == nchi_in, all needed C4 have been read in!
                                if (count_nchi==nchi_in) return;
                            }
                        }
                    }
                }
            }// end SCAN_BEGIN
        }// end MY_RAN==0
    } // end RESTART

#ifdef __MPI
	Parallel_Common::bcast_double(C4.ptr, C4.getSize());
	Parallel_Common::bcast_double(C4_old.ptr, C4_old.getSize());
	ofs_running << "\n bcast C4 and C4_old done." << endl;
#endif
	
    return;
}

void Coefficients::copy_c4_to_old(void)
{
    for (int i=0; i<this->C4.getSize(); i++)
    {
        this->C4_old.ptr[i] = this->C4.ptr[i];
    }
    return;
}


void Coefficients::trial_c4( const int &t, const int &l, const int &n, const int &ie )
{
    // (1)
    this->C4_old(t, l, n, ie) = this->C4(t, l, n, ie);
    // update between -0.1~0.1

    // (2)
    this->C4(t, l, n, ie) += Random::betweenMinus1and1()/10.0* this->accept_rate(t, l, n, ie);

#ifdef __MPI
	// don't need to bcast ! it's a random seed.
//	ofs_running << "\n C4=" << this->C4(t, l, n, ie);
//	Parallel_Common::bcast_double(  this->C4(t, l, n, ie) );
#endif

//	if( abs(this->C4(t, l, n, ie)) < 1.0e-3) //mohan modify 2009-10-10
//	{
//		this->C4(t, l, n, ie) = 0.0;
//	}
//	this->C4(t, l, n, ie) += Random::betweenMinus1and1();

    //  this->C4(t, l, n, ie) *= Random::betweenMinus2and2();
    //  cout << "\n ic = " << ic << " ie = " << ie << "  C4 = " << C4(t, l, n, ie);
    return;
}

/*
void Coefficients::accumulate_num_zero(void)
{
	cout << "\n Accumulate_num = " << this->accumulate_num;
	this->accumulate_num = 0;
}
*/

/*
void Coefficients::accumulating_C4( const int &t, const int &l, const int &n, const int &ie )
{
	// mohan add 2009-09-25
	double mix = this->C4_accumulate(t, l, n, ie) * this->acc_count(t,l,n,ie) + this->C4_old(t, l, n, ie);
	++acc_count(t,l,n,ie);
	++accumulate_num;
    this->C4_accumulate(t, l, n, ie) = mix / acc_count(t,l,n,ie);
//	cout << "\n Accumulate_num = " << this->accumulate_num;

	return;
}
*/

void Coefficients::update_c4( const int &t, const int &l, const int &n, const int &ie )
{
    //===============================================
    // store the new c4 parameters(random accepted!!
    // then next c4 will replace by a new set.
    //===============================================

    this->C4_old(t, l, n, ie) = this->C4(t, l, n, ie);
    ++accept_number(t,l,n,ie);

    return;
}

void Coefficients::go_back_c4( const int &t, const int &l, const int &n, const int &ie )
{
    //==================================
    // give up the random added number;
    //==================================
    this->C4(t, l, n, ie) = this->C4_old(t, l, n, ie);
    return;
}

/*
void Coefficients::kinetic_energy( const int &t, const int &l, const int &n, const int &ie)
{
	if(!ke_flip_flag)
	{
		if(ke_up_flag)
		{
			// C4>0.0 and C4_old>0.0
			this->ke_total_up += ( C4(t, l, n, ie) - C4_old(t, l, n, ie) ) * ie * ie;
		}
		else
		{
			// C4<0.0 and C4_old>0.0
			this->ke_total_down += ( C4_old(t, l, n, ie) - C4(t, l, n, ie) ) * ie * ie;
		}
	}
	else
	{
		// (1) C4>0.0 and C4_old<0.0,
		// so ke_total_up energy increase
		// ke_total_down energy decrease
		//
		// (2) C4<0.0 and C4_old>0.0,
		// so ke_total_up energy decrease
		// and ke_total_down energy increase
		this->ke_total_up += C4(t, l, n, ie) * ie * ie;
		this->ke_total_down += C4_old(t, l, n, ie) * ie * ie;
	}

	this->ke_total_old = this->ke_total_new;
	this->ke_total_new = abs(abs(ke_total_up) - abs(ke_total_down));

//	cout << "\n ke_total_up = " << ke_total_up;
//	cout << "\n ke_total_down = " << ke_total_down;
//	cout << "\n ke_total_old = " << ke_total_old;
//	cout << "\n ke_total_new = " << ke_total_new << endl;

	return;
}
*/

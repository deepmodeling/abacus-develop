#include "SpillageStep.h"
#include "../src_parallel/parallel_reduce.h"
#include "tools.h"
#include "Inverse_Matrix_S.h"
#include "Inverse_Matrix_S_inline.h"
#include "../src_tools/complexmatrix_inline.h"
#include "../src_unittest/src_tools/common_test.h"

#ifdef _OPENMP
#include<omp.h>
#endif

SpillageStep::SpillageStep()
{
    test = 0;
    info = new Type_Information[1];
    data = new Step_Data[1];
    wayc = new Way2c4[1];
    wayd = new Way2Data[1];
}
SpillageStep::~SpillageStep()
{
    delete[] info;
    delete[] data;
    delete[] wayc;
    delete[] wayd;
}


void SpillageStep::set_info(const int &ntype_in)
{
    this->ntype = ntype_in;
    cout << "\n Element type = " << ntype;
    assert(ntype > 0);
    delete[] info;
    info = new Type_Information[ntype];
    return;
}

void SpillageStep::set_nwfc()
{
    TITLE("SpillageStep","set_nwfc");
    assert( ntype > 0);
    //===============================
    // (1) generate nwfc, nwfc2
    //===============================
    this->nwfc = 0;
    this->nwfc2 = 0;
    for (int i=0; i<ntype; i++)
    {
        // where to set info?
        // answer : in Multi_Zeta
        // realized in SpillageStep::set_info;
        info[i].cal_nw();
        nwfc += info[i].nw;

        // nwfc2 consider (it, ia, l, n, m) for this level.
        nwfc2 += info[i].nw2;
        cout << " Type=" << info[i].id << " Readin_Orbital=" << info[i].nw << " Generate_Orbital=" << info[i].nw2 << endl;
    }

    cout << " Total Localized Orbitals in this Level = " << nwfc2 << endl;

    //===============================
    // (2) generate iw0, type,
    // i, L, N, m
    //===============================
    delete[] wayd;
    wayd = new Way2Data[nwfc2];
    int iw=0;
    int iw2=0;
    for (int it=0; it<ntype; it++)
    {
        for (int i=0; i<info[it].na; i++)
        {
            for (int l=0; l<info[it].lmax+1; l++)
            {
                // In case we use multi-zeta,
                // the change then be false.
                bool change = true;
                int increase = 0;
                //cout << "\n" << " l=" << l << " n=" << info[it].n[l];
                for (int n=0; n<info[it].n[l]; n++)
                {
                    int copy = iw;
                    if (info[it].fa[l] == "f") //full
                    {
                        for (int m=0; m<2*l+1; m++)
                        {
                            wayd[iw2].type = it;
                            wayd[iw2].i = i;
                            wayd[iw2].L = l;
                            wayd[iw2].N = n + info[it].nbase[l]; // mohan modify 2009-05-22
                            wayd[iw2].m = m; // mohan add 2010-06-16
                            wayd[iw2].iw0 = copy;
                            wayd[iw2].average = false;
                            if (test == 2)
                            {
                                cout <<"\n" << " L=" << l << " N=" << wayd[iw2].N
                                     << " m=" << m << " iw2iw[" << iw2 << "]=" << copy;
                            }
                            ++iw2; // iw2 is the index of local basis.
                            ++copy;
                        }
                    }
                    else if (info[it].fa[l] == "a") //average
                    {
                        wayd[iw2].type = it;
                        wayd[iw2].i = i;
                        wayd[iw2].L = l;
                        wayd[iw2].N = n + info[it].nbase[l];
                        wayd[iw2].m = 0; // mohan add 2010-06-16, maybe still a bug
                        wayd[iw2].iw0 = copy;
                        wayd[iw2].average = true;
                        ++iw2;
                        ++copy;
                    }
                    increase = 2*l+1; // mohan modified 2009-07-12
                    // be careful,if n[l]==0, iw doesn't increase.
                }
                iw+=increase; // iw is the index of using PW data in this step.
            }
        }
    }

    cout << " iw = " << iw << " iw2 = " << iw2 << endl;

    if (iw!=nwfc)
    {
        cout << "\n iw = " << iw;
        cout << "\n nwfc = " << nwfc;
        exit(0);
    }
    assert(iw2==nwfc2);

    //===========================
    // (4) generate ic
    //===========================
    int ic=0;
    iw2=0;
    for (int it=0; it<ntype; it++)
    {
        int nw_this_type;
        int index_copy;
        for (int i=0; i<info[it].na; i++)
        {
            if (i==0)
            {
                index_copy = iw2;
                for (int l=0; l<info[it].lmax+1; l++)
                {
                    for (int n=0; n<info[it].n[l]; n++)
                    {
                        // condier full or average.
                        if (info[it].fa[l] == "f")
                        {
                            for (int m=0; m<2*l+1; m++)
                            {
                                this->wayd[iw2].ic = ic;
                                if (test==2) cout << "\n iw2=" << iw2 << " ic=" << ic;
                                ++iw2;
                            }// end m
                        }
                        else if (info[it].fa[l] == "a")
                        {
                            this->wayd[iw2].ic = ic;
                            ++iw2;
                        }
                        ++ic;
                    }// end n
                }// end l
                nw_this_type = iw2 - index_copy;
            }
            else
            {
                for (int j=index_copy; j<index_copy+nw_this_type; j++)
                {
                    this->wayd[iw2].ic = this->wayd[j].ic;
                    if (test==2) cout << "\n iw2=" << iw2 << " ic=" << wayd[iw2].ic;
                    iw2++;
                }
            }
        }
    }

    //cout << "\n iw2 = " << iw2;
    //cout << "\n nwfc2 = " << nwfc2;
    assert(iw2==nwfc2);
    return;
}


void SpillageStep::set_iw00(void)
{
    //=========================
    // (3) generate iw00
    //=========================

    int iw00 = 0;
    int iw2 = 0;
    for (int it=0; it<ntype; it++)
    {
        if ( info[it].na != NA[it] )
        {
            cout << "\n info[" << it << "].na=" << info[it].na;
            cout << "\n input.QS_data[0].na[" << it << "]=" << NA[it];
            WARNING_QUIT("SpillageStep::set_iw00","in <OPTIMIZE> not equal.");
        }
		//cout << " type=" << it << " lmax=" << info[it].lmax << endl;  
        for (int i=0; i< NA[it]; i++)
        {
            //for (int l=0; l<= mz.lmax_type[it]; l++)
            for (int l=0; l<= LMAXALL; l++) // mohan fix bug 2010-08-06
            {
                if ( l <= info[it].lmax )
                {
                    for (int n=0; n<info[it].n[l]; n++)
                    {
                        if ( info[it].fa[l] == "f" )
                        {
                            int copy = iw00;
                            for (int m=0; m<2*l+1; m++)
                            {
                                this->wayd[iw2].iw00 = copy;// position to get overlap.
                                if (test == 2)
                                {
                                    cout << "\n iw2=" << iw2 << " iw00=" << copy;
                                }
                                ++copy;
                                ++iw2;
                            }
                        }
                        else if ( info[it].fa[l] == "a" )
                        {
                            this->wayd[iw2].iw00 = iw00; // just the start position to get overlap.
                            ++iw2;
                        }
                    }
                }
                iw00 += 2*l + 1;
            }
        }
    }

    if (iw00!=NWFCALL)
    {
        cout << "\n iw00=" << iw00 << " NWFCALL=" << NWFCALL;
        WARNING_QUIT("SpillageStep::set_iw00","iw00!=input.QS_data[0].nwfc");
    }
    assert(iw2==nwfc2);

    return;
}


void SpillageStep::set_nchi()
{
    if (test==1)TITLE("SpillageStep","set_nchi");
    assert( ntype >0 );
    //=========================
    // (1) generate nchi
    //=========================
    this->nchi = 0;
    for (int i=0; i<ntype; i++)
    {
        for (int l=0; l<info[i].lmax+1; l++)
        {
            for (int j=0; j<info[i].n[l]; j++)
            {
                nchi++;
            }
        }
    }

    if (test == 2)
    {
        cout << "\n nchi for this step = " << nchi;
    }

    //==========================
    // (2) genreate wayc.type
    // L, N, used in Metropolis
    //==========================
    delete[] wayc;
    this->wayc = new Way2c4[nchi];

    int ichi=0;
    for (int i=0; i<ntype; i++)
    {
        for (int l=0; l<info[i].lmax+1; l++)
        {
            for (int j=0; j<info[i].n[l]; j++)
            {
                wayc[ichi].type = i;
                wayc[ichi].L = l;
                wayc[ichi].N = j+info[i].nbase[l];
                if (l==0) wayc[ichi].spd = 's';
                else if (l==1) wayc[ichi].spd = 'p';
                else if (l==2) wayc[ichi].spd = 'd';
                else if (l==3) wayc[ichi].spd = 'f';
                else if (l==4) wayc[ichi].spd = 'g';
                else if (l==5) wayc[ichi].spd = 'h'; // right?
                else if (l==6) wayc[ichi].spd = 'i'; // right?
                else wayc[ichi].spd = 'x'; // !!!
                ichi++;

                if (test==2)
                {
                    cout << "\n ichi=" << ichi
                         << " type=" << i
                         << " L=" << l
                         << " N=" << j;
                }
            }
        }
    }
    assert(ichi == nchi);

    return;
}

// be called in Multizeta::init(), before each step Metropolis begins.
void SpillageStep::allocate_data( const int &istep )
{
    if (test==1)TITLE("SpillageStep","allocate_data");
    assert( nwfc > 0 );

    delete[] this->data;
    this->data = new Step_Data[ STRNUM ];
    this->ne = NE;

    for ( int is=0; is<STRNUM; is++)
    {
        //==================================
        // copy the same thing from
        // QS_data to data :
        // nks, weight, nbands, ne
        // but nwfc and nwfc2 are different
        //==================================
        data[is].init(
            NKS,
            input.QS_data[is].weight,
            NBANDS,
            NE,
            this->nwfc,
            this->nwfc2,
            NWFCALL // all read in wave functions number
        );

        if (istep == 0)
        {
            for (int ik=0; ik<data[is].nks; ik++)
            {
                for (int ib=0; ib<data[is].nbands; ib++)
                {
                    for (int ie=0; ie<data[is].ne; ie++)
                    {
                        for (int iw=0; iw<NWFCALL; iw++)
                        {
                            data[is].Qin(ik,ib,iw,ie)=input.QS_data[is].Qin(ik,ib,iw,ie);
                        }
                    }
                }
            }
        }

        //========================================================
        // copy the needed part, but not enough for orthogonal
        //========================================================
        /*
        int iw=0;
        int iw_step=0;
        for(int it=0; it<ntype; it++)
        {
        	assert( info[it].na == input.QS_data[is].na[it] );
        	for(int i=0; i< input.QS_data[is].na[it]; i++)
        	{
        		for(int l=0; l< input.QS_data[is].lmax+1; l++)
        		{
        			for(int m=0; m<2*l+1; m++)
        			{
        				// a break..
        				if( l <= info[it].lmax )
        				{
        					if( info[it].n[l] > 0 )
        					{
        						for(int ik=0; ik< data[is].nks; ik++)
        						{
        							for(int ib=0; ib< data[is].nbands; ib++)
        							{
        								for(int ie=0; ie< data[is].ne; ie++)
        								{
        									data[is].Qin(ik, ib, iw_step, ie) = input.QS_data[is].Qin(ik, ib, iw, ie);
        								}
        							}
        						}
        						if(test==2)cout << "\n iw_step=" << iw_step << " iw=" << iw;
        						iw_step++;
        					}
        				}
        				iw++;
        			}
        		}
        	}
        }
        */
//		cout << "\n iw_step = " << iw_step;
//		cout << "\n iw = " << iw;
    }

    return;
}

void SpillageStep::init_QS_matrix( const int &istr)
{
    if (test==1)TITLE("SpillageStep","init_QS_matrix");

    if (USEPW)
    {
        // notice! calculate_psi1d must consistent
        // with psi3d in the same init_QS_matrix
        // function. because the init_QS_matrix
        // may be called in two differnt places.
        PW.calculate_psi1d();
    }

    for (int ik = 0; ik < this->data[istr].nks; ik++)
    {
        if (USEPW)
        {
            PW.calculate_psi3d(ilevel, ik);
        }

        this->data[istr].Soverlap[ik].zero_out();
        this->data[istr].Qoverlap[ik].zero_out();
        // init Q matrix
        for (int ib = 0; ib < this->data[istr].nbands; ib++)
        {
            for (int iw = 0; iw < this->data[istr].nwfc2; iw++)
            {
                const int jw = this->wayd[iw].iw00;
                const int T = this->wayd[iw].type;
                const int L = this->wayd[iw].L;
                const int N = this->wayd[iw].N;
                for (int ie=0; ie<this->data[istr].ne; ie++)
                {
                    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    // if use 'full orbitals'
                    if ( this->wayd[iw].average == false)
                    {
                        data[istr].Qoverlap[ik](iw,ib) +=
                            input.Coef.C4( T, L, N, ie ) *
                            data[istr].Qin(ik, ib, jw, ie);
                    }
                    else// if use 'average orbitals'
                    {
                        const int nofm = 2*L+1; // number of m
                        complex<double> avalue = complex<double>(0.0, 0.0);

                        for (int m=0; m< nofm; m++)
                        {
                            avalue += data[istr].Qin(ik, ib, jw+m, ie);
                        }
                        avalue /= nofm;

                        data[istr].Qoverlap[ik](iw,ib) +=
                            input.Coef.C4( T, L, N, ie ) * avalue;
                    }
                    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                }
            }
        }

        // init S matrix
        for (int iw=0; iw<nwfc2; iw++)
        {
            const int T1 = this->wayd[iw].type;
            const int L1 = this->wayd[iw].L;
            const int N1 = this->wayd[iw].N;
            const int iw001 = this->wayd[iw].iw00;

            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            for (int iw2 = iw; iw2 < nwfc2; iw2++)
            {
                const int T2 = this->wayd[iw2].type;
                const int L2 = this->wayd[iw2].L;
                const int N2 = this->wayd[iw2].N;
                const int iw002 = this->wayd[iw2].iw00;

                if (USEPW)
                {
                    this->data[istr].Soverlap[ik](iw, iw2) = PW.calculateS(iw, iw2, ik);
#ifdef __MPI
                    Parallel_Reduce::reduce_complex_double_pool(this->data[istr].Soverlap[ik](iw, iw2));
#endif
                }
                else
                {
//				double value = 0.0;
                    for (int ie = 0; ie < ne; ie++)
                    {
                        for (int ie2 = 0; ie2 < ne; ie2++)
                        {
                            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            // case 1
                            if ( !wayd[iw].average && !wayd[iw2].average)
                            {
                                this->data[istr].Soverlap[ik](iw, iw2) +=
                                    input.Coef.C4( T1, L1, N1, ie ) *
                                    input.Coef.C4( T2, L2, N2, ie2) *
                                    input.QS_data[istr].Sq1q2[ik](iw001, iw002, ie, ie2);
                            }
                            // case 2
                            else
                            {
                                complex<double> avalue = complex<double>(0.0, 0.0);
                                // case 21
                                if ( wayd[iw].average && wayd[iw2].average)
                                {
                                    for (int m1=0; m1<2*L1+1; m1++)
                                    {
                                        for (int m2=0; m2<2*L2+1; m2++)
                                        {
                                            avalue += input.QS_data[istr].Sq1q2[ik](iw001+m1, iw002+m2, ie, ie2);
                                        }
                                    }
                                    avalue /= ((2*L1+1) * (2*L2+1));
                                }
                                // case 22
                                if ( !wayd[iw].average && wayd[iw2].average)
                                {
                                    for (int m2=0; m2<2*L2+1; m2++)
                                    {
                                        avalue += input.QS_data[istr].Sq1q2[ik](iw001, iw002+m2, ie, ie2);
                                    }
                                    avalue /= (2*L2+1);
                                }
                                // case 23
                                if ( wayd[iw].average && !wayd[iw2].average)
                                {
                                    for (int m1=0; m1<2*L1+1; m1++)
                                    {
                                        avalue += input.QS_data[istr].Sq1q2[ik](iw001+m1, iw002, ie, ie2);
                                    }
                                    avalue /= (2*L1+1);
                                }
                                this->data[istr].Soverlap[ik](iw, iw2) +=
                                    input.Coef.C4( T1, L1, N1, ie ) *
                                    input.Coef.C4( T2, L2, N2, ie2) *
                                    avalue;
                                //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            }
                        }// ie2
                    }// ie
                }// iw

//				if(iw==3)
//				{
//					ofs_running << "\n iw=" << iw << " iw2=" << iw2 << " S[" << ik << "]=" << this->data[istr].Soverlap[ik](iw, iw2);
//				}
            }// ib
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        }//ik


//		if(ik==0)
//		{
//			PRINTCM("S", data[istr].Soverlap[ik]);
//			PRINTCM("Q", data[istr].Qoverlap[ik]);
//		}

        // because S has symmetry : S^{+} = S
//        for(int iw=0; iw<nwfc2; iw++)
//        {
//           for(int iw2=0; iw2<iw; iw2++)
        //          {
        //              Soverlap[ik](iw, iw2) = conj( Soverlap[ik](iw2, iw) );
        //          }
        //     }


        //====================
        // get inverse matrix
		// Peize Lin update 2015-12-05
        //====================
		this->data[istr].inverse_S[ik].reset_loop();
		this->data[istr].inverse_S[ik].cal_inverse(this->data[istr].Soverlap[ik]);
        // notice : only half A is stored!!!
        this->data[istr].Sinv[ik] = this->data[istr].inverse_S[ik].get_inverse();
        //if(ik==0)
        //PRINTCM("init_QS_matrix :: Sinv",data[istr].inverse.A);
    }

    return;
}

// (1) Calculate The Spillage Value
// (2) For each Coefficient(type,l,n,ie) changed orbital.
double SpillageStep::get_spillage(
    const int &istr, // index of structure
    const int &ic, // index of orbitals: (type,l,n)
    const int &ie // index of eigenvalue of Jlq
)
{
    double spillage; // not need to initialized.
    double fill = 0.0;
    Step_Data* Stru = &this->data[istr];

    // mohan add 2010-05-02
    // (1) this part is to control the number of
    // bands needed to be performed.
    static int bands_start;
    static int bands_end;
    static int bands_number;

    if (BANDS_CONTROL)
    {
        bands_start = BANDS_START;
        bands_end = BANDS_END;
        assert(BANDS_END <= Stru->nbands);
        bands_number = bands_end - bands_start;
        assert(bands_number <= Stru->nbands);
    }
    else
    {
        bands_start = 0;
        bands_end = Stru->nbands;
        bands_number = bands_end - bands_start;
    }

    // (2) get new c4 coefficient from ic.
    const double c4_new = input.Coef.C4(
                              this->wayc[ic].type,
                              this->wayc[ic].L,
                              this->wayc[ic].N,
                              ie);

    // (3) get old c4 coefficient.
    const double c4_old = input.Coef.C4_old(
                              this->wayc[ic].type,
                              this->wayc[ic].L,
                              this->wayc[ic].N,
                              ie);

    int unstable = 0;

    if (USEPW)
    {
        PW.nwfc2 = nwfc2;
        PW.update_psi1d(ilevel, ic, ie, c4_new, c4_old);
    }

    //ofs_running << "\n ic_now = " << ic << " ie_now=" << ie;
	
	#ifdef _OPENMP			// Peize Lin add 2016-01-18	
	omp_lock_t unstable_lock;
	omp_init_lock(&unstable_lock);
	#endif
	
	//#ifdef _OPENMP
	//omp_set_num_threads(8);		// Peize Lin test
	//#endif

	#pragma omp parallel for			// Peize Lin add 2016-01-18	
    for (int ik = 0; ik < Stru->nks ; ik++)
    {		
        //	cout << "ik=" << ik << endl;
        //	cout << "\n ik = " << ik <<" c4_new=" << c4_new << " c4_old=" << c4_old << " nks=" << Stru->nks << endl;

        // If new Q,S corespodding to the new C4
        // is not accepted.
        Stru->Qtrial[ik] = Stru->Qoverlap[ik];
        Stru->Strial[ik] = Stru->Soverlap[ik];

        this->newQ(istr, ic, ie, ik, c4_new, c4_old);

        if (USEPW)
        {
            PW.update_psi3d(ilevel, ic, ik);

            timer::tick("PW_Basis","calculateS");

            // old version need to change pw_basis too
            if (NPROC_IN_POOL==1)
            {
                for (int mu=0; mu<nwfc2; mu++)
                {
                    for (int nu=mu; nu<nwfc2; nu++)
                    {
                        const int ic1 = this->wayd[mu ].ic;
                        const int ic2 = this->wayd[nu].ic;
                        if (ic1==ic || ic2==ic)
                        {
                            this->data[istr].Strial[ik](mu, nu) = PW.calculateS(mu, nu, ik);
                        }
                    }
                }
            }
            else				
            {
				// may error when openmp? unchecked		// Peize Lin note 2016-01-18
				#ifdef _OPENMP	
				if(0==omp_get_thread_num())
					cout<<"src_spillage/SpillageStep.cpp line1103 openmp "<<endl;	
				exit(0);
				#endif
				
                int count = 0;
                for (int mu=0; mu<nwfc2; mu++)
                {
                    for (int nu=mu; nu<nwfc2; nu++)
                    {
                        const int ic1 = this->wayd[mu].ic;
                        const int ic2 = this->wayd[nu].ic;

                        // a small trick
                        if (ic1==ic || ic2==ic)
                        {
                            PW.save[count] = PW.calculateS(mu, nu, ik);
                            PW.posmu[count] = mu;
                            PW.posnu[count] = nu;
                            //this->data[istr].Strial[ik](mu, nu) = PW.calculateS(mu, nu, ik);
                            //ofs_running << "\n save=" << save[count] << " mu=" << mu << " nu=" << nu;
                            ++count;
                        }
                    }
                }
#ifdef __MPI
                Parallel_Reduce::reduce_complex_double_pool(PW.save, count);
#endif
                for (int index=0; index<count; index++)
                {
                    this->data[istr].Strial[ik](PW.posmu[index], PW.posnu[index]) = PW.save[index];
                    //	ofs_running << "\n mu=" << pos[index*2] << " nu=" << pos[index*2+1] << " save=" << save[index];
                }
            }
            timer::tick("PW_Basis","calculateS");
        }
        else
        {
            this->newS(istr, ic, ie, ik, c4_new, c4_old);
        }

		Stru->inverse_S[ik].cal_inverse(Stru->Strial[ik]);		// Peize Lin update 2015-12-05
		
//		if(ik==3)
//		{
//	  		PRINTCM("S",data[istr].Strial[ik]);
//	  		PRINTCM("Q",data[istr].Qtrial[ik]);
//	  		PRINTCM("Sinv",Stru->inverse.A);
//		}

//		BLOCK_HERE("spillage");

        // assert all 'Stru' have same k points and bands.
        // mohan 2010-05-02
        Stru->Mk[ik] = 0.00;
//>>>>>
// Keep the flexibility to calculate the spillage value
// of different number of bands for different structures
        //for (int ib = 0; ib < Stru->nbands; ib++)
//<<<<<
        //but now I want to calculate the spillage of given number of bands
        //(1)from input, for all steps. (2) for each step.
        for (int ib=bands_start; ib<bands_end; ib++)
        {
            Stru->Mkb_used(ik,ib) = 0.0;
            for (int iw = 0; iw < Stru->nwfc2; iw++)
            {
                for (int iw2 = 0; iw2 < Stru->nwfc2; iw2++)
                {
                    //Mk += (conj(Qtrial[ik](iw, ib)) * sinv(iw, iw2) * Qtrial[ik](iw2, ib)).real();
                    if (iw <= iw2)
                    {
                        Stru->Mkb_used(ik,ib) += (conj( Stru->Qtrial[ik](iw, ib))
                                                  * Stru->inverse_S[ik].get_inverse()(iw, iw2)			// Peize Lin update 2015-12-05
                                                  * Stru->Qtrial[ik](iw2, ib)).real();
                    }
                    else
                    {
                        Stru->Mkb_used(ik,ib) += (conj( Stru->Qtrial[ik](iw, ib))
                                                  * conj( Stru->inverse_S[ik].get_inverse()(iw2, iw))	// Peize Lin update 2015-12-05
                                                  * Stru->Qtrial[ik](iw2, ib)).real();
                    }
                }
            }
            Stru->Mk[ik] += Stru->Mkb_used(ik,ib);
            //ofs_running << "ik=" << ik << " ib=" << ib << " Mk=" << Stru->Mk[ik] << endl;
        }

        Stru->Mk[ik] /= (double)bands_number;
        //cout << endl;

        // in fact this is caused by numerical unstability.
		#ifdef _OPENMP
		omp_set_lock(&unstable_lock);					// Peize Lin add 2016-01-18	
		#endif
        //if (Stru->Mk[ik]*Stru->weight[ik] > upbound)//mohan fix bug 2010-06-14
        if (fill > Stru->spillage0 )//mohan refix bug 2010-06-18
        {
            ++unstable;
        }
        else if (Stru->Mk[ik] < 0.0)
        {
            ++unstable;
        }
        else
        {
            fill += Stru->Mk[ik] * Stru->weight[ik];
        }
		#ifdef _OPENMP
		omp_unset_lock(&unstable_lock);					// Peize Lin add 2016-01-18	
		#endif		
    }//end ik

#ifdef __MPI
//   Parallel_Reduce::reduce_int_all(unstable);
    if (unstable>0)
    {
        ofs_running << "\n unstable=" << unstable;
        //return  Stru->spillage0;
    }
    Parallel_Reduce::reduce_double_allpool(fill);
#endif

    spillage = Stru->spillage0 - fill;
	
    if (spillage < 0.0 || spillage > 1.0 || unstable)
    {

// Peize Lin test			
for (int ik = 0; ik < Stru->nks ; ik++)
{
	//cout<<"ik:\t"<<ik<<endl;
	cout_matrix( Stru->Strial[ik]);
	cout_matrix( Stru->inverse_S[ik].get_inverse());
	for( size_t i=0; i<Stru->Strial[ik].nr; ++i)
		for( size_t j=0; j<i; ++j)
			Stru->Strial[ik](i,j) = conj(Stru->Strial[ik](j,i));
	ComplexMatrix inverse_tmp(Stru->inverse_S[ik].get_inverse());
	for( size_t i=0; i<inverse_tmp.nr; ++i)
		for( size_t j=0; j<i; ++j)
			inverse_tmp(i,j) = conj(inverse_tmp(j,i));	
	cout_matrix( Stru->Strial[ik] * inverse_tmp);
}

		cout << "\n spillage0 = " << Stru->spillage0;
        cout << "\n Type=" << this->wayc[ic].type
             << " L=" << this->wayc[ic].L
             << " N=" << this->wayc[ic].N;
        cout << "\n bands_number = " << bands_number;

        for (int ib=0; ib<bands_number; ib++)
        {
            double mm = 0.0;
            for (int ik=0; ik<NKS; ik++)
            {
                mm += Stru->Mkb_used(ik,ib) * Stru->weight[ik]/bands_number;
            }
            cout << "\n ib=" << ib << " new fill=" << mm;
        }

        double sum = 0.0;
        for (int ik=0; ik<NKS; ik++)
        {
            cout << "\n M[" << ik << "]=" << Stru->Mk[ik] * Stru->weight[ik] << " " << Stru->Mk[ik];
            sum += Stru->Mk[ik] * Stru->weight[ik];
        }
        cout << "\n sum fill = " << sum;
        cout << endl;
        QUIT();
        return Stru->spillage0;
    }

    return spillage;
}


void SpillageStep::newQ( const int &istr, const int &ic, const int &ie, const int &ik,
                         const double &c4_now ,const double &c4_old)
{
//	TITLE("SpillageStep","newQ");
    const double delta = c4_now - c4_old;
    for (int iw = 0; iw < nwfc2; iw++)
    {
        const int jw = this->wayd[iw].iw00;
        if ( this->wayd[iw].ic == ic)
        {
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if ( !wayd[iw].average )
            {
                // the change of dQ=<d Jlq(iw)|psi(k,b)> is due to
                // the change of Jlq(iw)=delta(Jlq(iw))
                for (int ib = 0; ib < this->data[istr].nbands; ib++)
                {
                    this->data[istr].Qtrial[ik](iw, ib) += delta
                                                         * this->data[istr].Qin(ik, ib, jw, ie);
                }
            }
            else // if use average orbital
            {
                const int L = wayd[iw].L;
                for (int ib = 0; ib < this->data[istr].nbands; ib++)
                {
                    complex<double> avalue = complex<double>(0.0, 0.0);
                    for (int m=0; m<2*L+1; m++)
                    {
                        avalue += this->data[istr].Qin(ik, ib, jw+m, ie);
                    }
                    avalue /= (2*L+1);
                    this->data[istr].Qtrial[ik](iw, ib) += delta * avalue;
                }
            }
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        }
    }
    return;
}


void SpillageStep::newS(
    const int &istr, const int &ic, const int &ie, const int &ik,
    const double &c4_now ,const double &c4_old)
{
//	TITLE("SpillageStep","newS");
    // example c = c2, ie2 = ie
    //==========================================================================================================
    //              iw0-s           iw1p1       iw2p2       iw3p3   iw4s
    //              c1              *c2         *c2         *c2     c3
    //              ea eb ec        ea eb ec    ea eb ec ea eb ec   ea eb ec
    //          ea                      h           h           h
    // iw0-s    eb                      h           h           h
    // c1       ec                      h           h           h
    //----------------------------------------------------------------------------------------------------------
    //          ea                      h           h           h
    // iw1-p1   *eb h  h  h          ?  h ?      ?  h ?         h   h  h   h    (? must do because JLq contain k vector)
    // *c2      ec                      h           h           h
    //----------------------------------------------------------------------------------------------------------
    //          ea                      h           h           h
    // iw2-p2   *eb h  h  h          ?  h ?         h           h   h  h   h
    // *c2      ec                      h           h           h
    //----------------------------------------------------------------------------------------------------------
    //          ea                      h           h           h
    // iw3-p3   *eb h  h  h             h           h           h   h  h   h
    // *c2      ec                      h           h           h
    //-----------------------------------------------------------------------------------------------------------
    //          ea                      h           h           h
    // iw4-s    eb                      h           h           h
    // c3       ec                      h           h           h
    //------------------------------------------------------------------------------------------------------------
    //          ea                      h           h           h
    // iw5-p1   eb                      .           .           .
    // c4       ec                      .           .           .
    //-----------------------------------------------------------------------------------------------------------
    //          ea
    // iw6-p2   eb
    // c4       ec
    //--------------------------------------------------------------------------
    //          ea
    // iw7-p3   eb
    // c4       ec
    //
    //          (e1)
    //==========================================================================
    //=================================================
    // search in all wavefunction ( type, i, l, m, n )
    //=================================================
    for (int iw = 0; iw < this->data[istr].nwfc2; iw++)
    {
        const int jw = this->wayd[iw].iw00;
        const int ic1 = this->wayd[iw].ic;
        const int T1 = this->wayd[iw].type;
        const int L1 = this->wayd[iw].L;
        const int N1 = this->wayd[iw].N;
        // belong to which group of c, eg. p1,p2,p3 may
        // have the same ic1.
        //=============================================
        // search in all other eigenvalues (ic1 != ic).
        //=============================================
        for (int ie1 = 0; ie1 < this->data[istr].ne; ie1++)
        {
            // update for the 'ic' row
            if ( ie1 != ie || ic1 !=ic )
            {
                for (int iw2 = iw; iw2 < this->data[istr].nwfc2; iw2++)
                {
                    const int kw = this->wayd[iw2].iw00;
                    if (this->wayd[iw2].ic == ic)
                    {
                        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        // case 1
                        if ( !wayd[iw].average && !wayd[iw2].average)
                        {
                            this->data[istr].Strial[ik](iw, iw2) +=
                                (c4_now * input.Coef.C4( T1, L1, N1, ie1 ) -
                                 c4_old * input.Coef.C4_old( T1, L1, N1, ie1 ) )
                                * input.QS_data[istr].Sq1q2[ik](jw, kw, ie1, ie);
                        }
                        else // case 2
                        {
                            const int L2 = this->wayd[iw2].L;
                            complex<double> avalue = complex<double>(0.0, 0.0);
                            // case 2.1
                            if ( wayd[iw].average && !wayd[iw2].average)
                            {
                                for (int m1=0; m1<2*L1+1; m1++)
                                {
                                    avalue += input.QS_data[istr].Sq1q2[ik](jw+m1, kw, ie1, ie);
                                }
                                avalue /= (2*L1+1);
                            }
                            // case 2.2
                            if ( !wayd[iw].average && wayd[iw2].average)
                            {
                                for (int m2=0; m2<2*L2+1; m2++)
                                {
                                    avalue += input.QS_data[istr].Sq1q2[ik](jw, kw+m2, ie1, ie);
                                }
                                avalue /= (2*L2+1);
                            }
                            // case 2.3
                            if ( wayd[iw].average && wayd[iw2].average)
                            {
                                for (int m1=0; m1<2*L1+1; m1++)
                                {
                                    for (int m2=0; m2<2*L2+1; m2++)
                                    {
                                        avalue += input.QS_data[istr].Sq1q2[ik](jw+m1, kw+m2, ie1, ie);
                                    }
                                }
                                avalue /= ( (2*L1+1)*(2*L2+1) );
                            }
                            this->data[istr].Strial[ik](iw, iw2) +=
                                (c4_now * input.Coef.C4( T1, L1, N1, ie1 ) -
                                 c4_old * input.Coef.C4_old( T1, L1, N1, ie1 ) )
                                * avalue;
                        }
                        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    }
                }
            }
            else
            {
                for (int iw2 = iw; iw2 < this->data[istr].nwfc2; iw2++)
                {
                    const int ic2 = this->wayd[iw2].ic;
                    const int kw = this->wayd[iw2].iw00;
                    const int T2 = this->wayd[iw2].type;
                    const int L2 = this->wayd[iw2].L;
                    const int N2 = this->wayd[iw2].N;
                    for (int ie2 = 0; ie2 < this->data[istr].ne; ie2++)
                    {
                        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        if ( !wayd[iw].average && !wayd[iw2].average)
                        {
                            this->data[istr].Strial[ik](iw, iw2) +=
                                (c4_now * input.Coef.C4( T2, L2, N2, ie2 ) -
                                 c4_old * input.Coef.C4_old( T2, L2, N2, ie2 ) )
                                * input.QS_data[istr].Sq1q2[ik](jw, kw, ie, ie2);
                        }
                        else // case 2
                        {
                            const int L2 = this->wayd[iw2].L;
                            complex<double> avalue = complex<double>(0.0, 0.0);
                            // case 2.1
                            if ( wayd[iw].average && !wayd[iw2].average)
                            {
                                for (int m1=0; m1<2*L1+1; m1++)
                                {
                                    avalue += input.QS_data[istr].Sq1q2[ik](jw+m1, kw, ie, ie2);
                                }
                                avalue /= (2*L1+1);
                            }
                            // case 2.2
                            if ( !wayd[iw].average && wayd[iw2].average)
                            {
                                for (int m2=0; m2<2*L2+1; m2++)
                                {
                                    avalue += input.QS_data[istr].Sq1q2[ik](jw, kw+m2, ie, ie2);
                                }
                                avalue /= (2*L2+1);
                            }
                            // case 2.3
                            if ( wayd[iw].average && wayd[iw2].average)
                            {
                                for (int m1=0; m1<2*L1+1; m1++)
                                {
                                    for (int m2=0; m2<2*L2+1; m2++)
                                    {
                                        avalue += input.QS_data[istr].Sq1q2[ik](jw+m1, kw+m2, ie, ie2);
                                    }
                                }
                                avalue /= ( (2*L1+1)*(2*L2+1) );
                            }
                            this->data[istr].Strial[ik](iw, iw2) +=
                                (c4_now * input.Coef.C4( T2, L2, N2, ie2 ) -
                                 c4_old * input.Coef.C4_old( T2, L2, N2, ie2 ) )
                                * avalue;
                        }
                        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    }
                }
            }
        }
    }

    //  PRINTCM("S",Strial[ik]);
    return;
}



void SpillageStep::updateQS( const int &istr )
{
//    if(test==1)TITLE("SpillageStep","updateQS");
    for (int ik=0; ik< data[istr].nks; ik++)
    {
        data[istr].Qoverlap[ik] = data[istr].Qtrial[ik];
        data[istr].Soverlap[ik] = data[istr].Strial[ik];
		data[istr].inverse_S[ik].update();				// Peize Lin update 2015-12-05
	}
    return;
}


#include "numerical_basis.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry.h"
#include "winput.h"
#include "../src_global/math_ylmreal.h"

bool Numerical_Basis::init_label = false;
Bessel_Basis Numerical_Basis::bessel_basis;
IntArray* Numerical_Basis::mu_index;

Numerical_Basis::Numerical_Basis() {}
Numerical_Basis::~Numerical_Basis() {}

//============================================================
// MEMBER FUNCTION :
// NAME : init
// DESCRIPTION : Two main functions:
// (1) start_from_file = true;
// Firstly, use check(1) to call bessel_basis.init
// to generate TableOne.
// Secondly readin C4 from file.
// Thirdly generate 3D atomic wfc in G space, put the
// results in psi.
//
// (2) If output overlap Q, start_from_file = false;
// Firstly, use check(0) to call bessel_basis,init
// to generate TableOne
// Secondly output overlap, use psi(evc) and jlq3d.
//============================================================
void Numerical_Basis::start_from_file_k( const int &ik, ComplexMatrix &psi)
{
    TITLE("Numerical_Basis","start_from_file_k");

    if (!Numerical_Basis::init_label)
    {
        // 1 stands for : start_from_file
        Numerical_Basis::bessel_basis.init( 1, pw.ecutwfc, ucell.ntype, ucell.lmax );
        Numerical_Basis::init_mu_index();
        Numerical_Basis::init_label = true;
    }
    this->numerical_atomic_wfc(ik, kv.ngk[ik], psi);
    return;
}


// The function is called in run_fp.cpp.
void Numerical_Basis::output_overlap( const ComplexMatrix *psi)
{
    TITLE("Numerical_Basis","output_overlap");
    NEW_PART("Overlap Data For Spillage Minimization");
    
	//---------------------------------------------------------
	// if the numerical_basis hasn't been initialized yet,
    // then we initial here.
	//---------------------------------------------------------
    if (!Numerical_Basis::init_label)
    {
        // 0 stands for : 'Faln' is not used.
        Numerical_Basis::bessel_basis.init( 0, pw.ecutwfc, ucell.ntype, ucell.lmax );
        Numerical_Basis::init_mu_index();
        Numerical_Basis::init_label = true;
    }

    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/" << ucell.latName << "." << ucell.lat0 << ".dat";
    if (MY_RANK==0)
    {
        ofs.open(ss.str().c_str());
    }

    const int nks = kv.nks;
    const int ne = Numerical_Basis::bessel_basis.get_ecut_number();

	// OVERLAP : < J_mu | Psi >
    realArray overlap_Q1(nks, NBANDS, NLOCAL, ne );
    realArray overlap_Q2(nks, NBANDS, NLOCAL, ne );

	// OVERLAP : < J_mu | J_nu >
    realArray *Sq_real = new realArray[nks];
    realArray *Sq_imag = new realArray[nks];

    // (1) allocate Sq matrix.
    if (winput::out_spillage == 2)
    {
        for (int ik=0; ik<nks; ik++) Sq_real[ik].create( NLOCAL, NLOCAL, ne, ne );
        for (int ik=0; ik<nks; ik++) Sq_imag[ik].create( NLOCAL, NLOCAL, ne, ne );
    }

    ZEROS(overlap_Q1.ptr, overlap_Q1.getSize() );
    ZEROS(overlap_Q2.ptr, overlap_Q2.getSize() );

    for (int ik=0; ik<nks; ik++)
    {
        ZEROS(Sq_real[ik].ptr, Sq_real[ik].getSize() );
        ZEROS(Sq_imag[ik].ptr, Sq_imag[ik].getSize() );
    }

	OUT(ofs_running,"number of k points",overlap_Q1.getBound1());
	OUT(ofs_running,"number of bands",overlap_Q1.getBound2());
	OUT(ofs_running,"number of local orbitals",overlap_Q1.getBound3());
	OUT(ofs_running,"number of eigenvalues of Jl(x)",overlap_Q1.getBound4());

    // nks now is the reduced k-points.
    for (int ik=0; ik<nks; ik++)
    {
        const int npw= kv.ngk[ik];
		ofs_running << " --------------------------------------------------------" << endl;
		ofs_running << " Print the overlap matrixs Q and S for this kpoint";
        ofs_running << "\n " << setw(8) << "ik" << setw(8) << "npw";
        ofs_running << "\n " << setw(8) << ik+1 << setw(8) << npw << endl;
		ofs_running << " --------------------------------------------------------" << endl;
        // search for all k-points.
        this->jlq3d_overlap(overlap_Q1, overlap_Q2, ik, ik, npw, psi[ik]);
        DONE(ofs_running,"jlq3d_overlap");

        // (2) generate Sq matrix if necessary.
        if (winput::out_spillage == 2)
        {
            this->Sq_overlap( Sq_real[ik], Sq_imag[ik], ik, npw );
            DONE(ofs_running,"Sq_overlap");
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_double_pool( overlap_Q1.ptr, overlap_Q1.getSize() );
    Parallel_Reduce::reduce_double_pool( overlap_Q2.ptr, overlap_Q2.getSize() );
    for (int ik=0; ik<nks; ik++)
    {
        Parallel_Reduce::reduce_double_pool( Sq_real[ik].ptr, Sq_real[ik].getSize() );
        Parallel_Reduce::reduce_double_pool( Sq_imag[ik].ptr, Sq_imag[ik].getSize() );
    }
#endif

	// only print out to the information by the first processor
    if (MY_RANK==0)
    {
        ofs.precision(10);
        ofs << ucell.lat0 << endl;

        ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        ofs << ucell.ntype << " ntype" << endl;
        for (int it=0; it<ucell.ntype; it++)
        {
            ofs << ucell.atoms[it].label << " label" << endl; // mohan add 2009-07-23
            ofs << ucell.atoms[it].na << " na" << endl;
            for (int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                ofs << ucell.atoms[it].tau[ia].x
                << " " << ucell.atoms[it].tau[ia].y
                << " " << ucell.atoms[it].tau[ia].z << endl;
            }
        }
        // ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
        ofs << pw.ecutwfc << " ecutwfc" << endl; // mohan add 2009-09-08

        // this parameter determine the total number of jlq.
        ofs << Numerical_Basis::bessel_basis.get_ecut() << " ecutwfc_jlq" << endl;//mohan modify 2009-09-08
        ofs << Numerical_Basis::bessel_basis.get_rcut() << " rcut_Jlq" << endl;

        // mohan add 'smooth' and 'sigma' 2009-08-28
        ofs << Numerical_Basis::bessel_basis.get_smooth() << " smooth" << endl;
        ofs << Numerical_Basis::bessel_basis.get_sigma() << " sigma" << endl;

        ofs << Numerical_Basis::bessel_basis.get_tolerence() << " tolerence" << endl;

        ofs << ucell.lmax << " lmax" << endl;
    }

    ofs << scientific;

    ofs << setprecision(8);
    // NOTICE: ofs_warning << "\n The precison may affect the optimize result.";
    
    this->output_overlap_Q( ofs, overlap_Q1, overlap_Q2 );

    if (winput::out_spillage == 2)
    {
        this->output_overlap_Sq(ss.str(), ofs, Sq_real, Sq_imag);
    }

    delete[] Sq_real;
    delete[] Sq_imag;

    if (MY_RANK==0) ofs.close();
    return;
}

void Numerical_Basis::output_overlap_Sq(
    const string &name,
    ofstream &ofs,
    const realArray *Sq_real,
    const realArray *Sq_imag)
{
    if (MY_RANK==0)
    {
        ofs << "\n<OVERLAP_Sq>";
        ofs.close();
    }

    int count = 0;
    for (int ik=0; ik< kv.nkstot; ik++)
    {
        if ( MY_POOL == Pkpoints.whichpool[ik] )
        {
            if ( RANK_IN_POOL == 0)
            {
                ofs.open(name.c_str(), ios::app);
                const int ik_now = ik - Pkpoints.startk_pool[MY_POOL];
                for (int i=0; i< Sq_real[ik_now].getSize(); i++)
                {
                    if (count%2==0) ofs << "\n";
                    ofs << " " << Sq_real[ik_now].ptr[i] << " " << Sq_imag[ik_now].ptr[i];
                    ++count;
                }

                ofs.close();
            }
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
        else
        {
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }

        /*
        if(MY_RANK==0)
        for(int i=0; i< Sq_real[ik].getSize(); i++)
        {
        	if(i%2==0) ofs << "\n";
        	ofs << " " << Sq_real[ik].ptr[i] << " " << Sq_imag[ik].ptr[i];
        }
        */

    }
    if (MY_RANK==0)
    {
        ofs.open(name.c_str(), ios::app);
        ofs << "\n</OVERLAP_Sq>" << endl;
    }
    return;
}

void Numerical_Basis::output_overlap_Q(
    ofstream &ofs,
    const realArray &overlap_Q1,
    const realArray &overlap_Q2)
{
    // (1)
    if (MY_RANK==0)
    {
        ofs << kv.nkstot << " nks" << endl;
        ofs	<< overlap_Q1.getBound2() << " nbands" << endl;
        ofs	<< overlap_Q1.getBound3() << " nwfc" << endl;
        ofs	<< overlap_Q1.getBound4() << " ne";
        ofs << "\n<WEIGHT_OF_KPOINTS>";
    }

    // (2)
    for (int ik=0; ik<kv.nkstot; ik++)
    {
        double kx, ky, kz, wknow;
#ifdef __MPI
        const int pool = Pkpoints.whichpool[ik];
        const int iknow = ik - Pkpoints.startk_pool[MY_POOL];
        if (RANK_IN_POOL==0)
        {
            if (MY_POOL==0)
            {
                if (pool==0)
                {
                    kx = kv.kvec_c[ik].x;
                    ky = kv.kvec_c[ik].y;
                    kz = kv.kvec_c[ik].z;
                    wknow = kv.wk[ik];
                }
                else
                {
                    MPI_Status ierror;
                    MPI_Recv(&kx, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&ky, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+1, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&kz, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+2, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&wknow, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+3, MPI_COMM_WORLD,&ierror);
                }
            }
            else
            {
                if (MY_POOL == pool)
                {
                    MPI_Send(&kv.kvec_c[iknow].x, 1, MPI_DOUBLE, 0, ik*4, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].y, 1, MPI_DOUBLE, 0, ik*4+1, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].z, 1, MPI_DOUBLE, 0, ik*4+2, MPI_COMM_WORLD);
                    MPI_Send(&kv.wk[iknow], 1, MPI_DOUBLE, 0, ik*4+3, MPI_COMM_WORLD);
                }
            }
        }
        // this barrier is very important
        MPI_Barrier(MPI_COMM_WORLD);
#else
        if (MY_RANK==0)
        {
            kx = kv.kvec_c[ik].x;
            ky = kv.kvec_c[ik].y;
            kz = kv.kvec_c[ik].z;
            wknow = kv.wk[ik];
        }
#endif

        if (MY_RANK==0)
        {
            ofs << "\n" << kx << " " << ky << " " << kz;
            ofs << " " << wknow * 0.5;
        }
    }

    // (3)
    if (MY_RANK==0)
    {
        ofs << "\n</WEIGHT_OF_KPOINTS>" << endl;
        ofs << "\n<OVERLAP_Q>";
    }

    // (4)
    /*
    if(MY_RANK==0)
    {
    //    	for( int i=0; i<overlap_Q1.getSize(); i++)
    //    	{
    //    		if( i%2==0 ) ofs << "\n";
    //    		ofs << " " << overlap_Q1.ptr[i] << " " << overlap_Q2.ptr[i];
    //    	}
    }
    */

    const int ne = overlap_Q1.getBound4();
    const int dim = NBANDS * NLOCAL * ne;
    double *Qtmp1 = new double[dim];
    double *Qtmp2 = new double[dim];
    int count = 0;

    for (int ik=0; ik<kv.nkstot; ik++)
    {
        ZEROS(Qtmp1, dim);
        ZEROS(Qtmp2, dim);
        Pkpoints.pool_collection(Qtmp1, Qtmp2, overlap_Q1, overlap_Q2, ik);
        if (MY_RANK==0)
        {
    //        ofs << "\n ik=" << ik;
            // begin data writing.
            for (int i=0; i<dim; i++)
            {
                if ( count%4==0 ) ofs << "\n";
                ofs << " " << Qtmp1[i] << " " << Qtmp2[i];
                ++count;
            }
            // end data writing.
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    delete[] Qtmp1;
    delete[] Qtmp2;

    // (5)
    if (MY_RANK==0)
    {
        ofs << "\n</OVERLAP_Q>" << endl;
    }
    return;
}

void Numerical_Basis::Sq_overlap(
    realArray &Sq_real,
    realArray &Sq_imag,
    const int &ik,
    const int &np)
{
    TITLE("Numerical_Basis","Sq_overlap");
    timer::tick("Numerical_Basis","Sq_overlap");

	ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS"  << endl;
	ofs_running << " S = < J_mu,q1 | J_nu,q2 >" << endl; 

	const double normalization = (4 * PI) * (4 * PI) / ucell.omega;			// Peize Lin add normalization 2015-12-29
	
    const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
    matrix ylm(total_lm, np);

    Vector3<double> *gk = new Vector3 <double> [np];
    for (int ig=0; ig<np; ig++)
    {
        gk[ig] = wf.get_1qvec_cartesian(ik, ig);
    }

    YlmReal::Ylm_Real(total_lm, np, gk, ylm);

    const int enumber = Numerical_Basis::bessel_basis.get_ecut_number();

    realArray flq(ucell.lmax+1, enumber, np);

    // get flq(G) = \int f(r)jl(G*r) from interpolation table.
    for (int l=0; l<ucell.lmax+1; l++)
    {
        for (int ie=0; ie<enumber; ie++)
        {
            for (int ig=0; ig<np; ig++)
            {
                flq(l,ie,ig) = Numerical_Basis::bessel_basis.
                               Polynomial_Interpolation2(l, ie, gk[ig].norm() * ucell.tpiba );
            }
        }
    }

    complex<double> *about_ig = new complex<double>[np];
	ZEROS(about_ig, np);

    ofs_running << "\n " << setw(5) << "ik"
    << setw(8) << "Type1"
    << setw(8) << "Atom1"
    << setw(8) << "L1"
    << setw(8) << "Type2"
    << setw(8) << "Atom2"
    << setw(8) << "L2" << endl;

    for (int T1 = 0; T1 < ucell.ntype; T1++) // 1.1
    {
        for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++) // 1.2
        {
            complex<double> *sk = wf.get_sk(ik, T1, I1);
            for (int T2=0; T2<ucell.ntype; T2++) // 2.1
            {
                for (int I2=0; I2<ucell.atoms[T2].na; I2++) // 2.2
                {
                    complex<double> *sk2 = wf.get_sk(ik, T2, I2);
                    for (int l = 0; l < ucell.atoms[T1].nwl+1; l++) // 1.3
                    {
                        complex<double> lphase = normalization * pow(IMAG_UNIT, l);			// Peize Lin add normalization 2015-12-29
                        for (int l2 = 0; l2 < ucell.atoms[T2].nwl+1; l2++) // 2.3
                        {
                            ofs_running << " " << setw(5) << ik+1
                            << setw(8) << ucell.atoms[T1].label
                            << setw(8) << I1+1
                            << setw(8) << l
                            << setw(8) << ucell.atoms[T2].label
                            << setw(8) << I2+1
                            << setw(8) << l2
                            << endl;

                            complex<double> lphase2 = pow(IMAG_UNIT, l2);
                            for (int ic=0; ic < ucell.nmax; ic++) // 1.5
                            {
                                for (int ic2=0; ic2 < ucell.nmax; ic2++) // 2.5
                                {
                                    for (int m=0; m<2*l+1; m++) // 1.6
                                    {
                                        const int lm = l*l+m;
                                        for (int ig=0; ig<np; ig++)
                                        {
                                            about_ig[ig] = conj( lphase * sk[ig] * ylm(lm, ig) );
                                        }
                                        for (int m2=0; m2<2*l2+1; m2++) // 2.6
                                        {
                                            const int lm2 = l2*l2+m2;
                                            const int mu = mu_index[T1](I1, l, ic, m);
                                            const int nu = mu_index[T2](I2,l2,ic2,m2);
                                            for (int ig=0; ig<np; ig++)
                                            {
                                                const complex<double> about_ig3= lphase2 * sk2[ig] * ylm(lm2, ig)
                                                                                 * about_ig[ig];

                                                for (int ie=0; ie < enumber; ie++) // 1.4
                                                {
                                                    for (int ie2=0; ie2 < enumber; ie2++) // 2.4
                                                    {
                                                        const complex<double> s =
                                                            about_ig3 * flq(l,ie,ig) * flq(l2,ie2,ig);
                                                        Sq_real( mu, nu, ie, ie2) += s.real();
                                                        Sq_imag( mu, nu, ie, ie2) += s.imag();

                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    delete[] sk2;
                }
            }
            delete[] sk;
        }
    }

    delete[] about_ig;
	delete[] gk; //mohan fix bug 2011-06-24
    timer::tick("Numerical_Basis","Sq_overlap");
    return;
}

void Numerical_Basis::jlq3d_overlap(
    realArray &overlap_Q1,
    realArray &overlap_Q2,
    const int &ik_ibz,
    const int &ik,
    const int &np,
    const ComplexMatrix &psi)
{
    TITLE("Numerical_Basis","jlq3d_overlap");
    timer::tick("Numerical_Basis","jlq3d_overlap");

	ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS AND BLOCH WAVE FUNCTIONS" << endl;
	ofs_running << " Q = < J_mu, q | Psi_n, k > " << endl;

	const double normalization = (4 * PI) / sqrt(ucell.omega);			// Peize Lin add normalization 2015-12-29

    const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
    matrix ylm(total_lm, np);

    Vector3<double> *gk = new Vector3 <double> [np];
    for (int ig=0; ig<np; ig++)
    {
        gk[ig] = wf.get_1qvec_cartesian(ik, ig);
    }

    YlmReal::Ylm_Real(total_lm, np, gk, ylm);

    ofs_running << "\n " << setw(5) << "ik"
    << setw(8) << "Type1"
    << setw(8) << "Atom1" 
	<< setw(8) << "L"
	<< endl;

    double *flq = new double[np];
    complex<double> overlapQ = ZERO;
    for (int T1 = 0; T1 < ucell.ntype; T1++)
    {
        //OUT("T1",T1);
        for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++)
        {
            //OUT("I1",I1);
            complex<double> *sk = wf.get_sk(ik, T1, I1);
            for (int L=0; L< ucell.atoms[T1].nwl+1; L++)
            {
                ofs_running << " " << setw(5) << ik+1
                            << setw(8) << ucell.atoms[T1].label
                            << setw(8) << I1+1 
							<< setw(8) << L
							<< endl;
                //OUT("l",l);
                complex<double> lphase = normalization * pow(IMAG_UNIT, L);			// Peize Lin add normalization 2015-12-29
                for (int ie=0; ie < Numerical_Basis::bessel_basis.get_ecut_number(); ie++)
                {
                    for (int ig=0; ig<np; ig++)
                    {
                        flq[ig] = Numerical_Basis::bessel_basis.Polynomial_Interpolation2
                                  (L, ie, gk[ig].norm() * ucell.tpiba );
                    }

                    const int N = 0;
                    assert( ucell.nmax == 1);
                    for (int m=0; m<2*L+1; m++)
                    {
                        const int lm = L*L+m;
                        for (int ib=0; ib<NBANDS; ib++)
                        {
                            complex<double> overlap_tmp = ZERO;
                            for (int ig=0; ig<np; ig++)
                            {
                                const complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                                overlap_tmp += conj( local_tmp ) * psi(ib, ig); // psi is bloch orbitals
                            }
                            overlap_Q1(ik_ibz, ib, mu_index[T1](I1, L, N, m), ie) = overlap_tmp.real();
                            overlap_Q2(ik_ibz, ib, mu_index[T1](I1, L, N, m), ie) = overlap_tmp.imag();
                        }
                    }
                }//end ie
            }//end l
            delete[] sk;
        }
    }

    delete[] flq;
    delete[] gk;
    timer::tick("Numerical_Basis","jlq3d_overlap");
    return;
}

void Numerical_Basis::init_mu_index(void)
{
	ofs_running << " Initialize the mu index" << endl;
    Numerical_Basis::mu_index = new IntArray[ucell.ntype];

    int mu = 0;
    for (int it=0; it<ucell.ntype; it++)
    {
        Numerical_Basis::mu_index[it].create(
            ucell.atoms[it].na,
            ucell.atoms[it].nwl+1,
            ucell.nmax,
            2*(ucell.atoms[it].nwl+1)+1); // m ==> 2*l+1

		// mohan added 2021-01-03
		ofs_running << "Type " << it+1 
		<< " number_of_atoms " << ucell.atoms[it].na
		<< " number_of_L " << ucell.atoms[it].nwl+1
		<< " number_of_n " << ucell.nmax
		<< " number_of_m " << 2*(ucell.atoms[it].nwl+1)+1 << endl;

        for (int ia=0; ia<ucell.atoms[it].na; ia++)
        {
            for (int l=0; l< ucell.atoms[it].nwl+1; l++)
            {
                for (int n=0; n< ucell.atoms[it].l_nchi[l]; n++)
                {
                    for (int m=0; m<2*l+1; m++)
                    {
                        Numerical_Basis::mu_index[it](ia,l,n,m) = mu;
                        mu++;
                    }
                }
            }
        }
    }
    return;
}

void Numerical_Basis::numerical_atomic_wfc(
    const int &ik,
    const int &np,
    ComplexMatrix &psi)
{
    TITLE("Numerical_Basis", "numerical_atomic_wfc");

    const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
    matrix ylm(total_lm, np);

    Vector3<double> *gk = new Vector3 <double> [np];
    for (int ig=0; ig<np; ig++)
    {
        gk[ig] = wf.get_1qvec_cartesian(ik, ig);
    }

    YlmReal::Ylm_Real(total_lm, np, gk, ylm);

    int index = 0;
    double *flq = new double[np];
    for (int it = 0; it < ucell.ntype; it++)
    {
        //OUT("it",it);
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            //OUT("ia",ia);
            complex<double> *sk = wf.get_sk(ik, it, ia);
            for (int l = 0; l < ucell.atoms[it].nwl+1; l++)
            {
                //OUT("l",l);
                complex<double> lphase = pow(IMAG_UNIT, l);
                for (int ic=0; ic < ucell.atoms[it].l_nchi[l]; ic++)
                {
                    //OUT("ic",ic);
                    for (int ig=0; ig<np; ig++)
                    {
                        flq[ig] = Numerical_Basis::bessel_basis.
                                  Polynomial_Interpolation(it, l, ic, gk[ig].norm() * ucell.tpiba );
                    }

                    for (int m=0; m<2*l+1; m++)
                    {
                        //OUT("m",m);
                        const int lm = l*l+m;
                        for (int ig=0; ig<np; ig++)
                        {
                            psi( Numerical_Basis::mu_index[it](ia,l,ic,m), ig) =
                                lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                        }
                    }
                }
            }
            delete[] sk;
        }
    }
    delete[] flq;
    delete[] gk;
}

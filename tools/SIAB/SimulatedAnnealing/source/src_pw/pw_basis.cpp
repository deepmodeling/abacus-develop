#include "pw_basis.h"
#include "pw_complement.h"
#include "../src_spillage/tools.h"
#include "ylm_real.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_tools/complexmatrix_inline.h"

PW_Basis::PW_Basis()
{
    // if not parallel, gg_global == gg, gdirect_global == gdirect
    // gcar_global == gcar
    gg_global = new double[1];
    gdirect_global = new Vector3<double>[1];
    gcar_global = new Vector3<double>[1];

#ifdef __MPI
    gg = new double[1];
    gdirect = new Vector3<double>[1];
    gcar = new Vector3<double>[1];
#endif

    ig1 = new int[1];
    ig2 = new int[1];
    ig3 = new int[1];

    psi3d = new ComplexMatrix[1];
    ylm = new matrix[1];
    Jlq = new complex<double>[1];

	save = new complex<double>[1];
	posmu = new int[1];
	posnu = new int[1];
}

PW_Basis::~PW_Basis()
{
    delete [] gcar_global;
    delete [] gdirect_global;
    delete [] gg_global;

#ifdef __MPI
    delete [] gcar;
    delete [] gdirect;
    delete [] gg;
#endif

    delete [] ig1;
    delete [] ig2;
    delete [] ig3;

    delete[] psi3d;
    delete[] ylm;
    delete[] Jlq;

	delete[] save;
	delete[] posmu;
	delete[] posnu;
}

void PW_Basis::init(void)
{
    TITLE("PW_Basis","init");
    timer::tick("PW_Basis","init");

    // (1) setup the energy cutoff, which is different in each pool,
    // because the number of k points is differnt.
    this->setup_gg();

    // (2) setup the FFT dimension according to energy cutoff.
    this->setup_FFT_dimension();

	OUT(ofs_running,"ggwfc2(Ry)",ggwfc2);

    // (3) get the number of plane waves for wave functions.
    this->ngmc_g = PW_complement::get_total_pw_number(0.0, ggchg, nx, ny, nz, GGT);
    this->ngmw_g = PW_complement::get_total_pw_number(0.0, ggwfc2, nx, ny, nz, GGT);

	OUT(ofs_running,"ngmc_g",ngmc_g);
	OUT(ofs_running,"ngmw_g",ngmw_g);

    // (4) get the local number of plane waves.
    assert(NPROC_IN_POOL>0);
    const int remain = ngmw_g % NPROC_IN_POOL;
    this->ngmw = this->ngmw_g / NPROC_IN_POOL;
    if (RANK_IN_POOL < remain)
    {
        this->ngmw++;
    }

    if (RANK_IN_POOL < remain)
    {
        ngmw_start = ngmw * RANK_IN_POOL;
    }
    else
    {
        ngmw_start = (ngmw+1) * remain + ngmw * (RANK_IN_POOL-remain);
    }
	
	OUT(ofs_running,"ngmw",ngmw);
    //ofs_running << "\n ngmw_start = " << ngmw_start;

    // (5) init |g|, g_direct, g_cartesian
    delete[] gdirect_global;
    delete[] gcar_global;
    delete[] gg_global;

    gdirect_global = new Vector3<double>[ngmc_g];// indices of G vectors
    gcar_global = new Vector3<double>[ngmc_g];
    gg_global = new double[ngmc_g];// store the |G|^2 of the 1d array

    PW_complement::get_total_pw(gg_global, gdirect_global, 0.0, ggchg, nx, ny, nz, this->GGT, ngmc_g);
    //PW_complement::get_total_pw(gg_global, gdirect_global, 0.0, ggwfc2, nx, ny, nz, this->GGT, ngmw_g);
    PW_complement::setup_GVectors(this->G, ngmc_g, gg_global, gdirect_global, gcar_global);

    // (6) get ig1, ig2, ig3
#ifdef __MPI
    this->get_MPI_GVectors();
#else
    this->get_GVectors();
#endif
	// mohan add on 2011-07-23
	stringstream ss;
	ss << "ORBITAL_PW_GK" << MY_RANK+1 << ".dat";
	ofstream ofs(ss.str().c_str());
	ofs << ggwfc2 << " (ggwfc2, Ry)" << endl; 
	ofs << ngmw << " (Number of plane wave basis for wave functions)" << endl;
	for(int ig=0; ig<ngmw; ++ig)
	{
		ofs << gcar[ig].x << " " << gcar[ig].y << " " << gcar[ig].z << endl;
	}
	ofs.close(); 

    this->set_igk();

    this->setup_structure_factor();
    timer::tick("PW_Basis","init");
    return;
}

void PW_Basis::setup_gg(void)
{
    TITLE("PW_Basis","setup_gg");

    //(1.1) lattice parameters
    this->GT = LATVEC.Inverse();
    this->G = GT.Transpose();
    this->GGT = G * GT;

    // (1.2) tpiba and ecut
    assert(LAT0 > 0.0);
    this->tpiba = TWO_PI/LAT0;
    this->tpiba2 = tpiba * tpiba;
    this->ggpsi = ECUT / tpiba2;

	double wfac = 4.0;
	this->ggchg = wfac * this->ggpsi;

    // (1.3) get ggwfc2
	assert( NKS > 0);
   	assert( CARKX != NULL );
   	assert( CARKY != NULL );
   	assert( CARKZ != NULL );

    for (int ik=0; ik<NKS; ik++)
    {
        const double k_mod = sqrt(CARKX[ik] * CARKX[ik] + CARKY[ik] * CARKY[ik] + CARKZ[ik] * CARKZ[ik]);
        const double tmp = sqrt(this->ggpsi) + k_mod;
        const double tmp2 = tmp * tmp ;
        if (this->ggwfc2 < tmp2) this->ggwfc2 = tmp2;
    }

	OUT(ofs_running,"NKS",NKS);
	OUT(ofs_running,"ggpsi(Ry)",ggpsi*tpiba2);
	OUT(ofs_running,"ggwfc2(Ry)",ggwfc2*tpiba2);

    // (1.4) initialize itia2iat
    this->itia2iat = new int*[NTYPE];
    int iat=0;
    for (int it=0; it<NTYPE; it++)
    {
        itia2iat[it] = new int[NA[it]];
        for (int ia=0; ia<NA[it]; ia++)
        {
            itia2iat[it][ia] = iat;
            ++iat;
        }
    }

    return;
}


void PW_Basis::set_igk(void)
{
    TITLE("PW_Basis","set_igk");

    // (1) set ngk and npwx
    this->npwx = 0;
    this->ngk = new int[NKS];
    ZEROS(this->ngk, NKS);

    assert(this->ngmw>0);

//    ofs_running << "\n\n ngmw = " << ngmw;
//    ofs_running << "\n gg[ngmw-1]=" << gg[ngmw-1];
    for (int ik=0; ik<NKS; ik++)
    {
        Vector3<double> kcar(CARKX[ik], CARKY[ik], CARKZ[ik]);
        const double k2 = kcar * kcar;
        int ng=0;
        for (int ig=0; ig<ngmw; ig++)
        {
            Vector3<double> f = this->gcar[ig] + kcar;
            const double gk2 = f * f;
            if (sqrt(this->gg[ig]) > sqrt(this->ggpsi) + sqrt(k2))
            {
                break;
            }
            if (gk2 <= this->ggpsi)
            {
                ++ng;
            }
        }
        this->ngk[ik] = ng;
        //ofs_running << "\n ngk[" << ik << "]=" << ngk[ik];
        if ( npwx < ng)
        {
            npwx = ng;
        }
    }

	OUT(ofs_running,"npwx",npwx);

    // set igk
    this->igk = new int*[NKS];
    for (int ik=0; ik<NKS; ik++)
    {
        igk[ik] = new int[this->ngk[ik]];
        ZEROS(igk[ik], ngk[ik]);
    }

    for (int ik=0; ik<NKS; ik++)
    {
        Vector3<double> kcar(CARKX[ik], CARKY[ik], CARKZ[ik]);
        int ng=0;
        const double k2 = kcar * kcar;
        for (int ig=0; ig<ngmw; ig++)
        {
            Vector3<double> f = this->gcar[ig] + kcar;
            const double gk2 = f * f;
            if (sqrt(this->gg[ig]) > sqrt(this->ggpsi) + sqrt(k2))
            {
                break;
            }
            if (gk2 <= this->ggpsi)
            {
                this->igk[ik][ng] = ig;
                ++ng;
            }
        }
    }

    return;
}

//  First stage Basis initialization.
//  Set up crystal structure parameters.
void PW_Basis::setup_FFT_dimension()
{
    TITLE("PW_Basis","setup_FFT_dimension");

    PW_complement::get_FFT_dimension(LATVEC, this->ggchg, nx, ny, nz);
    this->nxyz = nx * ny * nz;

	OUT(ofs_running,"FFT wavefunctions",nx,ny,nz);

    return;
}



#ifdef __MPI
void PW_Basis::get_MPI_GVectors(void)
{
    TITLE("PW_Basis","get_MPI_GVectors");

    assert(ngmw>0);

    // (1) first part
    delete[] gdirect;
    delete[] gcar;
    delete[] gg;
    this->gdirect = new Vector3<double>[ngmw];
    this->gcar  = new Vector3<double>[ngmw];
    this->gg = new double[ngmw];

    for (int ig=0; ig<ngmw; ig++)
    {
        const int igg = (ig * NPROC_IN_POOL) + RANK_IN_POOL;
        if (igg < ngmw_g)
        {
            gdirect[ig] = gdirect_global[igg];
            gcar[ig] = gcar_global[igg];
            gg[ig] = gg_global[igg];
        }
    }

    // (2) second part
    delete[] ig1;
    delete[] ig2;
    delete[] ig3;
    this->ig1 = new int[ngmw];
    this->ig2 = new int[ngmw];
    this->ig3 = new int[ngmw];

//    cout << "\n dim of ig1,2,3 = " << ngmw;

    for (int i = 0; i < ngmw; i++)
    {
        this->ig1[i] = int(this->gdirect[i].x) + nx;
        this->ig2[i] = int(this->gdirect[i].y) + ny;
        this->ig3[i] = int(this->gdirect[i].z) + nz;
    }

    return;
}//end setup_mpi_GVectors
#else
void PW_Basis::get_GVectors(void)
{
    TITLE("PW_Basis","get_GVectors");
    timer::tick("PW_Basis","get_GVectors");

    //************************************************************
    // g  : Store the G vectors in 1d array (Cartian coordinate)
    // ig : Store the G vectors in 1d array (Direct coordinate)
    // gg : store the |G|^2 of the 1d array
    //************************************************************

    //----------------------------------------------------------
    // EXPLAIN : if not parallel case, we use pointer
    // g and g_global are pointers to the same array
    //----------------------------------------------------------
    this->gcar = this->gcar_global;
    this->gdirect = this->gdirect_global;
    this->gg = this->gg_global;

    // (2) calculate ig1, ig2, ig3
    assert(ngmw>0);
    delete[] ig1;
    delete[] ig2;
    delete[] ig3;
    this->ig1 = new int[ngmw];
    this->ig2 = new int[ngmw];
    this->ig3 = new int[ngmw];

    for (int i = 0; i < ngmw; i++)
    {
        this->ig1[i] = int(this->gdirect[i].x) + nx;
        this->ig2[i] = int(this->gdirect[i].y) + ny;
        this->ig3[i] = int(this->gdirect[i].z) + nz;
    }

    timer::tick("PW_Basis","get_GVectors");
    return;
}//end get_GVectors;
#endif


//  Calculate structure factor
void PW_Basis::setup_structure_factor(void)
{
    TITLE("PW_Basis","setup_structure_factor");
    timer::tick("PW_Basis","setup_struc_factor");
    complex<double> ci_tpi = NEG_IMAG_UNIT * TWO_PI;
    complex<double> x;

    this->strucFac = new complex<double>*[NTYPE];
    for (int it=0; it<NTYPE; it++)
    {
        this->strucFac[it] = new complex<double>[ngmw];
        ZEROS( strucFac[it], ngmw);
    }

    for (int it=0; it< NTYPE; it++)
    {
        for (int ig=0; ig<ngmw; ig++)
        {
            double sum_cos = 0.0;
            double sum_sin = 0.0;
            for (int ia=0; ia< NA[it]; ia++)
            {
                //----------------------------------------------------------
                // EXPLAIN : Don't use Dot function until we can optimize
                // it, use the following x*x + y*y + z*z instead!
                //----------------------------------------------------------
                // e^{-i G*tau}

                const double theta = TWO_PI * (
                                         gcar[ig].x * CARPOSX[it][ia] +
                                         gcar[ig].y * CARPOSY[it][ia] +
                                         gcar[ig].z * CARPOSZ[it][ia] );
                sum_cos += cos( theta );
                sum_sin += sin( theta );
            }
            this->strucFac[it][ig] = complex<double>( sum_cos, -sum_sin );

            double tmpx = strucFac[it][ig].real() ;
            double tmpy = strucFac[it][ig].imag() ;
        }
    }

    int nat = 0;
    for (int it=0; it<NTYPE; it++)
    {
        nat += NA[it];
    }

    int i,j,ng;
    this->eigts1.create(nat, 2*this->nx + 1);
    this->eigts2.create(nat, 2*this->ny + 1);
    this->eigts3.create(nat, 2*this->nz + 1);

    Vector3<double> gtau;
    Vector3<double> tau;
    int iat = 0;
    for (int it = 0; it < NTYPE; it++)
    {
        for (int ia = 0; ia < NA[it]; ia++)
        {
            tau.x = CARPOSX[it][ia];
            tau.y = CARPOSY[it][ia];
            tau.z = CARPOSZ[it][ia];

            gtau = this->G * tau;  //HLX: fixed on 10/13/2006

            for (int n1 = -nx; n1 <= nx; n1++)
            {
                double arg = n1 * gtau.x;
                this->eigts1(iat, n1 + nx) = exp( ci_tpi*arg  );
            }
            for (int n2 = -ny; n2 <= ny; n2++)
            {
                double arg = n2 * gtau.y;
                this->eigts2(iat, n2 + ny) = exp( ci_tpi*arg );
            }
            for (int n3 = -nz; n3 <= nz; n3++)
            {
                double arg = n3 * gtau.z;
                this->eigts3(iat, n3 + nz) = exp( ci_tpi*arg );
            }
            iat++;
        }
    }
    timer::tick("PW_Basis","setup_struc_factor");
    return;
}

complex<double>* PW_Basis::get_sk(const int ik, const int it, const int ia)const
{
    timer::tick("PW_Basis","get_sk");
    const double arg = (CARKX[ik] * CARPOSX[it][ia]
                        + CARKY[ik] * CARPOSY[it][ia]
                        + CARKZ[ik] * CARPOSZ[it][ia] ) * TWO_PI;
    const complex<double> kphase = complex <double> ( cos(arg),  -sin(arg) );
    complex<double> *sk = new complex<double>[ this->ngk[ik] ];
    const int iat = this->itia2iat[it][ia];

    for (int ig=0; ig< this->ngk[ik]; ig++)
    {
//		cout << "\n ig=" << ig << endl;
        const int iig = this->igk[ik][ig];
        assert( iig < ngmw );

        if (!(ig1[iig] >= 0 && ig1[iig] < 2*nx + 1) )
        {
            cout << "\n ig1[iig] = " << ig1[iig] << endl;
			QUIT();
        }

        if (!(ig2[iig] >= 0 && ig2[iig] < 2*ny + 1) )
        {
            cout << "\n ig2[iig] = " << ig2[iig] << endl;
			QUIT();
        }

        if (!(ig3[iig] >= 0 && ig3[iig] < 2*nz + 1) )
        {
            cout << "\n iig = " << iig << endl;
            cout << "\n ig3[iig] = " << ig3[iig] << endl;
			QUIT();
        }

        sk[ig] = kphase
                 * this->eigts1(iat, this->ig1[iig])
                 * this->eigts2(iat, this->ig2[iig])
                 * this->eigts3(iat, this->ig3[iig]);
    }
    timer::tick("PW_Basis","get_sk");
    return sk;
}

Vector3<double> PW_Basis::get_1qvec_cartesian(const int ik,const int ig)const
{
    Vector3<double> kvec = Vector3<double>(CARKX[ik], CARKY[ik], CARKZ[ik]);
    Vector3<double> qvec = kvec + this->gcar[ this->igk[ik][ig] ];
    return qvec;
}

void PW_Basis::table(void)
{
    TITLE("PW_Basis", "table");
    NBasis.init_table();

    const int total_lm = ( LMAXUSED + 1) * ( LMAXUSED + 1);

    delete[] ylm;
    this->ylm = new matrix[NKS];

    for (int ik=0; ik<NKS; ik++)
    {
        const int npw = ngk[ik];
        this->ylm[ik].create(total_lm, npw);

        // (1) get the vector k+G
        Vector3<double> *gk = new Vector3 <double> [npw];
        for (int ig=0; ig<npw; ig++)
        {
            gk[ig] = this->get_1qvec_cartesian(ik, ig);
        }

        // (2) generate ylm according to different k.
        // (or we can save the ylm for only once but all k points.
        Ylm_Real(total_lm, npw, gk, this->ylm[ik]);

        delete[] gk;
    }


	// init index of iw00
	assert(NWFCALL>0);
	this->iwindex = new Way2iw[NWFCALL];

	int iw=0;
	for(int it=0; it<NTYPE; it++)
	{
		for(int ia=0; ia<NA[it]; ia++)
		{
			for(int l=0; l<LMAXALL+1; l++)
			{
				for(int m=0; m<2*l+1; m++)
				{
					iwindex[iw].type = it;
					iwindex[iw].i = ia;
					iwindex[iw].L = l;
					iwindex[iw].m = m;
					++iw;
				}
			}
		}
	}	
	assert(iw==NWFCALL);

    return;
}


// init psi1d for each Level.
void PW_Basis::allocate_psi1d(const int &il)
{
    TITLE("PW_Basis", "allocate_psi1d");

    this->nwfc2 = mz.Level[il].nwfc2;

    this->Dk = 0.01;
    this->kmesh = static_cast<int>(sqrt(ECUT) / Dk) + 1 + 4;
    if (kmesh % 2 == 0) ++kmesh;

    this->psi1d.create(NTYPE, LMAXUSED+1, NMAXUSED, kmesh);

	delete[] save;
	delete[] posmu;
	delete[] posnu;

	const int dim = nwfc2 * nwfc2;
	save = new complex<double>[dim];
	posmu = new int[dim];
	posnu = new int[dim];
	ZEROS(save, dim);
	ZEROS(posmu, dim);
	ZEROS(posnu, dim);
    return;
}

void PW_Basis::calculate_psi1d(void)
{
    timer::tick("PW_Basis","calculate_psi1d");
    psi1d.zero_out();
    for (int it=0; it<NTYPE; it++)
    {
        for (int l=0; l<LMAXUSED+1; l++)
        {
            for (int n=0; n<NMAXUSED; n++)
            {
                for (int ie=0; ie<NE; ie++)
                {
                    for (int ikm=0; ikm<kmesh; ikm++)
                    {
                        this->psi1d(it,l,n,ikm) += input.Coef.C4(it, l, n ,ie)
                                                   * Numerical_Basis::bessel_basis.TableOne(l, ie, ikm);
                    }
                }
            }
        }
    }
    timer::tick("PW_Basis","calculate_psi1d");
    return;
}

void PW_Basis::allocate_psi3d(const int &level)
{

    // allocate the number of wave functions.
    // consider(it, ia, l, n, m)
    delete[] psi3d;
    psi3d = new ComplexMatrix[NKS];

    for (int ik=0; ik<NKS; ik++)
    {
        this->psi3d[ik].create(nwfc2, npwx);
    }
    return;
}


void PW_Basis::calculate_psi3d(const int &ilevel, const int &ik)
{
    //TITLE("PW_Basis", "calculate_psi3d");
    timer::tick("PW_Basis", "calculate_psi3d");
    // (1) get the vector k+G

    ofs_running << " calculate_psi3d" << endl;
    ofs_running << " ilevel = " << ilevel << " ik = " << ik << endl;

    const int npw = ngk[ik];

    Vector3<double> *gk = new Vector3 <double> [npw];
    for (int ig=0; ig<npw; ig++)
    {
        gk[ig] = this->get_1qvec_cartesian(ik, ig);
    }

    double *flq = new double[npw];
    for (int iw=0; iw<this->nwfc2; iw++)
    {
        // for each wave function in nwfc2(it, i, l, n)
        const int it = mz.Level[ilevel].wayd[iw].type;
        const int ia = mz.Level[ilevel].wayd[iw].i;
        const int l  = mz.Level[ilevel].wayd[iw].L;
        const int n  = mz.Level[ilevel].wayd[iw].N;
        const int m  = mz.Level[ilevel].wayd[iw].m;
        complex<double> lphase = pow(IMAG_UNIT, l);

        // get flq from psi1d for each k point
        for (int ig=0; ig<npw; ig++)
        {
            flq[ig] = this->Polynomial_Interpolation(it, l, n, gk[ig].norm() * this->tpiba );
        }

        // get the structure wave functions for each k point.
        complex<double> *sk = this->get_sk(ik, it, ia);

        const int lm = l*l+m;
        for (int ig=0; ig<npw; ig++)
        {
            psi3d[ik](iw,ig) = lphase * sk[ig] * ylm[ik](lm, ig) * flq[ig];
        }

        delete[] sk;
    }

    delete[] flq;
    delete[] gk;
    timer::tick("PW_Basis", "calculate_psi3d");
    return;
}

complex<double> PW_Basis::calculateS(const int &iw, const int &iw2, const int &ik)
{
    complex<double> overlap = 0.0;
    for (int ig=0; ig<ngk[ik]; ig++)
    {
        overlap += conj(this->psi3d[ik](iw, ig)) * this->psi3d[ik](iw2, ig);
    }
// move this outside to get_spillage
#ifdef __MPI
//    Parallel_Reduce::reduce_complex_double_pool(overlap);
#endif
    return overlap;
}


void PW_Basis::update_psi1d(const int &il, const int &ic, const int &ie,
                            const double &c4_now, const double &c4_old)
{
    timer::tick("PW_Basis","update_psi1d");

    // a small trick
    // ALGORITHM : when we change psi1d, we don't know
    // this psi1d will be accepted or rejected.
    // However, it has been changed.
    // So, in the next time, we recalculate them
    // using C4.
    static int it_last = -1;
    static int l_last = -1;
    static int n_last = -1;

    const int it_now = mz.Level[il].wayc[ic].type;
    const int l_now = mz.Level[il].wayc[ic].L;
    const int n_now = mz.Level[il].wayc[ic].N;

    for (int it=0; it<NTYPE; it++)
    {
        for (int l=0; l<input.Coef.C4.getBound2(); l++)
        {
            for (int n=0; n<input.Coef.C4.getBound3(); n++)
            {
                if (it==it_last || l==l_last || n==n_last)
                {
                    for (int ikm=0; ikm<kmesh; ikm++)
                    {
						// old version
                        double sum = 0.0;
                        for (int ie2=0; ie2<NE; ie2++)
                        {
                            sum += input.Coef.C4(it, l, n ,ie2)
                                   * Numerical_Basis::bessel_basis.TableOne(l, ie2, ikm);
                        }
                        this->psi1d(it,l,n,ikm) = sum;
                    }
                }
				else if(it==it_now || l==l_now || n==n_now)
				{
					for (int ikm=0; ikm<kmesh; ikm++)
					{
						// a small trick
						this->psi1d(it,l,n,ikm) += (c4_now-c4_old) * Numerical_Basis::bessel_basis.TableOne(l, ie, ikm);
					}
				}
            }
        }
    }

    it_last = it_now;
    l_last = l_now;
    n_last = n_now;

    timer::tick("PW_Basis","update_psi1d");
    return;
}

void PW_Basis::update_psi3d( const int &il, const int &ic, const int &ik)
{
//    TITLE("PW_Basis","update_psi1d3d");
    timer::tick("PW_Basis","update_psi3d");

    const int npw = ngk[ik];
    static int it_old;
    static int ia_old;
    static int l_old;
    static int n_old;
    static bool allocate=false;
    static int *ic_last;
    if (!allocate)
    {
        ic_last = new int[NKS];
        for (int ik=0; ik<NKS; ik++) ic_last[ik] = -1;
        allocate=true;
    }

    it_old = -1;
    ia_old = -1;
    l_old = -1;
    n_old = -1;


    // (1) get the vector k+G
    Vector3<double> *gk = new Vector3 <double> [npw];
    for (int ig=0; ig<npw; ig++)
    {
        gk[ig] = this->get_1qvec_cartesian(ik, ig);
    }
    double *flq = new double[npw];
    complex<double> *sk = new complex<double>[1];
    complex<double> *samepart = new complex<double>[npw];
    for (int iw=0; iw<this->nwfc2; iw++)
    {
        const int ic1 = mz.Level[il].wayd[iw].ic;
        if (ic1 == ic || ic1 == ic_last[ik])
        {
            // for each wave function in nwfc2(it, i, l, n)
            const int it = mz.Level[il].wayd[iw].type;
            const int ia = mz.Level[il].wayd[iw].i;
            const int l  = mz.Level[il].wayd[iw].L;
            const int n  = mz.Level[il].wayd[iw].N;
            const int m  = mz.Level[il].wayd[iw].m;
            complex<double> lphase = pow(IMAG_UNIT, l);

            // a small trick: get the structure wave functions for each k point.
            if (it==it_old && ia==ia_old)
            {
                // do nothing
            }
            else
            {
                delete[] sk;
                sk = this->get_sk(ik, it, ia);
            }

            // small trick (2) : get flq from psi1d for each k point
            if (it==it_old && l==l_old && n==n_old)
            {
                // do nothing
            }
            else
            {
                for (int ig=0; ig<npw; ig++)
                {
                    flq[ig] = this->Polynomial_Interpolation(it, l, n, gk[ig].norm() * this->tpiba );
                }
            }

            // small trick (3)
            if (it==it_old && ia==ia_old && l==l_old && n==n_old)
            {
                //samepart is the same as the previous one.
            }
            else
            {
                for (int ig=0; ig<npw; ig++)
                {
                    samepart[ig] = lphase * sk[ig] * flq[ig];
                }
//				cout << "\n update samepart=" << samepart[35];
            }

            // final calculation
            const int lm = l*l+m;
            for (int ig=0; ig<npw; ig++)
            {
                //  psi3d[ik](iw,ig) = ylm[ik](lm, ig) * samepart[ig];
                psi3d[ik](iw,ig) = ylm[ik](lm, ig) * lphase * sk[ig] * flq[ig];

//				if(ig==35)
//				{
//					cout << "\n samepart=" << samepart[ig] << " haha=" << lphase * sk[ig] * flq[ig];
//					int ok; cin >> ok;
//				}
            }


            it_old = it;
            ia_old = ia;
            l_old = l;
            n_old = n;
        }
    }

    ic_last[ik] = ic;
    delete[] samepart;
    delete[] gk;
    delete[] flq;
	delete[] sk;
    // if not accepted, need to go back to origin psi3d.
    timer::tick("PW_Basis","update_psi3d");
    return;
}

double PW_Basis::Polynomial_Interpolation
(const int &it, const int &l, const int &n, const double &gnorm)const
{
    const double position =  gnorm / this->Dk;
    const int iq = static_cast<int>(position);
    assert(iq < kmesh-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y=
        this->psi1d(it, l, n, iq) * x1 * x2 * x3 / 6.0 +
        this->psi1d(it, l, n, iq) * x0 * x2 * x3 / 2.0 -
        this->psi1d(it, l, n, iq) * x1 * x0 * x3 / 2.0 +
        this->psi1d(it, l, n, iq) * x1 * x2 * x0 / 6.0 ;
    return y;
}


void PW_Basis::calculate_Jlq(const int &ik, const int &iw, const int &ie)
{
    //TITLE("PW_Basis", "calculate_Jlq");
	timer::tick("PW_Basis","calculate_Jlq");

    const int npw = this->ngk[ik];
	const int it= iwindex[iw].type; 
	const int ia= iwindex[iw].i;
	const int l = iwindex[iw].L;
	const int m = iwindex[iw].m;

    // generate 1d Jlq.
    Vector3<double> *gk = new Vector3 <double> [npw];
    double* flq=new double[npw];
    for (int ig=0; ig<npw; ig++)
    {
        gk[ig] = this->get_1qvec_cartesian(ik, ig);
    }
    for (int ig=0; ig<npw; ig++)
    {
        flq[ig] = Numerical_Basis::bessel_basis.Polynomial_Interpolation2(l, ie, gk[ig].norm() * this->tpiba );
    }

    // generate 3d Jlq.
    delete[] Jlq;
    Jlq = new complex<double>[npw];
    complex<double> *sk = this->get_sk(ik, it, ia);

    const int lm = l*l+m;
    complex<double> lphase = pow(IMAG_UNIT, l);
    for (int ig=0; ig<npw; ig++)
    {
        Jlq[ig] = lphase * sk[ig] * ylm[ik](lm, ig) * flq[ig];
    }

    delete[] sk;
	delete[] gk;
	delete[] flq;

	timer::tick("PW_Basis","calculate_Jlq");
    return;
}



complex<double> PW_Basis::calculate_Jlq_Phi(const int &ik, const int &mu)
{
	const int npw = this->ngk[ik];
	complex<double> overlap = complex<double>(0,0);
	for(int ig=0; ig<npw; ig++)
	{
		overlap += conj(Jlq[ig]) * psi3d[ik](mu,ig);

	}
#ifdef __MPI
    Parallel_Reduce::reduce_complex_double_pool(overlap);
#endif
	return overlap;
}

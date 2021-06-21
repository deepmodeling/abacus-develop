//==========================================================
// AUTHOR : Peize Lin
// DATE : 2015-03-10
//==========================================================

#include "exx_lip.h"
#include "src_global/global_function.h"
#include "src_global/vector3.h"
#include "src_pw/global.h"
#include "src_pw/klist.h"
#include "src_pw/wavefunc.h"
#include "src_pw/wavefunc_in_pw.h"
#include "src_global/lapack_connector.h"
#include <limits>

#include "src_external/src_test/test_function.h"

Exx_Lip::Exx_Lip( const Exx_Global::Exx_Info &info_global )
	:init_finish(false),
	 info(info_global),
	 exx_matrix(NULL),
	 exx_energy(0){}

Exx_Lip::Exx_Info::Exx_Info( const Exx_Global::Exx_Info &info_global )
	:hybrid_type(info_global.hybrid_type),
	 hse_omega(info_global.hse_omega){} 

void Exx_Lip::cal_exx()
{
	TITLE("Exx_Lip","cal_exx");
	auto my_time = [](timeval &t_begin) -> double
	{
		const double time_during = cal_time(t_begin);
		gettimeofday(&t_begin, NULL);
		return time_during;
	};
	auto cout_t = [](const string &name, const double t)
	{
		cout<<name<<"\t"<<t<<endl;
	};
	
timeval t;
gettimeofday(&t, NULL);
double t_phi_cal=0, t_qkg2_exp=0, t_b_cal=0, t_sum3_cal=0, t_b_sum=0, t_sum_all=0;
	wf_wg_cal();
cout_t("wf_wg_cal",my_time(t));
	psi_cal();
cout_t("psi_cal",my_time(t));
	for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
	{
		phi_cal(k_pack, ik);
t_phi_cal += my_time(t);

		judge_singularity(ik);
		for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
				sum1[iw_l*NLOCAL+iw_r] = (0.0,0.0);
		if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
		{			
			sum2_factor = 0.0;
			if(gzero_rank_in_pool==RANK_IN_POOL)
				for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
					for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
						sum3[iw_l][iw_r] = (0.0,0.0);
		}
		
		for( int iq_tmp=iq_vecik; iq_tmp<iq_vecik+q_pack->kv_ptr->nks/NSPIN; ++iq_tmp)					// !!! k_point parallel incompleted. need to loop iq in other pool
		{
			int iq = (ik<(k_pack->kv_ptr->nks/NSPIN)) ? (iq_tmp%(q_pack->kv_ptr->nks/NSPIN)) : (iq_tmp%(q_pack->kv_ptr->nks/NSPIN)+(q_pack->kv_ptr->nks/NSPIN));
			qkg2_exp(ik, iq);
t_qkg2_exp += my_time(t);
			for( int ib=0; ib<NBANDS; ++ib)
			{
				b_cal(ik, iq, ib);
t_b_cal += my_time(t);
				if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
					if(iq==iq_vecik)
						sum3_cal(iq,ib);
t_sum3_cal += my_time(t);
				b_sum(iq, ib);
t_b_sum += my_time(t); 
			}
		}
		sum_all(ik);
t_sum_all += my_time(t);
	}
	exx_energy_cal();
cout_t("exx_energy_cal",my_time(t));
cout_t("phi_cal",t_phi_cal);
cout_t("qkg2_exp",t_qkg2_exp);
cout_t("b_cal",t_b_cal);
cout_t("sum3_cal",t_sum3_cal);
cout_t("b_sum",t_b_sum);
cout_t("sum_all",t_sum_all);

	auto print_Hexxk = [&]()
	{
		static int istep=1;
		for(int ik=0; ik!=kv.nks; ++ik)
		{
			ofstream ofs("Hexxk_"+TO_STRING(istep++)+"_"+TO_STRING(ik)+"_"+TO_STRING(MY_RANK));
			for(int i=0; i!=NLOCAL; ++i)
			{
				for(int j=0; j!=NLOCAL; ++j)
					ofs<<exx_matrix[ik][i][j]<<"\t";
				ofs<<endl;
			}
		};
	};
}

/*
void Exx_Lip::cal_exx()
{
	TITLE("Exx_Lip","cal_exx");
	wf_wg_cal();
	psi_cal();
	for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
	{
		phi_cal(k_pack, ik);

		judge_singularity(ik);
		for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
				sum1[iw_l*NLOCAL+iw_r] = (0.0,0.0);
		if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
		{			
			sum2_factor = 0.0;
			if(gzero_rank_in_pool==RANK_IN_POOL)
				for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
					for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
						sum3[iw_l][iw_r] = (0.0,0.0);
		}
		
		for( int iq_tmp=iq_vecik; iq_tmp<iq_vecik+q_pack->kv_ptr->nks/NSPIN; ++iq_tmp)					// !!! k_point parallel incompleted. need to loop iq in other pool
		{
			int iq = (ik<(k_pack->kv_ptr->nks/NSPIN)) ? (iq_tmp%(q_pack->kv_ptr->nks/NSPIN)) : (iq_tmp%(q_pack->kv_ptr->nks/NSPIN)+(q_pack->kv_ptr->nks/NSPIN));
			qkg2_exp(ik, iq);
			for( int ib=0; ib<NBANDS; ++ib)
			{
				b_cal(ik, iq, ib);
				if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
					if(iq==iq_vecik)
						sum3_cal(iq,ib);
					b_sum(iq, ib);
			}
		}
		sum_all(ik);
	}
	exx_energy_cal();
}
*/

void Exx_Lip::init(K_Vectors *kv_ptr_in, wavefunc *wf_ptr_in, PW_Basis *pw_ptr_in, Use_FFT *UFFT_ptr_in, UnitCell_pseudo *ucell_ptr_in)
{
	TITLE("Exx_Lip","init");
	try
	{
		k_pack = new k_package;
		k_pack->kv_ptr = kv_ptr_in;
		k_pack->wf_ptr = wf_ptr_in;
		pw_ptr = pw_ptr_in;
		UFFT_ptr = UFFT_ptr_in;
		ucell_ptr = ucell_ptr_in;

		int gzero_judge(-1);
		if (pw_ptr->gcar[0]==Vector3<double>(0.0,0.0,0.0))
		{
			gzero_judge = RANK_IN_POOL;
		}
		MPI_Allreduce(&gzero_judge, &gzero_rank_in_pool, 1, MPI_INT, MPI_MAX, POOL_WORLD);

		k_pack->wf_wg.create(k_pack->kv_ptr->nks,NBANDS);

		k_pack->hvec_array = new ComplexMatrix [k_pack->kv_ptr->nks];
		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
		{
			k_pack->hvec_array[ik].create(NLOCAL,NBANDS);
		}

		if (pot.start_pot=="atomic")
		{
			q_pack = k_pack;
		}
		else if(pot.start_pot=="file")
		{
			read_q_pack();
		}

		phi = new complex<double>*[NLOCAL];
		for( int iw=0; iw<NLOCAL; ++iw)
		{
			phi[iw] = new complex<double>[pw_ptr->nrxx];
		}

		psi = new complex<double>**[q_pack->kv_ptr->nks];
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			psi[iq] = new complex<double> *[NBANDS];
			for( int ib=0; ib<NBANDS; ++ib)
			{
				psi[iq][ib] = new complex<double>[pw_ptr->nrxx];
			}
		}

		recip_qkg2 = new double [pw_ptr->ngmc];

		b = new complex<double> [NLOCAL*pw_ptr->ngmc];

		sum1 = new complex<double> [NLOCAL*NLOCAL];

		if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
			if(gzero_rank_in_pool==RANK_IN_POOL)
			{
				b0 = new complex<double> [NLOCAL];
				sum3 = new complex<double> *[NLOCAL];
				for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
				{
					sum3[iw_l] = new complex<double> [NLOCAL];
				}
			}
			else
			{
				b0 = NULL;
				sum3 = NULL;
			}
		else
		{
			b0 = NULL;
			sum3 = NULL;
		}

		exx_matrix = new complex<double> **[k_pack->kv_ptr->nks];
		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
		{
			exx_matrix[ik] = new complex<double>*[NLOCAL];
			for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			{
				exx_matrix[ik][iw_l] = new complex<double>[NLOCAL];
			}
		}
	}
	catch(const std::bad_alloc &ex)
	{
		WARNING_QUIT("exact_exchange","Memory");
	}

	init_finish = true;
}

Exx_Lip::~Exx_Lip()
{
	TITLE("Exx_Lip","~Exx_Lip");
	if( init_finish)
	{
		for( int iw=0; iw<NLOCAL; ++iw)
		{
			delete[] phi[iw];	phi[iw]=NULL;
		}
		delete[] phi;		phi=NULL;

		for( int iq=0;iq<q_pack->kv_ptr->nks;++iq)
		{
			for( int ib=0;ib<NBANDS;++ib)
			{
				delete[] psi[iq][ib];	psi[iq][ib]=NULL;
			}
			delete[] psi[iq];	psi[iq]=NULL;
		}
		delete[] psi;		psi=NULL;

		delete[] recip_qkg2;	recip_qkg2=NULL;

		delete[] b;		b=NULL;

		delete[] sum1;		sum1=NULL;

		if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
			if(gzero_rank_in_pool==RANK_IN_POOL)
			{
				delete[] b0;	b0=NULL;
				for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
				{
					delete[] sum3[iw_l];	sum3[iw_l]=NULL;
				}
				delete[] sum3;		sum3=NULL;
			}

		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
		{
			for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			{
				delete[] exx_matrix[ik][iw_l];	exx_matrix[ik][iw_l]=NULL;
			}
			delete[] exx_matrix[ik];	exx_matrix[ik]=NULL;
		}
		delete[] exx_matrix;	exx_matrix=NULL;

		delete[] k_pack->hvec_array;	k_pack->hvec_array=NULL;
		delete k_pack;

		if (pot.start_pot=="atomic")
		{
			q_pack = NULL;
		}
		else if(pot.start_pot=="file")
		{
			delete q_pack->kv_ptr;	q_pack->kv_ptr=NULL;
			delete q_pack->wf_ptr;	q_pack->wf_ptr=NULL;
			delete[] q_pack->hvec_array;	q_pack->hvec_array=NULL;
			delete q_pack;	q_pack=NULL;
		}
	}
}

void Exx_Lip::wf_wg_cal()
{
	TITLE("Exx_Lip","wf_wg_cal");
	if(NSPIN==1)
		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
			for( int ib=0; ib<NBANDS; ++ib)
				k_pack->wf_wg(ik,ib) = k_pack->wf_ptr->wg(ik,ib)/2;
	else if(NSPIN==2)
		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
			for( int ib=0; ib<NBANDS; ++ib)
				k_pack->wf_wg(ik,ib) = k_pack->wf_ptr->wg(ik,ib);
}

void Exx_Lip::phi_cal(k_package *kq_pack, int ikq)
{
	for( int iw=0; iw< NLOCAL; ++iw)
	{
		ZEROS( UFFT_ptr->porter, pw_ptr->nrxx );
		for( int ig=0; ig<kq_pack->kv_ptr->ngk[ikq]; ++ig)
			UFFT_ptr->porter[ pw_ptr->ig2fftw[kq_pack->wf_ptr->igk(ikq,ig)] ] = kq_pack->wf_ptr->wanf2[ikq](iw,ig);
		pw_ptr->FFT_wfc.FFT3D(UFFT_ptr->porter,1);
		int ir(0);
		for( int ix=0; ix<pw_ptr->ncx; ++ix)
		{
			const double phase_x = kq_pack->kv_ptr->kvec_d[ikq].x * ix / pw_ptr->ncx;
			for( int iy=0; iy<pw_ptr->ncy; ++iy)
			{
				const double phase_xy = phase_x + kq_pack->kv_ptr->kvec_d[ikq].y * iy / pw_ptr->ncy;
				for( int iz=pw_ptr->nczp_start; iz<pw_ptr->nczp_start+pw_ptr->nczp; ++iz)
				{
					const double phase_xyz = phase_xy + kq_pack->kv_ptr->kvec_d[ikq].z * iz / pw_ptr->ncz;
					const complex<double> exp_tmp = exp(phase_xyz*TWO_PI*IMAG_UNIT);
					phi[iw][ir] = UFFT_ptr->porter[ir]*exp_tmp;
					++ir;
				}
			}
		}
	}
}

void Exx_Lip::psi_cal()
{
	TITLE("Exx_Lip","psi_cal");
	if (pot.start_pot=="atomic")
	{
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			for( int ib = 0; ib < NBANDS; ++ib)
			{
				ZEROS( UFFT_ptr->porter, pw_ptr->nrxx );
				for( int ig = 0; ig < q_pack->kv_ptr->ngk[iq] ; ++ig)
				{
					UFFT_ptr->porter[ pw_ptr->ig2fftw[q_pack->wf_ptr->igk(iq,ig)] ] = q_pack->wf_ptr->evc[iq](ib,ig);
				}
				pw_ptr->FFT_wfc.FFT3D(UFFT_ptr->porter,1);
				int ir(0);
				for( int ix=0; ix<pw_ptr->ncx; ++ix)
				{
					const double phase_x = q_pack->kv_ptr->kvec_d[iq].x * ix / pw_ptr->ncx;
					for( int iy=0; iy<pw_ptr->ncy; ++iy)
					{
						const double phase_xy = phase_x + q_pack->kv_ptr->kvec_d[iq].y * iy / pw_ptr->ncy;
						for( int iz=pw_ptr->nczp_start; iz<pw_ptr->nczp_start+pw_ptr->nczp; ++iz)
						{
							const double phase_xyz = phase_xy + q_pack->kv_ptr->kvec_d[iq].z * iz / pw_ptr->ncz;
							const complex<double> exp_tmp = exp(phase_xyz*TWO_PI*IMAG_UNIT);
							psi[iq][ib][ir] = UFFT_ptr->porter[ir]*exp_tmp;
							++ir;
						}
					}
				}
			}
		}
	}
	else if(pot.start_pot=="file")
	{
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			phi_cal( q_pack, iq);
			for( int ib=0; ib<NBANDS; ++ib)
			{
				ZEROS(psi[iq][ib],pw_ptr->nrxx);
				for( int iw=0; iw<NLOCAL; ++iw)
				{
					for( int ir=0; ir<pw_ptr->nrxx; ++ir)
					{
						psi[iq][ib][ir] += q_pack->hvec_array[iq](iw,ib) * phi[iw][ir];
					}
				}
			}
		}
	}
///////////////////////////////////////////////////////////////////////////
//			!!! k_point parallel incompleted. need to loop iq in other pool
///////////////////////////////////////////////////////////////////////////
}

void Exx_Lip::judge_singularity( int ik)
{
	if (pot.start_pot=="atomic")
	{
		iq_vecik = ik;
	}
	else if(pot.start_pot=="file")
	{
		double min_q_minus_k(numeric_limits<double>::max());
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			const double q_minus_k ( (q_pack->kv_ptr->kvec_c[iq] - k_pack->kv_ptr->kvec_c[ik]).norm2() );
			if(q_minus_k < min_q_minus_k)
			{
				min_q_minus_k = q_minus_k;
				iq_vecik = iq;
			}
		}
	}
}


void Exx_Lip::qkg2_exp(int ik, int iq)
{
	for( int ig=0; ig<pw_ptr->ngmc; ++ig)
	{
		const double qkg2 = ( (q_pack->kv_ptr->kvec_c[iq] - k_pack->kv_ptr->kvec_c[ik] + pw_ptr->gcar[ig]) *(TWO_PI/ucell_ptr->lat0)).norm2();
		if( (Exx_Global::Hybrid_Type::PBE0==info.hybrid_type) || (Exx_Global::Hybrid_Type::HF==info.hybrid_type) )
		{
			if( abs(qkg2)<1e-10 )
				recip_qkg2[ig] = 0.0;												// 0 to ignore bb/qkg2 when qkg2==0
			else
				recip_qkg2[ig] = 1.0/qkg2;
			sum2_factor += recip_qkg2[ig] * exp(-info.lambda*qkg2) ;
			recip_qkg2[ig] = sqrt(recip_qkg2[ig]);
		}
		else if(Exx_Global::Hybrid_Type::HSE==info.hybrid_type)
		{
			if( abs(qkg2)<1e-10 )
				recip_qkg2[ig] = 1.0/(2*info.hse_omega);
			else
				recip_qkg2[ig] = sqrt( (1-exp(-qkg2/(4*info.hse_omega*info.hse_omega))) /qkg2);
		}
	}
}


void Exx_Lip::b_cal( int ik, int iq, int ib)
{
	const Vector3<double> q_minus_k = q_pack->kv_ptr->kvec_d[iq] - k_pack->kv_ptr->kvec_d[ik];
	vector<complex<double> > mul_tmp(pw_ptr->nrxx);
	for( size_t ir=0,ix=0; ix<pw_ptr->ncx; ++ix)
	{
		const double phase_x = q_minus_k.x*ix/pw_ptr->ncx;
		for( size_t iy=0; iy<pw_ptr->ncy; ++iy)
		{
			const double phase_xy = phase_x + q_minus_k.y*iy/pw_ptr->ncy;
			for( size_t iz=pw_ptr->nczp_start; iz<pw_ptr->nczp_start+pw_ptr->nczp; ++iz)
			{
				const double phase_xyz = phase_xy + q_minus_k.z*iz/pw_ptr->ncz;
				mul_tmp[ir] = exp(-phase_xyz*TWO_PI*IMAG_UNIT);
				mul_tmp[ir] *= psi[iq][ib][ir];
				++ir;
			}
		}
	}

	complex<double> * const porter = UFFT_ptr->porter;
	const int * const ig2fftc = pw_ptr->ig2fftc;
	for(size_t iw=0; iw< NLOCAL; ++iw)
	{
		const complex<double> * const phi_w = phi[iw];
		for( size_t ir=0; ir<pw_ptr->nrxx; ++ir)
		{
			porter[ir] = conj(phi_w[ir]) * mul_tmp[ir] ;
//			porter[ir] = phi_w[ir] * psi_q_b[ir] *exp_tmp[ir] ;
		}
		pw_ptr->FFT_chg.FFT3D( porter, -1);
		if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
			if((iq==iq_vecik) && (gzero_rank_in_pool==RANK_IN_POOL))							/// need to check while use k_point parallel
				b0[iw] = porter[ pw_ptr->ig2fftc[0] ];
		complex<double> * const b_w = b+iw*pw_ptr->ngmc;
		for( size_t ig=0; ig<pw_ptr->ngmc; ++ig)
			b_w[ig] = porter[ ig2fftc[ig] ] * recip_qkg2[ig];
	}
}


void  Exx_Lip::sum3_cal(int iq, int ib)
{
	if( gzero_rank_in_pool == RANK_IN_POOL )
		for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
				sum3[iw_l][iw_r] += b0[iw_l] * conj(b0[iw_r]) * q_pack->wf_wg(iq,ib);
}


void Exx_Lip::b_sum( int iq, int ib)			// Peize Lin change 2019-04-14
{
	// sum1[iw_l,iw_r] += \sum_{ig} b[iw_l,ig] * conj(b[iw_r,ig]) * q_pack->wf_wg(iq,ib)
	LapackConnector::zherk(
		'U','N',
		NLOCAL, pw_ptr->ngmc,
		q_pack->wf_wg(iq,ib), b, pw_ptr->ngmc,
		1.0, sum1, NLOCAL);
//	cblas_zherk( CblasRowMajor, CblasUpper, CblasNoTrans,
//				NLOCAL, pw_ptr->ngmc,
//				q_pack->wf_wg(iq,ib), static_cast<void*>(b), pw_ptr->ngmc,
//				1.0, static_cast<void*>(sum1), NLOCAL);
}

void Exx_Lip::sum_all(int ik)
{
	double sum2_factor_g(0.0);
	if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
		MPI_Reduce( &sum2_factor, &sum2_factor_g, 1, MPI_DOUBLE, MPI_SUM, gzero_rank_in_pool, POOL_WORLD);

	for( size_t iw_l=1; iw_l<NLOCAL; ++iw_l)
		for( size_t iw_r=0; iw_r<iw_l; ++iw_r)
			sum1[iw_l*NLOCAL+iw_r] = conj(sum1[iw_r*NLOCAL+iw_l]);		// Peize Lin add conj 2019-04-14

	for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
	{
		for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
		{
			exx_matrix[ik][iw_l][iw_r] = 2.0* (-4*PI/ucell_ptr->omega *sum1[iw_l*NLOCAL+iw_r]);
			if( Exx_Global::Hybrid_Type::HF==info.hybrid_type || Exx_Global::Hybrid_Type::PBE0==info.hybrid_type )
				if(gzero_rank_in_pool==RANK_IN_POOL)
				{
					exx_matrix[ik][iw_l][iw_r] += 2.0* (4*PI/ucell_ptr->omega *sum3[iw_l][iw_r] *sum2_factor_g );
					exx_matrix[ik][iw_l][iw_r] += 2.0* (-1/sqrt(info.lambda*PI)*(q_pack->kv_ptr->nks/NSPIN) * sum3[iw_l][iw_r]);
				}
		}
	}
}

void Exx_Lip::exx_energy_cal()
{
	TITLE("Exx_Lip","exx_energy_cal");
	
	double exx_energy_tmp = 0.0;

	for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
	{
		for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
		{
			for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
			{
				for( int ib=0; ib<NBANDS; ++ib)
				{
					exx_energy_tmp += (exx_matrix[ik][iw_l][iw_r] *conj(k_pack->hvec_array[ik](iw_l,ib)) *k_pack->hvec_array[ik](iw_r,ib) ).real() *k_pack->wf_wg(ik,ib);
				}
			}
		}
	}
	MPI_Allreduce( &exx_energy_tmp, &exx_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);				// !!! k_point parallel incompleted. different pools have different kv.nks => deadlock
	exx_energy *= (NSPIN==1) ? 2 : 1;
	exx_energy /= 2;										// ETOT = E_band - 1/2 E_exx

	#if TEST_EXX==1
	{
		ofstream ofs("exx_matrix.dat",ofstream::app);
		static int istep=0;
		ofs<<"istep:\t"<<istep++<<endl;
		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
		{
			ofs<<"ik:\t"<<ik<<endl;
			for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			{
				for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
				{
					ofs<<exx_matrix[ik][iw_l][iw_r]<<"\t";
				}
				ofs<<endl;
			}
			ofs<<endl;
		}
		ofs.close();
	}
	{
		ofstream ofs("DM.dat",ofstream::app);
		static int istep=0;
		ofs<<"istep:\t"<<istep++<<endl;
		for( int ik=0; ik<k_pack->kv_ptr->nks; ++ik)
		{
			ofs<<"ik:\t"<<ik<<endl;
			for( int iw_l=0; iw_l<NLOCAL; ++iw_l)
			{
				for( int iw_r=0; iw_r<NLOCAL; ++iw_r)
				{
					complex<double> DM = {0,0};
					for( int ib=0; ib<NBANDS; ++ib )
						DM += conj(k_pack->hvec_array[ik](iw_l,ib)) *k_pack->hvec_array[ik](iw_r,ib) *k_pack->wf_wg(ik,ib);
					ofs<<DM<<"\t";
				}
				ofs<<endl;
			}
			ofs<<endl;
		}
		ofs.close();
	}
	#elif TEST_EXX==-1
		#error
	#endif

	return;
}

void Exx_Lip::write_q_pack() const
{
    if (CHR.out_charge==0)
		return;

	if(!RANK_IN_POOL)
	{
		const string exx_q_pack = "exx_q_pack/";

		const string command_mkdir = "test -d " + global_out_dir + exx_q_pack + " || mkdir " + global_out_dir + exx_q_pack;
		system( command_mkdir.c_str() );	// Need to check

		const string command_kpoint = "test -f " + global_out_dir + exx_q_pack + global_kpoint_card + " || cp " + global_kpoint_card + " " + global_out_dir + exx_q_pack + global_kpoint_card;
		system( command_kpoint.c_str() );	// Need to check

		stringstream ss_wf_wg;
		ss_wf_wg << global_out_dir << exx_q_pack << "wf_wg_" << MY_POOL;
		ofstream ofs_wf_wg(ss_wf_wg.str().c_str());
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			for( int ib=0; ib<NBANDS; ++ib)
			{
				ofs_wf_wg<<q_pack->wf_wg(iq,ib)<<"\t";
			}
			ofs_wf_wg<<endl;
		}
		ofs_wf_wg.close();

		stringstream ss_hvec;
		ss_hvec	<< global_out_dir << exx_q_pack << "hvec_" << MY_POOL;
		ofstream ofs_hvec(ss_hvec.str().c_str());
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			for( int iw=0; iw<NLOCAL; ++iw)
			{
				for( int ib=0; ib<NBANDS; ++ib)
				{
					ofs_hvec<<q_pack->hvec_array[iq](iw,ib).real()<<" "<<q_pack->hvec_array[iq](iw,ib).imag()<<" ";
				}
				ofs_hvec<<endl;
			}
		}
		ofs_hvec.close();
	}
	return;
}

void Exx_Lip::read_q_pack()
{
	const string exx_q_pack = "exx_q_pack/";

	q_pack = new k_package();

	q_pack->kv_ptr = new K_Vectors();
	const string exx_kpoint_card = global_out_dir + exx_q_pack + global_kpoint_card;
	q_pack->kv_ptr->set( symm, exx_kpoint_card, NSPIN, ucell_ptr->G, ucell_ptr->latvec );
//	q_pack->kv_ptr->set( symm, exx_kpoint_card, NSPIN, ucell_ptr->G, ucell_ptr->latvec, &Pkpoints );


	q_pack->wf_ptr = new wavefunc();
	q_pack->wf_ptr->allocate(q_pack->kv_ptr->nks); // mohan update 2021-02-25
//	q_pack->wf_ptr->init(q_pack->kv_ptr->nks,q_pack->kv_ptr,ucell_ptr,pw_ptr,&ppcell,&ORB,&hm,&Pkpoints);
	q_pack->wf_ptr->table_local.create(ucell.ntype, ucell.nmax_total, NQX);
//	q_pack->wf_ptr->table_local.create(q_pack->wf_ptr->ucell_ptr->ntype, q_pack->wf_ptr->ucell_ptr->nmax_total, NQX);
	Wavefunc_in_pw::make_table_q(ORB.orbital_file, q_pack->wf_ptr->table_local);
//	Wavefunc_in_pw::make_table_q(q_pack->wf_ptr->ORB_ptr->orbital_file, q_pack->wf_ptr->table_local, q_pack->wf_ptr);
	for(int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		Wavefunc_in_pw::produce_local_basis_in_pw(iq, q_pack->wf_ptr->wanf2[iq], q_pack->wf_ptr->table_local);
//		Wavefunc_in_pw::produce_local_basis_in_pw(iq, q_pack->wf_ptr->wanf2[iq], q_pack->wf_ptr->table_local, q_pack->wf_ptr);
	}

	q_pack->wf_wg.create(q_pack->kv_ptr->nks,NBANDS);
	if(!RANK_IN_POOL)
	{
		stringstream ss_wf_wg;
		ss_wf_wg << global_out_dir << exx_q_pack << "wf_wg_" << MY_POOL;
		ifstream ifs_wf_wg(ss_wf_wg.str().c_str());
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			for( int ib=0; ib<NBANDS; ++ib)
			{
				ifs_wf_wg>>q_pack->wf_wg(iq,ib);
			}
		}
		ifs_wf_wg.close();
	}
	MPI_Bcast( q_pack->wf_wg.c, q_pack->kv_ptr->nks*NBANDS, MPI_DOUBLE, 0, POOL_WORLD);

	q_pack->hvec_array = new ComplexMatrix [q_pack->kv_ptr->nks];
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		q_pack->hvec_array[iq].create(NLOCAL,NBANDS);
	}
	if(!RANK_IN_POOL)
	{
		stringstream ss_hvec;
		ss_hvec	<< global_out_dir << exx_q_pack << "hvec_" << MY_POOL;
		ifstream ifs_hvec(ss_hvec.str().c_str());
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			for( int iw=0; iw<NLOCAL; ++iw)
			{
				for( int ib=0; ib<NBANDS; ++ib)
				{
					double a,b;
					ifs_hvec>>a>>b;
					q_pack->hvec_array[iq](iw,ib) = {a,b};
				}
			}
		}
		ifs_hvec.close();
	}
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		MPI_Bcast( q_pack->hvec_array[iq].c, NLOCAL*NBANDS, mpicomplex, 0, POOL_WORLD);
	}

	return;
}

/*
void Exx_Lip::write_q_pack() const
{

	if( !RANK_IN_POOL )
	{
       	stringstream ssc;
        ssc << global_out_dir << "exx_q_pack_" << MY_POOL;
		ofstream ofs(ssc.str().c_str());
    	if (!ofs)
    	{
        	WARNING("Exx_Lip::write_q_pack","Can't create Exx_Lip File!");
    	}

		ofs<<q_pack->kv_ptr->nks<<endl;
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			ofs<<q_pack->kv_ptr->ngk[iq]<<" ";
		}
		ofs<<endl;
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			ofs<<q_pack.kvec_c[iq].x<<" "<<q_pack.kvec_c[iq].y<<" "<<q_pack.kvec_c[iq].z<<" "<<endl;
			ofs<<q_pack.kvec_d[iq].x<<" "<<q_pack.kvec_d[iq].y<<" "<<q_pack.kvec_d[iq].z<<" "<<endl;
		}
		ofs<<endl;
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			for( int ib=0; ib<NBANDS; ++ib)
			{
				ofs<<q_pack.wf_wg[iq][ib]<<" ";
			}
			ofs<<endl;
		}
		ofs<<endl;
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			for( int iw=0; iw<NLOCAL; ++iw)
			{
				for( int ib=0; ib<NBANDS; ++ib)
				{
					ofs<<q_pack.hvec_array[iq](iw,ib).real()<<" "<<q_pack.hvec_array[iq](iw,ib).imag()<<" ";
				}
				ofs<<endl;
			}
		}
		ofs.close();
	}
	return;
}
*/

/*
void Exx_Lip::read_q_pack()
{

	ifstream ifs;
	if( !RANK_IN_POOL )
	{
       	stringstream ssc;
        ssc << global_out_dir << "exx_q_pack_" << MY_POOL;
		ifs.open(ssc.str().c_str());
    	if (!ifs)
    	{
        	WARNING("Exx_Lip::write_q_pack","Can't read Exx_Lip File!");
    	}

		ifs >> q_pack->kv_ptr->nks;
	}

	MPI_Bcast( &q_pack->kv_ptr->nks, 1, MPI_INT, 0, POOL_WORLD);

	q_pack->kv_ptr->ngk = new int [q_pack->kv_ptr->nks];
	q_pack.kvec_c = new Vector3<double> [q_pack->kv_ptr->nks];
	q_pack.kvec_d = new Vector3<double> [q_pack->kv_ptr->nks];
	double *kvec_tmp = new double [q_pack->kv_ptr->nks*6];				// just for MPI
	q_pack.wf_wg = new double *[q_pack->kv_ptr->nks];
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		q_pack.wf_wg[iq] = new double[NBANDS];
	}
	q_pack.hvec_array = new ComplexMatrix [q_pack->kv_ptr->nks];
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		q_pack.hvec_array[iq].create(NLOCAL,NBANDS);
	}

	if( !RANK_IN_POOL )
	{
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			ifs>>q_pack->kv_ptr->ngk[iq];
		}
		for( int iq_tmp = 0; iq_tmp < q_pack->kv_ptr->nks*6; ++iq_tmp)
		{
			ifs>>kvec_tmp[iq_tmp];
		}
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			for( int ib=0; ib<NBANDS; ++ib)
			{
				ifs>>q_pack.wf_wg[iq][ib];
			}
		}
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			for( int iw=0; iw<NLOCAL; ++iw)
			{
				for( int ib=0; ib<NBANDS; ++ib)
				{
					ifs>>q_pack.hvec_array[iq](iw,ib).real()>>q_pack.hvec_array[iq](iw,ib).imag();
				}
			}
		}
		ifs.close();
	}

	MPI_Bcast( q_pack->kv_ptr->ngk, q_pack->kv_ptr->nks, MPI_INT, 0, POOL_WORLD);
	MPI_Bcast( kvec_tmp, q_pack->kv_ptr->nks*6, MPI_DOUBLE, 0, POOL_WORLD);
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		q_pack.kvec_c[iq].x = kvec_tmp[iq*6+0];
		q_pack.kvec_c[iq].y = kvec_tmp[iq*6+1];
		q_pack.kvec_c[iq].z = kvec_tmp[iq*6+2];
		q_pack.kvec_d[iq].x = kvec_tmp[iq*6+3];
		q_pack.kvec_d[iq].y = kvec_tmp[iq*6+4];
		q_pack.kvec_d[iq].z = kvec_tmp[iq*6+5];
	}
	delete[] kvec_tmp;
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		MPI_Bcast( q_pack.wf_wg[iq], NBANDS, MPI_DOUBLE, 0, POOL_WORLD);
	}
	for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
	{
		MPI_Bcast( q_pack.hvec_array[iq].c, NLOCAL*NBANDS, mpicomplex, 0, POOL_WORLD);
	}


	auto test_print = [&]()
	{
		stringstream sss;
		sss << global_out_dir << "exx_q_pack_tmp" << MY_RANK;
		ofstream ofs(sss.str().c_str());
		if (!ofs)
		{
			WARNING("Exx_Lip::write_q_pack","Can't create Exx_Lip File!");
		}

		ofs<<q_pack->kv_ptr->nks<<endl;
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			ofs<<q_pack->kv_ptr->ngk[iq]<<" ";
		}
		ofs<<endl;
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			ofs<<q_pack.kvec_c[iq].x<<" "<<q_pack.kvec_c[iq].y<<" "<<q_pack.kvec_c[iq].z<<" "<<endl;
			ofs<<q_pack.kvec_d[iq].x<<" "<<q_pack.kvec_d[iq].y<<" "<<q_pack.kvec_d[iq].z<<" "<<endl;
		}
		ofs<<endl;
		for( int iq = 0; iq < q_pack->kv_ptr->nks; ++iq)
		{
			for( int ib=0; ib<NBANDS; ++ib)
			{
				ofs<<q_pack.wf_wg[iq][ib]<<" ";
			}
			ofs<<endl;
		}
		ofs<<endl;
		for( int iq=0; iq<q_pack->kv_ptr->nks; ++iq)
		{
			for( int iw=0; iw<NLOCAL; ++iw)
			{
				for( int ib=0; ib<NBANDS; ++ib)
				{
					ofs<<q_pack.hvec_array[iq](iw,ib).real()<<" "<<q_pack.hvec_array[iq](iw,ib).imag()<<" ";
				}
				ofs<<endl;
			}
		}
		ofs.close();
	};

	return;
}
*/


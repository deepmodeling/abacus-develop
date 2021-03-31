#include "exx_abfs-screen-cauchy.h"

#include "src_pw/global.h"

#include "src_external/src_test/src_ri/exx_abfs-screen-test.h"
//double Exx_Abfs::Screen::Cauchy::num_screen1 = 0;
//double Exx_Abfs::Screen::Cauchy::num_screen2 = 0;
//double Exx_Abfs::Screen::Cauchy::num_screen3 = 0;
//double Exx_Abfs::Screen::Cauchy::num_cal = 0;
	
void Exx_Abfs::Screen::Cauchy::init(
	const bool flag_screen_cauchy_in,
	const double threshold_in,
	const Abfs::Vector3_Order<int> Born_von_Karman_period_in)
{
	TITLE("Exx_Abfs::Screen::Cauchy::init");
	flag_screen_cauchy = flag_screen_cauchy_in;
	threshold = threshold_in;
	Born_von_Karman_period = Born_von_Karman_period_in;
}

void Exx_Abfs::Screen::Cauchy::cal_norm_C_max( 
	const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> & Cs,
	const Element_Basis_Index::IndexLNM & index_lcaos,
	const Element_Basis_Index::IndexLNM & index_abfs)
{
	if(!flag_screen_cauchy)	return;
	TITLE("Exx_Abfs::Screen::Cauchy::cal_norm_C_max");
	
	for( const auto & Cs1 : Cs )
	{
		const size_t iat1 = Cs1.first;
		const size_t it1 = ucell.iat2it[iat1];
		for( const auto & Cs2 : Cs1.second )
		{
			const size_t iat2 = Cs2.first;
			const size_t it2 = ucell.iat2it[iat2];
			for( const auto & Cs3 : Cs2.second )
			{
				const Abfs::Vector3_Order<int> &box2 = Cs3.first;
				const matrix & C = *Cs3.second;

				const size_t nw1 = index_lcaos[it1].count_size,
				             nw2 = index_lcaos[it2].count_size,
				             nw3 = index_abfs[it1].count_size;

				{ //    || C(i) ||    || C(i) C^+(i) ||
					valarray<double> C_norm2_outer(nw1);
					valarray<double> C_norm4_outer(nw1);
					for( size_t iw1=0; iw1!=nw1; ++iw1 )
					{
						// C_outer( iw2, iw3 ) = C( iw1, iw2, iw3 )
						matrix C_outer( nw2, nw3 );
//						for( size_t iw2=0; iw2!=nw2; ++iw2 )
//							memcpy( C_outer.c+iw2*nw3, C.c+(iw1*nw2+iw2)*nw3, sizeof(double)*nw3 );
						memcpy( C_outer.c, C.c+iw1*nw2*nw3, sizeof(double)*nw2*nw3 );

						C_norm2_outer[iw1] = C_outer.norm();		// || C ||
						C_norm4_outer[iw1] = m_mT(C_outer).norm();	// || C C^+ ||
					}
					C_norm2_outer_max[iat1][iat2][box2] = C_norm2_outer.max();				// || C ||
					C_norm4_outer_max[iat1][iat2][box2] = sqrt( C_norm4_outer.max() );		// sqrt{ || C C^+ || }
				}

				{ //    || C(j) ||    || C(j) C^+(j) ||
					valarray<double> C_norm2_inner(nw2);
					valarray<double> C_norm4_inner(nw2);
					for( size_t iw2=0; iw2!=nw2; ++iw2 )
					{
						// C_inner( iw1, iw3 ) = C( iw1, iw2, iw3 )
						matrix C_inner( nw1, nw3 );
//						for( size_t iw1=0; iw1!=nw1; ++iw1 )
//							memcpy( C_inner.c+iw1*nw3, C.c+(iw1*nw2+iw2)*nw3, sizeof(double)*nw3 );
						const size_t nw23 = nw2*nw3;
						const auto length = sizeof(double)*nw3;
						double * C_inner_ptr = C_inner.c;
						const double * C_ptr = C.c+iw2*nw3;
						for( size_t iw1=0; iw1!=nw1; ++iw1, C_inner_ptr+=nw3, C_ptr+=nw23 )
							memcpy( C_inner_ptr, C_ptr, length );

						C_norm2_inner[iw2] = C_inner.norm();		// || C ||
						C_norm4_inner[iw2] = m_mT(C_inner).norm();	// || C C^+ ||
					}
					C_norm2_inner_max[iat1][iat2][box2] = C_norm2_inner.max();				// || C ||
					C_norm4_inner_max[iat1][iat2][box2] = sqrt( C_norm4_inner.max() );		// sqrt{ || C C^+ || }
				}

				/*
				matrix C_norm2_12( nw1, nw2 );
				for( size_t iw1=0; iw1!=nw1; ++iw1 )
					for( size_t iw2=0; iw2!=nw2; ++iw2 )
						C_norm2_12(iw1,iw2) = nrm2( nw3, C+iw1+iw2, 1 );
				valarray<double> C_norm2_1(nw1);
				for( size_t iw1=0; iw1!=nw1; ++iw1 )
					C_norm2_1(iw1) = nrm2( nw2, C_norm2_12+iw1, 1 );
				valarray<double> C_norm2_2(nw2);
				for( size_t iw2=0; iw2!=nw2; ++iw2 )
					C_norm2_2(iw2) = nrm2( nw1, C_norm2_12+iw2, 1 );

				C_norm2_1_max[iat1][iat2][box2] = C_norm2_1.max();
				C_norm2_2_max[iat1][iat2][box2] = C_norm2_2.max();
				*/
			}
		}
	}
}

void Exx_Abfs::Screen::Cauchy::cal_norm_V( const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> & Vs )
{
	if(!flag_screen_cauchy)	return;
	TITLE("Exx_Abfs::Screen::Cauchy::cal_norm_V");
	
	for( const auto & Vs1 : Vs )
	{
		const size_t iat1 = Vs1.first;
		for( const auto & Vs2 : Vs1.second )
		{
			const size_t iat2 = Vs2.first;
			for( const auto & Vs3 : Vs2.second )
			{
				const Abfs::Vector3_Order<int> &box2 = Vs3.first;
				const matrix & V = *Vs3.second;

				V_norm4[iat1][iat2][box2] = sqrt( m_mT(V).norm() );		// \sqrt{ || V V^+ || }
			}
		}
	}
}

void Exx_Abfs::Screen::Cauchy::cal_norm_D_max( const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> & Ds )
{
	if(!flag_screen_cauchy)	return;
	TITLE("Exx_Abfs::Screen::Cauchy::cal_norm_D_max");
	
	for( size_t is=0; is!=NSPIN; ++is )
	{
		for( const auto & DsA : Ds[is] )
		{
			const size_t iat1 = DsA.first;
			for( const auto & DsB : DsA.second )
			{
				const size_t iat2 = DsB.first;
				for( const auto & DsC : DsB.second )
				{
					const Abfs::Vector3_Order<int> &box2 = DsC.first;
					const matrix & D = DsC.second;
					const double D_norm4 = sqrt( m_mT(D).norm() );											// \sqrt{ || D D^+ || }
					D_norm4_max[iat1][iat2][box2] = max( D_norm4, D_norm4_max[iat1][iat2][box2] );			// \sqrt{ || D D^+ || }
				}
			}
		}
	}
}

Exx_Abfs::Screen::Cauchy::Info_Step Exx_Abfs::Screen::Cauchy::input_info(
	const size_t iat1, const size_t iat2, const size_t iat3, const size_t iat4,
	const Abfs::Vector3_Order<int> & box2,
	const Abfs::Vector3_Order<int> & box3,
	const Abfs::Vector3_Order<int> & box4 ) const
{
	if(!flag_screen_cauchy)	return Info_Step();

	Info_Step info_step;
	
	info_step.C1_norm4_max = C_norm4_outer_max.at(iat1).at(iat3).at(box3);
	info_step.C3_norm4_max = C_norm4_inner_max.at(iat1).at(iat3).at(box3);
	info_step.C2_norm4_max = C_norm4_outer_max.at(iat2).at(iat4).at(box4);
	info_step.C4_norm4_max = C_norm4_inner_max.at(iat2).at(iat4).at(box4);

	info_step.C1_norm2_max = C_norm2_outer_max.at(iat1).at(iat3).at(box3);
	info_step.C3_norm2_max = C_norm2_inner_max.at(iat1).at(iat3).at(box3);

	info_step.V_norm4 = V_norm4.at(iat1).at(iat2).at(box2);

	if( const double*const D_ptr = static_cast<const double*const>(MAP_EXIST( D_norm4_max, iat3, iat4, Abfs::Vector3_Order<int>(box2-box3+box4)%Born_von_Karman_period )) )
		info_step.D34_norm4_max = *D_ptr;
	else
		info_step.D34_norm4_max = 0;
	if( const double*const D_ptr = static_cast<const double*const>(MAP_EXIST( D_norm4_max, iat3, iat2, Abfs::Vector3_Order<int>(box2-box3     )%Born_von_Karman_period )) )
		info_step.D32_norm4_max = *D_ptr;
	else
		info_step.D32_norm4_max = 0;
	if( const double*const D_ptr = static_cast<const double*const>(MAP_EXIST( D_norm4_max, iat1, iat4, Abfs::Vector3_Order<int>(box2     +box4)%Born_von_Karman_period )) )
		info_step.D14_norm4_max = *D_ptr;
	else
		info_step.D14_norm4_max = 0;
	if( const double*const D_ptr = static_cast<const double*const>(MAP_EXIST( D_norm4_max, iat1, iat2, Abfs::Vector3_Order<int>(box2          )%Born_von_Karman_period )) )
		info_step.D12_norm4_max = *D_ptr;
	else
		info_step.D12_norm4_max = 0;

	return info_step;
}

int Exx_Abfs::Screen::Cauchy::postcalA( const Info_Step & info_step ) const
{
	if( !flag_screen_cauchy )	return 4;
	if( info_step.C1_norm4_max * info_step.V_norm4 * info_step.C2_norm4_max * info_step.D34_norm4_max > threshold )		return 4;
	if( info_step.C3_norm4_max * info_step.V_norm4 * info_step.C2_norm4_max * info_step.D14_norm4_max > threshold )		return 3;
	if( info_step.C1_norm4_max * info_step.V_norm4 * info_step.C4_norm4_max * info_step.D32_norm4_max > threshold )		return 2;
	if( info_step.C3_norm4_max * info_step.V_norm4 * info_step.C4_norm4_max * info_step.D12_norm4_max > threshold )		return 1;
//	++num_screen1;
	return 0;
}

int Exx_Abfs::Screen::Cauchy::postcalB(
	const Info_Step & info_step,
	const matrix & VC_T,			// iw2, \mu1, iw4
	const size_t range_iw2,
	const size_t range_mu1,
	const size_t range_iw4,
	const int postcal_in) const
{
	if( !flag_screen_cauchy )	return 4;
	switch(postcal_in)							// Attention: case and go on calculating
	{
		case 4: case 3:
		{
			const double VC2_norm2_max = cal_matrix_outer_max( VC_T, range_iw2, range_mu1 * range_iw4 );
			if( info_step.C1_norm4_max * VC2_norm2_max * info_step.D34_norm4_max > threshold )		return 4;		// H_{12}
			if( info_step.C3_norm4_max * VC2_norm2_max * info_step.D14_norm4_max > threshold )		return 3;		// H_{32}
		}
		case 2: case 1:
		{
			const double VC4_norm2_max = cal_matrix_inner_max( VC_T, range_iw2 * range_mu1, range_iw4 );
			if( info_step.C1_norm4_max * VC4_norm2_max * info_step.D32_norm4_max > threshold )		return 2;		// H_{14}
			if( info_step.C3_norm4_max * VC4_norm2_max * info_step.D12_norm4_max > threshold )		return 1;		// H_{34}
		}
	}
//	++num_screen2;
	return 0;
}

bool Exx_Abfs::Screen::Cauchy::postcalC(
	const Info_Step & info_step,
	const matrix & DVC,				// iw1/iw3, \mu1, iw2/iw4
	const size_t range_iw13,
	const size_t range_mu1,
	const size_t range_iw24,
	const size_t iw13H) const
{
	if( !flag_screen_cauchy )	return true;
	const double DVC24_norm2_max = cal_matrix_inner_max( DVC, range_iw13 * range_mu1, range_iw24 );
	switch(iw13H)
	{
		case 1:
			if( info_step.C1_norm2_max * DVC24_norm2_max > threshold )
			{
//				num_cal += 0.25/NSPIN;
				return true;
			}
			else 
			{
//				num_screen3 += 0.25/NSPIN; 
				return false; 
			}
			break;
		case 3:
			if( info_step.C3_norm2_max * DVC24_norm2_max > threshold )	
			{ 
//				num_cal += 0.25/NSPIN; 
				return true; 
			}
			else
			{ 
//				num_screen3 += 0.25/NSPIN; 
				return false; 
			}	
			break;
		default:
			throw domain_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}
}

// max_j \sqrt{ \sum_i m(i,j)^2 }
double Exx_Abfs::Screen::Cauchy::cal_matrix_inner_max( const matrix & m, const size_t ni, const size_t nj ) const
{
	assert( m.nr*m.nc == ni*nj );
	valarray<double> m_inner(0.0,nj);
	double * const m_inner_ptr = &m_inner[0];
	for( size_t i=0; i<ni; ++i )
	{
		const double * const m_iw1 = m.c + i*nj;
		#pragma ivdep
		for( size_t j=0; j<nj; ++j )
			m_inner_ptr[j] += m_iw1[j] * m_iw1[j];
	}
	return sqrt(m_inner.max());
}

// max_i \sqrt{ \sum_j m(i,j)^2 }
double Exx_Abfs::Screen::Cauchy::cal_matrix_outer_max( const matrix & m, const size_t ni, const size_t nj ) const
{
	assert( m.nr*m.nc == ni*nj );
	valarray<double> m_outer(ni);
	for( size_t i=0; i<ni; ++i )
		m_outer[i] = LapackConnector::nrm2( nj, m.c+i*nj, 1 );
	return m_outer.max();
}

// m m^+
matrix Exx_Abfs::Screen::Cauchy::m_mT( const matrix & m ) const
{
	matrix mm( m.nr, m.nr, false );
	LapackConnector::gemm( 'N','T', m.nr,m.nr,m.nc, 1, m.c,m.nc, m.c,m.nc, 0, mm.c,mm.nc );
	return mm;
}
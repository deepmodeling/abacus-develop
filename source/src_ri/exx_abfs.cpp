#include "exx_abfs.h"

#include "exx_abfs-abfs_index.h"
#include "exx_abfs-jle.h"
#include "exx_abfs-inverse_matrix_double.h"
#include "exx_abfs-io.h"
#include "exx_abfs-construct_orbs.h"

#include "exx_abfs-matrix_orbs11.h"
#include "exx_abfs-matrix_orbs21.h"
#include "exx_abfs-matrix_lcaoslcaos_lcaoslcaos.h"

#include "../module_orbital/ORB_read.h"
#include "conv_coulomb_pot.h"
#include "conv_coulomb_pot-inl.h"

#include "../module_base/global_function.h"
#include "../src_pw/global.h"
#include<sys/time.h>				// Peize Lin test

int Exx_Abfs::Lmax = 0;		// Peize Lin test

// <a|a> .I
// &
// <a|a> <a|b>
// <b|a> <b|b> .I
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> Exx_Abfs::cal_I(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms,
	const ModuleBase::Element_Basis_Index::IndexLNM &index )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> ms_I;

	for( const auto &m1 : ms )
	{
		const size_t TA = m1.first;
		for( const auto &m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto &m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto &m4 : m3.second )
				{
					const size_t IB = m4.first;

					Inverse_Matrix_Double ms_tmp;

					if( TA==TB && IA==IB )
					{
						const size_t size_A = index[TA].count_size;

						ms_tmp.init( size_A );
						ms_tmp.input(
							ms.at(TA).at(IA).at(TA).at(IA) );
					}
					else
					{
						const size_t size_A = index[TA].count_size;
						const size_t size_B = index[TB].count_size;

						ms_tmp.init( size_A + size_B );
						ms_tmp.input(
							ms.at(TA).at(IA).at(TA).at(IA),
							ms.at(TA).at(IA).at(TB).at(IB),
							ms.at(TB).at(IB).at(TA).at(IA),
							ms.at(TB).at(IB).at(TB).at(IB));
					}

					// Peize Lin test
//					std::cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<std::endl;
//					const ModuleBase::matrix matrix_origin( ms_tmp.A );
//					std::cout<<matrix_origin<<std::endl;

					ms_tmp.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dpotrf );

					// Peize Lin test
//					std::cout<<ms_tmp.A<<std::endl;
//					std::cout<<ms_tmp.A * matrix_origin <<std::endl;
//					std::cout<<matrix_origin * ms_tmp.A<<std::endl;

					if( TA==TB && IA==IB )
					{
						const size_t size_A = index[TA].count_size;

						ms_I[TA][IA][TB][IB].resize( 1, std::vector<ModuleBase::matrix>(1) );
						ms_I[TA][IA][TB][IB][0][0].create( size_A, size_A );

						ms_tmp.output(
							ms_I[TA][IA][TB][IB][0][0]);
					}
					else
					{
						const size_t size_A = index[TA].count_size;
						const size_t size_B = index[TB].count_size;

						ms_I[TA][IA][TB][IB].resize( 2, std::vector<ModuleBase::matrix>(2) );
						ms_I[TA][IA][TB][IB][0][0].create( size_A, size_A );
						ms_I[TA][IA][TB][IB][0][1].create( size_A, size_B );
						ms_I[TA][IA][TB][IB][1][0].create( size_B, size_A );
						ms_I[TA][IA][TB][IB][1][1].create( size_B, size_B );

						ms_tmp.output(
							ms_I[TA][IA][TB][IB][0][0],
							ms_I[TA][IA][TB][IB][0][1],
							ms_I[TA][IA][TB][IB][1][0],
							ms_I[TA][IA][TB][IB][1][1]);
					}
				}
			}
		}
	}
	return ms_I;
}

// <ij|P> * <P|P>.I
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> Exx_Abfs::cal_C(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_abfs,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_abfs_abfs_I )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> ms_C;

	for( const auto & m1 : ms_lcaos2_abfs )
	{
		const size_t TA = m1.first;
		for( const auto & m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto &m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto &m4 : m3.second )
				{
					const size_t IB = m4.first;

					if( TA==TB && IA==IB )
					{
						ms_C[TA][IA][TB][IB].resize(1);
						ms_C[TA][IA][TB][IB][0].create( m4.second[0].nr, m4.second[0].nc );

						const auto &m_lcaos2_abfs = m4.second;
						const auto &m_abfs_abfs_I = ms_abfs_abfs_I.at(TA).at(IA).at(TB).at(IB);

						ms_C[TA][IA][TB][IB][0] = m_lcaos2_abfs[0] * m_abfs_abfs_I[0][0];
					}
					else
					{
						ms_C[TA][IA][TB][IB].resize(2);
						for( size_t i=0; i<2; ++i )
							ms_C[TA][IA][TB][IB][i].create( m4.second[i].nr, m4.second[i].nc );

						const auto &m_lcaos2_abfs = m4.second;
						const auto &m_abfs_abfs_I = ms_abfs_abfs_I.at(TA).at(IA).at(TB).at(IB);
						ms_C[TA][IA][TB][IB][0] = m_lcaos2_abfs[0] * m_abfs_abfs_I[0][0]
													+ m_lcaos2_abfs[1] * m_abfs_abfs_I[1][0];
						ms_C[TA][IA][TB][IB][1] = m_lcaos2_abfs[0] * m_abfs_abfs_I[0][1]
													+ m_lcaos2_abfs[1] * m_abfs_abfs_I[1][1];
					}
				}
			}
		}
	}
	return ms_C;
}

// (<ij|P>*<P|P>.I) * <P|P> * (<P|P>.I*<P|kl>)
void Exx_Abfs::cal_CVC(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_C,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms_abfs_abfs ) const
{
	std::ofstream ofs("ms_CVC");

	for( const auto & m11 : ms_C )
	{
		const size_t TA = m11.first;
		for( const auto & m12 : m11.second )
		{
			const size_t IA = m12.first;
			for( const auto & m13 : m12.second )
			{
				const size_t TB = m13.first;
				for( const auto & m14 : m13.second )
				{
					const size_t IB = m14.first;
					for( const auto & m21 : ms_C )
					{
						const size_t TC = m21.first;
						for( const auto & m22 : m21.second )
						{
							const size_t IC = m22.first;
							for( const auto & m23 : m22.second )
							{
								const size_t TD = m23.first;
								for( const auto & m24 : m23.second )
								{
									const size_t ID = m24.first;

									auto m_00 = [&](){ return m14.second[0] * ms_abfs_abfs.at(TA).at(IA).at(TC).at(IC) * transpose(m24.second[0]); };
									auto m_01 = [&](){ return m14.second[0] * ms_abfs_abfs.at(TA).at(IA).at(TD).at(ID) * transpose(m24.second[1]); };
									auto m_10 = [&](){ return m14.second[1] * ms_abfs_abfs.at(TB).at(IB).at(TC).at(IC) * transpose(m24.second[0]); };
									auto m_11 = [&](){ return m14.second[1] * ms_abfs_abfs.at(TB).at(IB).at(TD).at(ID) * transpose(m24.second[1]); };

									//matrix_CVC[TA][IA][TB][IB]|[TC][IC][TD][ID] =
									ModuleBase::matrix mm;		// Peize Lin test

									if( TA==TB && IA==IB )
									{
										if( TC==TD && IC==ID )
										{
											mm = m_00();
										}
										else
										{
											mm = m_00()+m_01();
										}
									}
									else
									{
										if( TC==TD && IC==ID )
										{
											mm = m_00()+m_10();
										}
										else
										{
											mm = m_00()+m_01()+m_10()+m_11();
										}
									}

									// Peize Lin test
									ofs<<IA<<"\t"<<IB<<"\t"<<IC<<"\t"<<ID<<std::endl;
									mm.print(ofs, 1E-10)<<std::endl;
								}
							}
						}
					}
				}
			}
		}
	}
	ofs.close();
}

// <ij|P> * <P|P>.I * <P|ij>
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> Exx_Abfs::cal_lcaos2_lcaos2_proj_asa(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_asa,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_asa_asa_I,
	const ModuleBase::Element_Basis_Index::Range &range,
	const ModuleBase::Element_Basis_Index::IndexLNM &index)
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> ms_lcaos2_lcaos2_proj_asa;
	for( const auto &m1 : ms_lcaos2_asa )
	{
		const size_t TA = m1.first;
		for( const auto & m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto & m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto & m4 : m3.second )
				{
					const size_t IB = m4.first;
					const auto & m_lcaos2_asa = m4.second;
					const auto & m_abfs_abfs_I = ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB);

					std::vector<ModuleBase::matrix> mql(2);
					size_t matrix_num;
					if( TA==TB && IA==IB )
					{
						matrix_num = 1;
						mql[0] = m_lcaos2_asa[0] * m_abfs_abfs_I[0][0];
					}
					else
					{
						matrix_num = 2;
						mql[0] = m_lcaos2_asa[0] * m_abfs_abfs_I[0][0] + m_lcaos2_asa[1] * m_abfs_abfs_I[1][0];
						mql[1] = m_lcaos2_asa[0] * m_abfs_abfs_I[0][1] + m_lcaos2_asa[1] * m_abfs_abfs_I[1][1];
					}

					ms_lcaos2_lcaos2_proj_asa[TA][IA][TB][IB].create(
						index[TA].count_size,
						index[TB].count_size );
					for( size_t LA=0; LA!=range[TA].size(); ++LA )
						for( size_t NA=0; NA!=range[TA][LA].N; ++NA )
							for( size_t MA=0; MA!=range[TA][LA].M; ++MA )
								for( size_t LB=0; LB!=range[TB].size(); ++LB )
									for( size_t NB=0; NB!=range[TB][LB].N; ++NB )
										for( size_t MB=0; MB!=range[TB][LB].M; ++MB )
										{
											double mv_sas = 0.0;
											for( size_t im=0; im!=matrix_num; ++im )
											{
												assert( mql[im].nc==m_lcaos2_asa[im].nc );
												for( size_t ic=0; ic!=mql[im].nc; ++ic )
												{
													mv_sas
													+= mql[im](
														Exx_Abfs::Abfs_Index::get_index_index(
															index,TA,LA,NA,MA,
															index,TB,LB,NB,MB ),
														ic)
													* m_lcaos2_asa[im](
														Exx_Abfs::Abfs_Index::get_index_index(
															index,TA,LA,NA,MA,
															index,TB,LB,NB,MB ),
														ic);
												}
											}
											ms_lcaos2_lcaos2_proj_asa[TA][IA][TB][IB](
												index[TA][LA][NA][MA],
												index[TB][LB][NB][MB] )
											= mv_sas;
										}
				}
			}
		}
	}
	return ms_lcaos2_lcaos2_proj_asa;
}

// <ij|P> * <P|P>.I * <P|jY>
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> Exx_Abfs::cal_lcaos2_jys_proj_asa(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_asa,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_asa_asa_I,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms_asa_jys)
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> ms_lcaos2_jys_proj_asa;

	for( const auto & m1 : ms_lcaos2_asa )
	{
		const size_t TA = m1.first;
		for( const auto & m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto &m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto &m4 : m3.second )
				{
					const size_t IB = m4.first;

					if( TA==TB && IA==IB )
					{
						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB].resize(1);

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB][0]
						= ms_lcaos2_asa.at(TA).at(IA).at(TB).at(IB)[0]
						* ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB)[0][0]
						* ms_asa_jys.at(TA).at(IA).at(TA).at(IA);
					}
					else
					{
						std::vector<ModuleBase::matrix> ms_tmp(2);
						for( size_t i=0; i<2; ++i )
							ms_tmp[i]
							= ms_lcaos2_asa.at(TA).at(IA).at(TB).at(IB)[0]
							* ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB)[0][i]
							+ ms_lcaos2_asa.at(TA).at(IA).at(TB).at(IB)[1]
							* ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB)[1][i];

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB].resize(2);

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB][0]
						= ms_tmp[0]
						* ms_asa_jys.at(TA).at(IA).at(TA).at(IA)
						+ ms_tmp[1]
						* ms_asa_jys.at(TB).at(IB).at(TA).at(IA);

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB][1]
						= ms_tmp[0]
						* ms_asa_jys.at(TA).at(IA).at(TB).at(IB)
						+ ms_tmp[1]
						* ms_asa_jys.at(TB).at(IB).at(TB).at(IB);
					}
				}
			}
		}
	}
	return ms_lcaos2_jys_proj_asa;
}

/*
void cal_R_supercell()
{
	std::vector<Vector3_Exx> R_supercell;
	for( size_t x=0; x!=GlobalC::kv.nmp[0]; ++x)
		for( size_t y=0; y!=GlobalC::kv.nmp[1]; ++y )
			for( size_t z=0; z!=GlobalC::kv.nmp[2]; ++z )
				R_supercell.push_back(Vector3_Exx(x,y,z));
}
*/
/*
void density_matrix()
{
	std::vector<ModuleBase::matrix> DM_k(GlobalC::kv.nks);
	for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
	{
		for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
		{
			for( size_t iw1=0; iw1!=GlobalV::NLOCAL; ++iw1 )
			{
				for( size_t iw2=0; iw2!=GlobalV::NLOCAL; ++iw2 )
				{
					DM_k[ik](iw1,iw2) += GlobalC::wf.wg(ik,ib) * conj(GlobalC::LOWF.wfc_k_grid[ik][ib][iw1]) * GlobalC::LOWF.wfc_k_grid[ik][ib][iw2];
				}
			}
		}
	}

	std::vector<size_t,std::vector<ModuleBase::matrix>> DM_R( GlobalV::NSPIN, std::vector<ModuleBase::matrix>(R_supercell.size()) );
	for( size_t is=0; is!=GlobalV::NSPIN; ++is )
	{
		const size_t k_start = (GlobalV::NSPIN==1) ? 0 : ((is==0) ? 0 : (GlobalC::kv.nks/2));
		const size_t k_end = (GlobalV::NSPIN==1) ? GlobalC::kv.nks : ((is==0) ? (GlobalC::kv.nks/2) : GlobalC::kv.nks);
		for( size_it iR=0; iR!=R_supercell.size(); ++iR )
		{
			for( size_t ik=k_start; ik!=k_end; ++ik )
			{
				DM_R[is][iR] += exp(ModuleBase::TWO_PI*IMAG_UNIT*dot(GlobalC::kv.kvec_d[ik],R_supercell[iR])) * DM[ik];
			}
		}
	}
}
*/

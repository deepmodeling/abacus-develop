#include <fstream>
#include <stdexcept>

#include "exx_abfs-io.h"
#include "exx_abfs-jle.h"
#include "exx_abfs-abfs_index.h"
#include "src_pw/global.h"
#include "src_lcao/ORB_read.h"
#include "src_global/global_function.h"
#include "src_global/math_integral.h" // mohan add 2021-04-03


vector<vector<vector<Numerical_Orbital_Lm>>> Exx_Abfs::IO::construct_abfs(
	const LCAO_Orbitals &orbs,
	const vector<string> &files_abfs,
	const double kmesh_times )
{
	vector<vector<vector<Numerical_Orbital_Lm>>> abfs( files_abfs.size() );
	for( size_t T=0; T!=files_abfs.size(); ++T )
		abfs[T] = construct_abfs_T( 
			files_abfs[T],
			T,
			static_cast<int>(orbs.get_kmesh() * kmesh_times) | 1,			// Nk must be odd
//			orbs.get_dk() / kmesh_times,
			orbs.get_dk(),								// Peize Lin change 2017-04-16
			orbs.get_dr_uniform() );		
	
	return abfs;
}

vector<vector<vector<Numerical_Orbital_Lm>>> Exx_Abfs::IO::construct_abfs( 
	const vector<vector<vector<Numerical_Orbital_Lm>>> & abfs_pre,
	const LCAO_Orbitals &orbs,
	const vector<string> &files_abfs,
	const double kmesh_times )
{
	vector<vector<vector<Numerical_Orbital_Lm>>> 
		&&abfs = construct_abfs( orbs, files_abfs, kmesh_times );
			
	assert( abfs.size() == abfs_pre.size() );
	for( size_t T=0; T!=abfs.size(); ++T )
	{
		if( abfs[T].size() < abfs_pre[T].size() )
			abfs[T].resize( abfs_pre[T].size() );
		for( size_t L=0; L!=abfs_pre[T].size(); ++L )
		{
			abfs[T][L].insert( abfs[T][L].begin(), abfs_pre[T][L].begin(), abfs_pre[T][L].end() );
		}
	}	
	
	return abfs;
}

vector<vector<Numerical_Orbital_Lm>> Exx_Abfs::IO::construct_abfs_T( 
	const string & file_name,
	const int &T,
	const int &nk,
	const double &dk,
	const double &dr_uniform)
{
	string label;
	size_t L_size;
	map<size_t,size_t> N_size;
	size_t meshr;
	double dr;
	map<size_t,map<size_t,vector<double>>> psis;
	
	/*----------------------
	  1.read abfs
	----------------------*/
	string word;
	
	ifstream ifs( file_name.c_str() );
	if(!ifs)
		throw runtime_error(" Can't find the abfs ORBITAL file.");	
	
	while( ifs.good() )
	{
		ifs >> word;
		
		if( "Element"==word )
		{
			READ_VALUE( ifs, label );
		}
		else if ( "Lmax"==word )
		{
			READ_VALUE( ifs, L_size );
			
			if( L_size>=9 )
			{
				stringstream ss;
				ss<<"Lmax>=9 error in "<<__FILE__<<" line "<<__LINE__;
				throw invalid_argument(ss.str());
			}
		}
		else if ( "Sorbital-->"==word )
		{
			READ_VALUE( ifs, N_size[0] );
		}
		else if ( "Porbital-->"==word )
		{
			READ_VALUE( ifs, N_size[1] );
		}
		else if ( "Dorbital-->"==word )
		{
			READ_VALUE( ifs, N_size[2] );
		}
		else if ( "Forbital-->"==word )
		{
			READ_VALUE( ifs, N_size[3] );
		}
		else if ( "Gorbital-->"==word )
		{
			READ_VALUE( ifs, N_size[4] );
		}
		else if ( "Horbital-->"==word )
		{
			READ_VALUE( ifs, N_size[5] );
		}
		else if ( "Iorbital-->"==word )
		{
			READ_VALUE( ifs, N_size[6] );
		}
		else if ( "Jorbital-->"==word )
		{
			READ_VALUE( ifs, N_size[7] );
		}
		else if ( "Korbital-->"==word )
		{
			READ_VALUE( ifs, N_size[8] );
		}
		else if ( "END"==word )
		{
			break;
		}		
	}
	
	CHECK_NAME(ifs, "Mesh");
	ifs >> meshr;
	
	CHECK_NAME(ifs, "dr");
	ifs >> dr;

	while(ifs.good())
	{
		ifs >> word;
		if(word=="Type")
		{
			string s_L, s_N;
			ifs >> s_L >> s_N;
			
			size_t T,L,N;
			ifs >> T >> L >> N;
			
			psis[L][N].resize(meshr);
			for(int ir=0; ir!=meshr; ir++)
			{
				ifs >> psis[L][N][ir];
			}		
		}

	}
	ifs.close();

	 
	/*----------------------
	  2.check L,N orbital
	----------------------*/
	for( size_t L=0; L<=L_size; ++L )
		if( N_size.find(L) == N_size.end() )
		{
			stringstream ss;
			ss<<"Can't find N of L="<<L<<" in "<<file_name;
			throw domain_error(ss.str());
		}
	for( size_t L=0; L<=L_size; ++L )
		for( size_t N=0; N!=N_size[L]; ++N )
			if( psis.find(L)==psis.end() || psis[L].find(N)==psis[L].end() )
			{
				stringstream ss;
				ss<<"Can't find abf of L="<<L<<" T="<<T<<" in "<<file_name;
				throw domain_error(ss.str());
			}

			
	/*----------------------
	  3.rab, radial
	----------------------*/		
	if(meshr%2==0)	++meshr;
	
	vector<double> rab(meshr);
	vector<double> radial(meshr);
	for( int ir=0; ir!=meshr; ++ir )
	{
		rab[ir] = dr;
		radial[ir] = ir*dr; //mohan 2010-04-19
	}

	
	/*----------------------
	  4.normalize psi
	----------------------*/
	for( size_t L=0; L<=L_size; ++L )
	{	
		for( size_t N=0; N!=N_size[L]; ++N )
		{
			vector<double> psir(meshr);
			vector<double> inner(meshr);
			psis[L][N].resize(meshr);
			for( int ir=0; ir!=meshr; ++ir )
			{
				psir[ir] = psis[L][N][ir] * radial[ir];
				inner[ir] = psir[ir] * psir[ir];
			}
			double unit = 0.0;	
			Integral::Simpson_Integral(meshr, VECTOR_TO_PTR(inner), VECTOR_TO_PTR(rab), unit);
			for( int ir=0; ir!=meshr; ++ir )
			{
				psis[L][N][ir] /= sqrt(unit);
			}
		}
	}

	
	/*----------------------
	  5.construct abfs
	----------------------*/	
	vector<vector<Numerical_Orbital_Lm>> abfs_T;

	abfs_T.resize(L_size+1);
	for( size_t L=0; L<=L_size; ++L )
	{
		abfs_T[L].resize(N_size[L]);
		for( size_t N=0; N!=N_size[L]; ++N )
		{
			abfs_T[L][N].set_orbital_info(
				label,
				T, //type
				L, //angular momentum L
				N, // number of orbitals of this L
				meshr, // number of radial mesh
				VECTOR_TO_PTR(rab),
				VECTOR_TO_PTR(radial),// radial mesh value(a.u.)
				Numerical_Orbital_Lm::Psi_Type::Psi,
				VECTOR_TO_PTR(psis[L][N]), // radial wave function
				nk,
				dk,
				dr_uniform,
				false,
				true);		
		}
	}
	
	return abfs_T;
}

void Exx_Abfs::IO::print_matrix( 
		const string &file_name_prefix, 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &matrixes_Q, 
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &matrixes_S,
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &matrixes_V,
		const Element_Basis_Index::Range &range_jles, 
		const Element_Basis_Index::IndexLNM &index_jles, 
		const Element_Basis_Index::Range &range_lcaos,
		const Element_Basis_Index::IndexLNM &index_lcaos )
{
	auto print_header = [&]( ofstream &ofs, size_t TA, size_t IA, size_t TB, size_t IB )
	{
		ofs << ucell.lat0 << endl;

		ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;
		
		if( TA==TB )
		{
			ofs << 1 << " ntype" << endl;
			ofs << ucell.atoms[TA].label << " label" << endl;
			if( IA==IB )
			{
				ofs << 1 << " na" << endl;
				ofs << ucell.atoms[TA].tau[IA].x << " " 
					<< ucell.atoms[TA].tau[IA].y << " " 
					<< ucell.atoms[TA].tau[IA].z << endl;
			}
			else
			{
				ofs << 2 << " na" << endl;
				ofs << ucell.atoms[TA].tau[IA].x << " " 
					<< ucell.atoms[TA].tau[IA].y << " "
					<< ucell.atoms[TA].tau[IA].z << endl;
				ofs << ucell.atoms[TB].tau[IB].x << " " 
					<< ucell.atoms[TB].tau[IB].y << " " 
					<< ucell.atoms[TB].tau[IB].z << endl;
			}
		}
		else
		{
			ofs << 2 << " ntype" << endl;
			ofs << ucell.atoms[TA].label << " label" << endl;
			ofs << 1 << " na" << endl;
			ofs << ucell.atoms[TA].tau[IA].x << " " 
				<< ucell.atoms[TA].tau[IA].y << " " 
				<< ucell.atoms[TA].tau[IA].z << endl;
			ofs << ucell.atoms[TB].label << " label" << endl;
			ofs << 1 << " na" << endl;
			ofs << ucell.atoms[TB].tau[IB].x << " " 
				<< ucell.atoms[TB].tau[IB].y << " " 
				<< ucell.atoms[TB].tau[IB].z << endl;
		}
		
		// ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
		ofs << Exx_Abfs::Jle::Ecut_exx << " ecutwfc" << endl; // mohan add 2009-09-08

		// this parameter determine the total number of jlq.
		ofs << Exx_Abfs::Jle::Ecut_exx << " ecutwfc_jlq" << endl;//mohan modify 2009-09-08

		assert( ORB.Phi[TA].getRcut() == ORB.Phi[TB].getRcut() );
		ofs << ORB.Phi[TA].getRcut() << " rcut_Jlq" << endl;

		// mohan add 'smooth' and 'sigma' 2009-08-28
		ofs << 0 << " smooth" << endl;
		ofs << 0 << " sigma" << endl;

		ofs << Exx_Abfs::Jle::tolerence << " tolerence" << endl;

		ofs << Exx_Abfs::Jle::Lmax << " lmax" << endl;

		ofs << kv.nkstot << " nks" << endl;
		ofs	<< index_lcaos[TA].count_size * index_lcaos[TB].count_size << " nbands" << endl;
		
		auto cal_sum_M = [&range_jles](size_t T) -> size_t
		{
			size_t sum_M = 0;
			for( size_t L = 0; L!=range_jles[T].size(); ++L )
				sum_M += range_jles[T][L].M;
			return sum_M;
		};
		const size_t nwfc = (TA==TB && IA==IB) ? cal_sum_M(TA) : cal_sum_M(TA)+cal_sum_M(TB);
		ofs	<< nwfc << " nwfc" << endl;
		
		const size_t ecut_numberA = static_cast<size_t>( sqrt( Exx_Abfs::Jle::Ecut_exx ) * ORB.Phi[TA].getRcut() / PI ); // Rydberg Unit
		const size_t ecut_numberB = static_cast<size_t>( sqrt( Exx_Abfs::Jle::Ecut_exx ) * ORB.Phi[TB].getRcut() / PI ); // Rydberg Unit
		assert( ecut_numberA == ecut_numberB );
		ofs	<< ecut_numberA << " ne" << endl;
		
		ofs << "<WEIGHT_OF_KPOINTS>" << endl;
		for( int ik=0; ik!=kv.nkstot; ++ik )		
		{
			ofs << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z;
			ofs << " " << kv.wk[ik] * 0.5 << endl;
		}
		ofs << "</WEIGHT_OF_KPOINTS>" << endl;

		ofs << endl;
	};
	
	auto print_Q = [&]( ofstream &ofs, const size_t TA, const size_t IA, const size_t TB, const size_t IB, const double scale=1 )
	{
		/*---------------------
		  < jY | Psi >
		---------------------*/	

		ofs<< "<OVERLAP_Q>" << std::endl;
		
		for( int ik=0; ik!=kv.nkstot; ++ik )
			for( size_t LA=0; LA!=range_lcaos[TA].size(); ++LA )
				for( size_t NA=0; NA!=range_lcaos[TA][LA].N; ++NA )
					for( size_t MA=0; MA!=range_lcaos[TA][LA].M; ++MA )
						for( size_t LB=0; LB!=range_lcaos[TB].size(); ++LB )
							for( size_t NB=0; NB!=range_lcaos[TB][LB].N; ++NB )	
								for( size_t MB=0; MB!=range_lcaos[TB][LB].M; ++MB )
								{
									const vector<matrix> & matrix_Q = matrixes_Q.at(TA).at(IA).at(TB).at(IB);
									const size_t index_lcao 
										= Exx_Abfs::Abfs_Index::get_index_index( 
											index_lcaos,TA,LA,NA,MA, 
											index_lcaos,TB,LB,NB,MB );
											
									auto f = [&]( size_t Tj, size_t Ij )
									{
										for( size_t Lj=0; Lj!=range_jles[Tj].size(); ++Lj )
											for( size_t Mj=0; Mj!=range_jles[Tj][Lj].M; ++Mj )				// attention! M first N second
												for( size_t Nj=0; Nj!=range_jles[Tj][Lj].N; ++Nj )
													ofs << matrix_Q[Ij]( 
															index_lcao, 
															index_jles[Tj][Lj][Nj][Mj] )
														 * scale << "\t" 
														<< 0.0 * scale << std::endl;
									};

									if( TA==TB && IA==IB )
									{
										f(TA,0);
									}
									else
									{
										f(TA,0);
										f(TB,1);
									}
								}

		
		ofs<< "</OVERLAP_Q>" << std::endl << std::endl;
	};


	auto print_S = [&]( ofstream &ofs, const size_t TA, const size_t IA, const size_t TB, const size_t IB, const double scale=1 )
	{
		/*---------------------
		  < jY | jY >
		---------------------*/		
		auto q1q2 = [&]( size_t T1, size_t I1, size_t T2, size_t I2, size_t L1, size_t M1, size_t L2, size_t M2)
		{
			const matrix & matrix_S = matrixes_S.at(T1).at(I1).at(T2).at(I2);
			for( size_t N1=0; N1!=range_jles[T1][L1].N; ++N1)
				for( size_t N2=0; N2!=range_jles[T2][L2].N; ++N2)
					ofs<< matrix_S( index_jles[T1][L1][N1][M1], index_jles[T2][L2][N2][M2] ) * scale << "\t" << 0.0 * scale << std::endl;
		};
		
		auto atom2 = [&]( size_t T1, size_t I1, size_t T2, size_t I2, size_t L1, size_t M1)
		{
			for( size_t L2=0; L2!=range_jles[T2].size(); ++L2 )	
				for( size_t M2=0; M2!=range_jles[T2][L2].M; ++M2 )
					q1q2(T1,I1,T2,I2,L1,M1,L2,M2);
		};
		
		auto atom1 = [&]( size_t T1, size_t I1)
		{
			for( size_t L1=0; L1!=range_jles[T1].size(); ++L1 )
				for( size_t M1=0; M1!=range_jles[T1][L1].M; ++M1 )
				{
					if( TA==TB && IA==IB )
					{
						atom2(T1,I1,TA,IA,L1,M1);
					}
					else
					{
						atom2(T1,I1,TA,IA,L1,M1);
						atom2(T1,I1,TB,IB,L1,M1);
					}
				}
		};	
		
		ofs<< "<OVERLAP_Sq>" <<endl;
		
		for( int ik=0; ik!=kv.nkstot; ++ik )
		{
			if( TA==TB && IA==IB )
			{
				atom1(TA,IA);
			}
			else
			{
				atom1(TA,IA);
				atom1(TB,IB);
			}
		}	
		
		ofs<< "</OVERLAP_Sq>" << std::endl << std::endl;
	};


	auto print_V = [&]( ofstream &ofs, const size_t TA, const size_t IA, const size_t TB, const size_t IB, const double scale=1 )
	{
		/*---------------------
		  < Psi | Psi >
		---------------------*/
		ofs << "<OVERLAP_V>" << std::endl;
		
		const matrix & matrix_V = matrixes_V.at(TA).at(IA).at(TB).at(IB);
		
		for( int ik=0; ik!=kv.nkstot; ++ik )
	//		for ib = 0 to NBANDS
			for( size_t LA=0; LA!=range_lcaos[TA].size(); ++LA )
				for( size_t NA=0; NA!=range_lcaos[TA][LA].N; ++NA )
					for( size_t MA=0; MA!=range_lcaos[TA][LA].M; ++MA )
						for( size_t LB=0; LB!=range_lcaos[TB].size(); ++LB )
							for( size_t NB=0; NB!=range_lcaos[TB][LB].N; ++NB )
								for( size_t MB=0; MB!=range_lcaos[TB][LB].M; ++MB )
									ofs << matrix_V( index_lcaos[TA][LA][NA][MA], index_lcaos[TB][LB][NB][MB] ) * scale << std::endl;
		
		ofs << "</OVERLAP_V>" << std::endl << std::endl;
	};

	
	auto cal_R2 = []( const size_t TA, const size_t IA, const size_t TB, const size_t IB ) ->double	
	{
		Numerical_Orbital::set_position( ucell.atoms[TA].tau[IA], ucell.atoms[TB].tau[IB] );
		const double R = Numerical_Orbital::get_distance()*ucell.lat0;
		return R*R;
	};
	
	for( size_t TA=0; TA!=ucell.ntype; ++TA )
	{
		for( size_t IA=0; IA!=ucell.atoms[TA].na; ++IA )
		{
			for( size_t TB=TA; TB!=ucell.ntype; ++TB )
			{
				for( size_t IB=((TB==TA)?IA:0); IB!=ucell.atoms[TB].na; ++IB )
				{
					ofstream ofs(( file_name_prefix+"matrix_"+TO_STRING(TA)+"_"+TO_STRING(IA)+"_"+TO_STRING(TB)+"_"+TO_STRING(IB) ).c_str());
					print_header( ofs, TA, IA, TB, IB );
//					const double scale = 1.0 / max( matrixes_V.at(TA).at(IA).at(TB).at(IB) );
					const double scale = 1;		// Peize Lin test
//					const double scale = cal_R2(TA,IA,TB,IB);
					print_Q( ofs, TA, IA, TB, IB, sqrt(scale) );
					print_S( ofs, TA, IA, TB, IB );
					print_V( ofs, TA, IA, TB, IB, scale );
				}
			}
		}
	}
}

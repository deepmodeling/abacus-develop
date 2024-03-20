#include "exx_abfs-construct_orbs.h"

#include "ABFs_Construct-PCA.h"
#include "module_base/gram_schmidt_orth-inl.h"
#include "module_base/gram_schmidt_orth.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"             //for ucell
#include "module_ri/test_code/exx_abfs-construct_orbs-test.h" // Peize Lin test

std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::Construct_Orbs::change_orbs(
	const LCAO_Orbitals &orbs_in,
	const double kmesh_times )
{
	ModuleBase::TITLE("Exx_Abfs::Construct_Orbs::change_orbs");

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs;
	orbs.resize( orbs_in.get_ntype() );
	for (int T = 0;  T < orbs_in.get_ntype() ; T++)
	{
		orbs[T].resize( orbs_in.Phi[T].getLmax()+1 );
		for (int L=0; L <= orbs_in.Phi[T].getLmax() ; L++)
		{
			orbs[T][L].resize( orbs_in.Phi[T].getNchi(L) );
			for (int N = 0; N < orbs_in.Phi[T].getNchi(L); ++N)
			{
				const auto &orb_origin = orbs_in.Phi[T].PhiLN(L,N);
				orbs[T][L][N].set_orbital_info(
					orb_origin.getLabel(),
					orb_origin.getType(),
					orb_origin.getL(),
					orb_origin.getChi(),
					orb_origin.getNr(),
					orb_origin.getRab(),
					orb_origin.getRadial(),
					Numerical_Orbital_Lm::Psi_Type::Psi,
					orb_origin.getPsi(),
					static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,		// Nk must be odd
					orb_origin.getDk(),							// Peize Lin change 2017-04-16
//					orb_origin.getDk() / kmesh_times,
					orb_origin.getDruniform(),
					false,
					true, GlobalV::CAL_FORCE);
			}
		}
	}
	return orbs;
}

std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::Construct_Orbs::change_orbs(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs_in,
	const double kmesh_times )
{
	ModuleBase::TITLE("Exx_Abfs::Construct_Orbs::change_orbs");
	return orbital( get_psi(orbs_in), orbs_in, kmesh_times );
}

// P = u/r * Y
/*
template<typename Orbs_Type>
std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::Construct_Orbs::abfs_same_atom(
	const Orbs_Type &orbs,
	const double kmesh_times,
	const double norm_threshold )
{
	const std::vector<std::vector<std::vector<std::vector<double>>>>
		abfs_same_atom_psir = psir_mult_psir( orbs );
	const std::vector<std::vector<std::vector<std::vector<double>>>>
		abfs_same_atom_psir_orth = orth( abfs_same_atom_psir, orbs, norm_threshold );
	const std::vector<std::vector<std::vector<std::vector<double>>>>
		abfs_same_atom_psi_orth = div_r( abfs_same_atom_psir_orth, orbs.get_r_radial );
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = orbital( abfs_same_atom_psi_orth, orbs, kmesh_times );
	return abfs_same_atom;
}
*/

// P = f * Y
std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::Construct_Orbs::abfs_same_atom(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs,
	const double kmesh_times_mot,
	const double times_threshold )
{
	ModuleBase::TITLE("Exx_Abfs::Construct_Orbs::abfs_same_atom");

	const std::vector<std::vector<std::vector<std::vector<double>>>>
		abfs_same_atom_psi = psi_mult_psi( orbs );

	const std::vector<std::vector<std::vector<std::vector<double>>>>
		abfs_same_atom_orth_psi = orth( abfs_same_atom_psi, orbs );
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = orbital( abfs_same_atom_orth_psi, orbs, 1 );

	#if TEST_EXX_LCAO==1
		print_orbs(abfs_same_atom_psi,"abfs_same_atom_psi.dat");
		print_orbs(abfs_same_atom_orth_psi,"abfs_same_atom_orth_psi.dat");
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

	const std::vector<std::vector<std::vector<std::vector<double>>>>
		abfs_same_atom_pca_psi = pca( abfs_same_atom, orbs, kmesh_times_mot, times_threshold );

	#if TEST_EXX_LCAO==1
		print_orbs(abfs_same_atom_pca_psi,"abfs_same_atom_pca_psi.dat");
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom_pca = orbital( abfs_same_atom_pca_psi, orbs, 1 );
	return abfs_same_atom_pca;
}

/*
template<>
std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::psi_mult_psi(
	const LCAO_Orbitals &orbs )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> psi_mult_psi( orbs.get_ntype() );
	for( int T=0; T!=orbs.get_ntype(); ++T )
	{
		psi_mult_psi[T].resize( 2*orbs.Phi[T].getLmax()+1 );
		for( int L1=0; L1<=orbs.Phi[T].getLmax(); ++L1 )
		{
			for( int N1=0; N1!=orbs.Phi[T].getNchi(L1); ++N1 )
			{
				for( int L2=L1; L2<=orbs.Phi[T].getLmax(); ++L2 )
				{
					for( int N2=((L2==L1)?N1:0); N2!=orbs.Phi[T].getNchi(L2); ++N2 )
					{
						assert( orbs.Phi[T].PhiLN(L1,N1).getNr()==orbs.Phi[T].PhiLN(L2,N2).getNr() );

						std::vector<double> mult_psir( orbs.Phi[T].PhiLN(L1,N1).getNr() );
						for( int ir=0; ir!=orbs.Phi[T].PhiLN(L1,N1).getNr(); ++ir)
						{
							mult_psir[ir] = orbs.Phi[T].PhiLN(L1,N1).getPsi(ir) * orbs.Phi[T].PhiLN(L2,N2).getPsi(ir) ;
						}
						for( int L_new=std::abs(L2-L1); L_new<=L1+L2; ++L_new )
						{
							psi_mult_psi[T][L_new].push_back(mult_psir);
						}
					}
				}
			}
		}
	}
	return psi_mult_psi;
}
*/

std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::psi_mult_psi(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> psi_mult_psi( orbs.size() );
	for( int T=0; T!=orbs.size(); ++T )
	{
		psi_mult_psi[T].resize( 2*orbs[T].size()-1 );
		for( int L1=0; L1!=orbs[T].size(); ++L1 )
		{
			for( int N1=0; N1!=orbs[T][L1].size(); ++N1 )
			{
				for( int L2=L1; L2!=orbs[T].size(); ++L2 )
				{
					for( int N2=((L2==L1)?N1:0); N2!=orbs[T][L2].size(); ++N2 )
					{
						assert( orbs[T][L1][N1].getNr()==orbs[T][L2][N2].getNr() );

						std::vector<double> mult_psir( orbs[T][L1][N1].getNr() );
						for( int ir=0; ir!=orbs[T][L1][N1].getNr(); ++ir)
						{
							mult_psir[ir] = orbs[T][L1][N1].getPsi(ir) * orbs[T][L2][N2].getPsi(ir) ;
						}
						for( int L_new=std::abs(L1-L2); L_new<=L1+L2; ++L_new )
						{
							psi_mult_psi[T][L_new].push_back(mult_psir);
						}
					}
				}
			}
		}
	}
	return psi_mult_psi;
}

/*
template<>
std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::psir_mult_psir(
	const LCAO_Orbitals &orbs )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> psir_mult_psir( orbs.get_ntype() );
	for( int T=0; T!=orbs.get_ntype(); ++T )
	{
		psir_mult_psir[T].resize( 2*orbs.Phi[T].getLmax()+1 );
		for( int L1=0; L1<=orbs.Phi[T].getLmax(); ++L1 )
		{
			for( int N1=0; N1!=orbs.Phi[T].getNchi(L1); ++N1 )
			{
				for( int L2=L1; L2<=orbs.Phi[T].getLmax(); ++L2 )
				{
					for( int N2=((L2==L1)?N1:0); N2!=orbs.Phi[T].getNchi(L2); ++N2 )
					{
						assert( orbs.Phi[T].PhiLN(L1,N1).getNr()==orbs.Phi[T].PhiLN(L2,N2).getNr() );

						std::vector<double> mult_psir( orbs.Phi[T].PhiLN(L1,N1).getNr() );
						for( int ir=0; ir!=orbs.Phi[T].PhiLN(L1,N1).getNr(); ++ir)
						{
							mult_psir[ir] = orbs.Phi[T].PhiLN(L1,N1).getPsi_r(ir) * orbs.Phi[T].PhiLN(L2,N2).getPsi_r(ir) ;
						}
						for( int L_new=std::abs(L1-L2); L_new<=L1+L2; ++L_new )
						{
							psir_mult_psir[T][L_new].push_back(mult_psir);
						}
					}
				}
			}
		}
	}
	return psir_mult_psir;
}
*/

std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::psir_mult_psir(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> psir_mult_psir( orbs.size() );
	for( int T=0; T!=orbs.size(); ++T )
	{
		psir_mult_psir[T].resize( 2*orbs[T].size()-1 );
		for( int L1=0; L1!=orbs[T].size(); ++L1 )
		{
			for( int N1=0; N1!=orbs[T][L1].size(); ++N1 )
			{
				for( int L2=L1; L2!=orbs[T].size(); ++L2 )
				{
					for( int N2=((L2==L1)?N1:0); N2!=orbs[T][L2].size(); ++N2 )
					{
						assert( orbs[T][L1][N1].getNr()==orbs[T][L2][N2].getNr() );

						std::vector<double> mult_psir( orbs[T][L1][N1].getNr() );
						for( int ir=0; ir!=orbs[T][L1][N1].getNr(); ++ir)
						{
							mult_psir[ir] = orbs[T][L1][N1].getPsi_r(ir) * orbs[T][L2][N2].getPsi_r(ir) ;
						}
						for( int L_new=std::abs(L1-L2); L_new<=L1+L2; ++L_new )
						{
							psir_mult_psir[T][L_new].push_back(mult_psir);
						}
					}
				}
			}
		}
	}
	return psir_mult_psir;
}

std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::pca(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs,
	const double kmesh_times_mot,
	const double times_threshold )
{
	if(times_threshold>1)
		return std::vector<std::vector<std::vector<std::vector<double>>>>(abfs.size());

	const std::vector<std::vector<std::pair<std::vector<double>,RI::Tensor<double>>>>
		eig = ABFs_Construct::PCA::cal_PCA( orbs, abfs, kmesh_times_mot );

	const std::vector<std::vector<std::vector<std::vector<double>>>> psis = get_psi( abfs );
	std::vector<std::vector<std::vector<std::vector<double>>>> psis_new( psis.size() );

	for( size_t T=0; T!=eig.size(); ++T )
	{
		double eig_value_max = 0;
		for( size_t L=0; L!=eig[T].size(); ++L )
			for( size_t M=0; M!=eig[T][L].first.size(); ++M )
			{
//ofs<<T<<"\t"<<L<<"\t"<<M<<"\t"<<eig[T][L].first[M]<<std::endl;
				eig_value_max = std::max( eig_value_max, eig[T][L].first[M] );
			}
		const double eig_value_threshold = eig_value_max * times_threshold;

//ofs<<"eig_value_max:\t"<<eig_value_max<<std::endl;
//ofs<<"eig_value_threshold:\t"<<eig_value_threshold<<std::endl;

		if(eig_value_max)
		{
			psis_new[T].resize( psis[T].size() );
			for( size_t L=0; L!=eig[T].size(); ++L )
			{
				const std::vector<double> &eig_value = eig[T][L].first;
				const RI::Tensor<double> &eig_vec = eig[T][L].second;
				for( size_t M=0; M!=eig_value.size(); ++M )
				{
					if( eig_value[M] > eig_value_threshold )
					{
						std::vector<double> psi_new( psis[T][L][0].size() );
						for( size_t N=0; N!=psis[T][L].size(); ++N )
							for( size_t ir=0; ir!=psi_new.size(); ++ir )
								psi_new[ir] += eig_vec(M,N) * psis[T][L][N][ir];
						psis_new[T][L].push_back( psi_new );
					}
				}
			}
		}
		else
		{
			ModuleBase::WARNING(ModuleBase::GlobalFunc::TO_STRING(__FILE__),
				"Element "+ModuleBase::GlobalFunc::TO_STRING(T)+" , all training data (lcao[i]*lcao[j]) are all the same. So PCA randomly choose an abf as the result.");
			psis_new[T].resize( psis[T].size() );
			for( size_t L=0; L!=psis[T].size(); ++L )
				if( !psis[T][L].empty() )
				{
					psis_new[T][L].push_back(psis[T][L][0]);
					break;
				}
		}
	}
//ofs.close();
	return psis_new;
}


std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::orth(
	const std::vector<std::vector<std::vector<std::vector<double>>>> &psis,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs,
	const double norm_threshold )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> psis_orth( psis.size() );
	for( int T=0; T!=psis.size(); ++T )
	{
		const Numerical_Orbital_Lm &orb = orbs[T][0][0];
		ModuleBase::Gram_Schmidt_Orth<double,double> gso(
			orb.get_rab(),
			ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Sphere );
		psis_orth[T].resize( psis[T].size() );
		for( int L=0; L!=psis[T].size(); ++L )
		{
			psis_orth[T][L] = gso.cal_orth( psis[T][L], norm_threshold );		// Peize Lin test 2016-10-08
//			psis_orth[T][L] = gso.cal_orth( psis[T][L] );		// Peize Lin test 2016-10-08
		}
	}
	return psis_orth;
}

std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::div_r(
	const std::vector<std::vector<std::vector<std::vector<double>>>> &psirs,
	const std::vector<double> &r_radial )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> psis( psirs.size() );
	for( auto T=0; T!=psirs.size(); ++T )
	{
		psis[T].resize( psirs[T].size() );
		for( auto L=0; L!=psirs[T].size(); ++L )
		{
			psis[T][L].resize( psirs[T][L].size() );
			for( auto N=0; N!=psirs[T][L].size(); ++N )
			{
				psis[T][L][N].resize( psirs[T][L][N].size() );
				psis[T][L][N][0] = 0;
				for( auto ir=1; ir!=psirs[T][L][N].size(); ++ir )
				{
					psis[T][L][N][ir] = psirs[T][L][N][ir] / r_radial[ir];
				}
			}
		}
	}
	return psis;
}

std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::get_psi(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> orbs_psi( orbs.size() );
	for( int T=0; T!=orbs.size(); ++T )
	{
		orbs_psi[T].resize( orbs[T].size() );
		for( int L=0; L!=orbs[T].size(); ++L )
		{
			orbs_psi[T][L].resize( orbs[T][L].size() );
			for( int N=0; N!=orbs[T][L].size(); ++N )
			{
				orbs_psi[T][L][N] = orbs[T][L][N].get_psi();
			}
		}
	}
	return orbs_psi;
}

std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::Construct_Orbs::orbital(
	const std::vector<std::vector<std::vector<std::vector<double>>>> &psis,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs_info,
	const double kmesh_times)
{
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs_new( psis.size() );
	for( int T=0; T!=psis.size(); ++T )
	{
		const Numerical_Orbital_Lm &orb_info = orbs_info[T][0][0];
		orbs_new[T].resize( psis[T].size() );
		for( int L=0; L!=psis[T].size(); ++L )
		{
			orbs_new[T][L].resize( psis[T][L].size() );
			for( int N=0; N!=psis[T][L].size(); ++N )
			{
				orbs_new[T][L][N].set_orbital_info(
					orb_info.getLabel(),
					T,
					L,
					N,
					orb_info.getNr(),
					orb_info.getRab(),
					orb_info.getRadial(),
					Numerical_Orbital_Lm::Psi_Type::Psi,
					ModuleBase::GlobalFunc::VECTOR_TO_PTR(psis[T][L][N]),
					static_cast<int>(orb_info.getNk() * kmesh_times) | 1,	// Nk must be odd
					orb_info.getDk(),					// Peize Lin test 2017-04-16
//					orb_info.getDk() / kmesh_times,
					orb_info.getDruniform(),
					false,
					true, GlobalV::CAL_FORCE);
			}
		}
	}
	return orbs_new;
}

/*
std::vector<std::vector<std::vector<std::vector<double>>>> Exx_Abfs::Construct_Orbs::get_psi(
	const LCAO_Orbitals &orbs )
{
	std::vector<std::vector<std::vector<std::vector<double>>>> orbs_psi( orbs.get_ntype() );
	for( int T=0; T!=orbs.get_ntype(); ++T )
	{
		orbs_psi[T].resize( orbs.Phi[T].getLmax()+1 );
		for( int L=0; L<=orbs.Phi[T].getLmax(); ++L )
		{
			orbs_psi[T][L].resize( orbs.Phi[T].getNchi(L) );
			for( int N=0; N!=orbs.Phi[T].getNchi(L); ++N )
			{
				orbs_psi[T][L][N] = orbs.Phi[T].PhiLN(L,N).get_psi();
			}
		}
	}
	return orbs_psi;
}
*/

/*
template<>
inline const Numerical_Orbital_Lm &Exx_Abfs::Construct_Orbs::get_orbital(
	const LCAO_Orbitals &orbs,
	const size_t T, const size_t L, const size_t N)
{
	return orbs.Phi[T].PhiLN(L,N);
}
*/

/*
template<>
inline const Numerical_Orbital_Lm &Exx_Abfs::Construct_Orbs::get_orbital(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs,
	const size_t T, const size_t L, const size_t N)
{
	return orbs[T][L][N];
}
*/

void Exx_Abfs::Construct_Orbs::print_orbs_size(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs,
	std::ostream &os)
{
	os<<" Auxiliary basis functions"<<std::endl;
	const std::vector<char> L_labels = {'s', 'p', 'd'};
	for(std::size_t T=0; T<orbs.size(); ++T)
	{
		os<<"\t\t"<<GlobalC::ucell.atoms[T].label<<"\t\t";
		for(std::size_t L=0; L<orbs[T].size(); ++L)
		{
			const char L_label =
				L < L_labels.size()
				? L_labels[L]
				: 'f' + (L-L_labels.size());
			os<<orbs[T][L].size()<<" "<<L_label<<"\t\t";
		}
		os<<std::endl;
	}
}
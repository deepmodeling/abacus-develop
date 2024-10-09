//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef MATRIX_ORB21_HPP
#define MATRIX_ORB21_HPP

#include "Matrix_Orbs21.h"
#include "RI_Util.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

template<typename Tdata>
RI::Tensor<Tdata> Matrix_Orbs21::cal_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
	const Matrix_Order &matrix_order) const
{
	RI::Tensor<Tdata> m;
	const size_t sizeA1 = index_A1[TA].count_size;
	const size_t sizeA2 = index_A2[TA].count_size;
	const size_t sizeB = index_B[TB].count_size;
	switch(matrix_order)
	{
		case Matrix_Order::A1A2B: m = RI::Tensor<Tdata>({sizeA1, sizeA2, sizeB});	break;
		case Matrix_Order::A1BA2: m = RI::Tensor<Tdata>({sizeA1, sizeB, sizeA2});	break;
		case Matrix_Order::BA1A2: m = RI::Tensor<Tdata>({sizeB, sizeA1, sizeA2});	break;
		case Matrix_Order::BA2A1: m = RI::Tensor<Tdata>({sizeB, sizeA2, sizeA1});	break;
		case Matrix_Order::A2A1B: m = RI::Tensor<Tdata>({sizeA2, sizeA1, sizeB});	break;
		case Matrix_Order::A2BA1: m = RI::Tensor<Tdata>({sizeA2, sizeB, sizeA1});	break;
		default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}

	for( const auto &co3 : center2_orb21_s.at(TA).at(TB) )
	{
		const int LA1 = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA1 = co4.first;
			for( size_t MA1=0; MA1!=2*LA1+1; ++MA1 )
			{
				for( const auto &co5 : co4.second )
				{
					const int LA2 = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NA2 = co6.first;
						for( size_t MA2=0; MA2!=2*LA2+1; ++MA2 )
						{
							for( const auto &co7 : co6.second )
							{
								const int LB = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB = co8.first;
									for( size_t MB=0; MB!=2*LB+1; ++MB )
									{
										const Tdata overlap = co8.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB );
										const size_t iA1 = index_A1[TA][LA1][NA1][MA1];
										const size_t iA2 = index_A2[TA][LA2][NA2][MA2];
										const size_t iB = index_B[TB][LB][NB][MB];
										switch(matrix_order)
										{
											case Matrix_Order::A1A2B:	m(iA1,iA2,iB) = overlap;	break;
											case Matrix_Order::A1BA2:	m(iA1,iB,iA2) = overlap;	break;
											case Matrix_Order::A2A1B:	m(iA2,iA1,iB) = overlap;	break;
											case Matrix_Order::A2BA1:	m(iA2,iB,iA1) = overlap;	break;
											case Matrix_Order::BA1A2:	m(iB,iA1,iA2) = overlap;	break;
											case Matrix_Order::BA2A1:	m(iB,iA2,iA1) = overlap;	break;
											default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return m;
}

template<typename Tdata>
std::array<RI::Tensor<Tdata>,3> Matrix_Orbs21::cal_grad_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
	const Matrix_Order &matrix_order) const
{
	std::array<RI::Tensor<Tdata>,3> m;
	const size_t sizeA1 = index_A1[TA].count_size;
	const size_t sizeA2 = index_A2[TA].count_size;
	const size_t sizeB = index_B[TB].count_size;
	for(int i=0; i<m.size(); ++i)
	{
		switch(matrix_order)
		{
			case Matrix_Order::A1A2B: m[i] = RI::Tensor<Tdata>({sizeA1, sizeA2, sizeB});	break;
			case Matrix_Order::A1BA2: m[i] = RI::Tensor<Tdata>({sizeA1, sizeB, sizeA2});	break;
			case Matrix_Order::BA1A2: m[i] = RI::Tensor<Tdata>({sizeB, sizeA1, sizeA2});	break;
			case Matrix_Order::BA2A1: m[i] = RI::Tensor<Tdata>({sizeB, sizeA2, sizeA1});	break;
			case Matrix_Order::A2A1B: m[i] = RI::Tensor<Tdata>({sizeA2, sizeA1, sizeB});	break;
			case Matrix_Order::A2BA1: m[i] = RI::Tensor<Tdata>({sizeA2, sizeB, sizeA1});	break;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	for( const auto &co3 : center2_orb21_s.at(TA).at(TB) )
	{
		const int LA1 = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA1 = co4.first;
			for( size_t MA1=0; MA1!=2*LA1+1; ++MA1 )
			{
				for( const auto &co5 : co4.second )
				{
					const int LA2 = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NA2 = co6.first;
						for( size_t MA2=0; MA2!=2*LA2+1; ++MA2 )
						{
							for( const auto &co7 : co6.second )
							{
								const int LB = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB = co8.first;
									for( size_t MB=0; MB!=2*LB+1; ++MB )
									{
										const std::array<double,3> grad_overlap = RI_Util::Vector3_to_array3(co8.second.cal_grad_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB ));
										const size_t iA1 = index_A1[TA][LA1][NA1][MA1];
										const size_t iA2 = index_A2[TA][LA2][NA2][MA2];
										const size_t iB = index_B[TB][LB][NB][MB];
										for(size_t i=0; i<m.size(); ++i)
										{
											switch(matrix_order)
											{
												case Matrix_Order::A1A2B:	m[i](iA1,iA2,iB) = grad_overlap[i];	break;
												case Matrix_Order::A1BA2:	m[i](iA1,iB,iA2) = grad_overlap[i];	break;
												case Matrix_Order::A2A1B:	m[i](iA2,iA1,iB) = grad_overlap[i];	break;
												case Matrix_Order::A2BA1:	m[i](iA2,iB,iA1) = grad_overlap[i];	break;
												case Matrix_Order::BA1A2:	m[i](iB,iA1,iA2) = grad_overlap[i];	break;
												case Matrix_Order::BA2A1:	m[i](iB,iA2,iA1) = grad_overlap[i];	break;
												default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return m;
}

template <typename Tdata> 
std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::vector<RI::Tensor<Tdata>>>>>> Matrix_Orbs21::cal_overlap_matrix_all(
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B) const
{
	ModuleBase::TITLE("Matrix_Orbs21","cal_overlap_matrix");

	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<RI::Tensor<Tdata>>>>>> matrixes;

	for( const auto &co1 : center2_orb21_s )
	{
		const size_t TA = co1.first;
		for( size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA )
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );

			for( const auto &co2 : co1.second )
			{
				const size_t TB = co2.first;
				for( size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					const ModuleBase::Vector3<double> &tauB( GlobalC::ucell.atoms[TB].tau[IB] );

					const RI::Tensor<Tdata> &&m = cal_overlap_matrix<Tdata>( TA, TB, tauA, tauB, index_A1, index_A2, index_B, Matrix_Order::A2BA1 );
					matrixes[TA][IA][TB][IB].resize(2);
					matrixes[TA][IA][TB][IB][0] = std::move(m);
					const RI::Tensor<Tdata> &&n = cal_overlap_matrix<Tdata>( TA, TB, tauA, tauB, index_A1, index_A2, index_B, Matrix_Order::BA2A1 );
					matrixes[TB][IB][TA][IA].resize(2);
					matrixes[TB][IB][TA][IA][1] = std::move(n);
                }
			}
		}
    }
	// matrixes[T][I][T][I][0] = matrixes[T][I][T][I][1], so delete repeat
    for (auto m1 : matrixes)
	{
		const size_t T = m1.first;
		for( auto m2 : m1.second )
		{
			const size_t I = m2.first;
			matrixes[T][I][T][I].resize(1);
		}
	}

	return matrixes;
}
#endif

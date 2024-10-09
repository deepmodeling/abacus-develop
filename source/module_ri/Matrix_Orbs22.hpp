//=======================
// AUTHOR : Peize Lin
// DATE :   2023-02-23
//=======================

#ifndef MATRIX_ORB22_HPP
#define MATRIX_ORB22_HPP

#include "Matrix_Orbs22.h"
#include "RI_Util.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

template<typename Tdata>
RI::Tensor<Tdata> Matrix_Orbs22::cal_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B2,
	const Matrix_Order &matrix_order) const
{
	RI::Tensor<Tdata> m;
	const size_t sizeA1 = index_A1[TA].count_size;
	const size_t sizeA2 = index_A2[TA].count_size;
	const size_t sizeB1 = index_B1[TB].count_size;
	const size_t sizeB2 = index_B2[TB].count_size;
	switch(matrix_order)
	{
		case Matrix_Order::A1A2B1B2: m = RI::Tensor<Tdata>({sizeA1, sizeA2, sizeB1, sizeB2});	break;
		case Matrix_Order::A1A2B2B1: m = RI::Tensor<Tdata>({sizeA1, sizeA2, sizeB2, sizeB1});	break;
		case Matrix_Order::A1B1A2B2: m = RI::Tensor<Tdata>({sizeA1, sizeB1, sizeA2, sizeB2});	break;
		case Matrix_Order::A1B1B2A2: m = RI::Tensor<Tdata>({sizeA1, sizeB1, sizeB2, sizeA2});	break;
		case Matrix_Order::A1B2A2B1: m = RI::Tensor<Tdata>({sizeA1, sizeB2, sizeA2, sizeB1});	break;
		case Matrix_Order::A1B2B1A2: m = RI::Tensor<Tdata>({sizeA1, sizeB2, sizeB1, sizeA2});	break;
		case Matrix_Order::A2A1B1B2: m = RI::Tensor<Tdata>({sizeA2, sizeA1, sizeB1, sizeB2});	break;
		case Matrix_Order::A2A1B2B1: m = RI::Tensor<Tdata>({sizeA2, sizeA1, sizeB2, sizeB1});	break;
		case Matrix_Order::A2B1A1B2: m = RI::Tensor<Tdata>({sizeA2, sizeB1, sizeA1, sizeB2});	break;
		case Matrix_Order::A2B1B2A1: m = RI::Tensor<Tdata>({sizeA2, sizeB1, sizeB2, sizeA1});	break;
		case Matrix_Order::A2B2A1B1: m = RI::Tensor<Tdata>({sizeA2, sizeB2, sizeA1, sizeB1});	break;
		case Matrix_Order::A2B2B1A1: m = RI::Tensor<Tdata>({sizeA2, sizeB2, sizeB1, sizeA1});	break;
		case Matrix_Order::B1A1A2B2: m = RI::Tensor<Tdata>({sizeB1, sizeA1, sizeA2, sizeB2});	break;
		case Matrix_Order::B1A1B2A2: m = RI::Tensor<Tdata>({sizeB1, sizeA1, sizeB2, sizeA2});	break;
		case Matrix_Order::B1A2A1B2: m = RI::Tensor<Tdata>({sizeB1, sizeA2, sizeA1, sizeB2});	break;
		case Matrix_Order::B1A2B2A1: m = RI::Tensor<Tdata>({sizeB1, sizeA2, sizeB2, sizeA1});	break;
		case Matrix_Order::B1B2A1A2: m = RI::Tensor<Tdata>({sizeB1, sizeB2, sizeA1, sizeA2});	break;
		case Matrix_Order::B1B2A2A1: m = RI::Tensor<Tdata>({sizeB1, sizeB2, sizeA2, sizeA1});	break;
		case Matrix_Order::B2A1A2B1: m = RI::Tensor<Tdata>({sizeB2, sizeA1, sizeA2, sizeB1});	break;
		case Matrix_Order::B2A1B1A2: m = RI::Tensor<Tdata>({sizeB2, sizeA1, sizeB1, sizeA2});	break;
		case Matrix_Order::B2A2A1B1: m = RI::Tensor<Tdata>({sizeB2, sizeA2, sizeA1, sizeB1});	break;
		case Matrix_Order::B2A2B1A1: m = RI::Tensor<Tdata>({sizeB2, sizeA2, sizeB1, sizeA1});	break;
		case Matrix_Order::B2B1A1A2: m = RI::Tensor<Tdata>({sizeB2, sizeB1, sizeA1, sizeA2});	break;
		case Matrix_Order::B2B1A2A1: m = RI::Tensor<Tdata>({sizeB2, sizeB1, sizeA2, sizeA1});	break;
		default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}

	for( const auto &co3 : center2_orb22_s.at(TA).at(TB) )
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
								const int LB1 = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB1 = co8.first;
									for( size_t MB1=0; MB1!=2*LB1+1; ++MB1 )
									{
										for( const auto &co9 : co8.second )
										{
											const int LB2 = co9.first;
											for( const auto &co10 : co9.second )
											{
												const size_t NB2 = co10.first;
												for( size_t MB2=0; MB2!=2*LB2+1; ++MB2 )
												{
													const Tdata overlap = co10.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB1, MB2 );
													const size_t iA1 = index_A1[TA][LA1][NA1][MA1];
													const size_t iA2 = index_A2[TA][LA2][NA2][MA2];
													const size_t iB1 = index_B1[TB][LB1][NB1][MB1];
													const size_t iB2 = index_B2[TB][LB2][NB2][MB2];
													switch(matrix_order)
													{
														case Matrix_Order::A1A2B1B2:	m(iA1,iA2,iB1,iB2) = overlap;	break;
														case Matrix_Order::A1A2B2B1:	m(iA1,iA2,iB2,iB1) = overlap;	break;
														case Matrix_Order::A1B1A2B2:	m(iA1,iB1,iA2,iB2) = overlap;	break;
														case Matrix_Order::A1B1B2A2:	m(iA1,iB1,iB2,iA2) = overlap;	break;
														case Matrix_Order::A1B2A2B1:	m(iA1,iB2,iA2,iB1) = overlap;	break;
														case Matrix_Order::A1B2B1A2:	m(iA1,iB2,iB1,iA2) = overlap;	break;
														case Matrix_Order::A2A1B1B2:	m(iA2,iA1,iB1,iB2) = overlap;	break;
														case Matrix_Order::A2A1B2B1:	m(iA2,iA1,iB2,iB1) = overlap;	break;
														case Matrix_Order::A2B1A1B2:	m(iA2,iB1,iA1,iB2) = overlap;	break;
														case Matrix_Order::A2B1B2A1:	m(iA2,iB1,iB2,iA1) = overlap;	break;
														case Matrix_Order::A2B2A1B1:	m(iA2,iB2,iA1,iB1) = overlap;	break;
														case Matrix_Order::A2B2B1A1:	m(iA2,iB2,iB1,iA1) = overlap;	break;
														case Matrix_Order::B1A1A2B2:	m(iB1,iA1,iA2,iB2) = overlap;	break;
														case Matrix_Order::B1A1B2A2:	m(iB1,iA1,iB2,iA2) = overlap;	break;
														case Matrix_Order::B1A2A1B2:	m(iB1,iA2,iA1,iB2) = overlap;	break;
														case Matrix_Order::B1A2B2A1:	m(iB1,iA2,iB2,iA1) = overlap;	break;
														case Matrix_Order::B1B2A1A2:	m(iB1,iB2,iA1,iA2) = overlap;	break;
														case Matrix_Order::B1B2A2A1:	m(iB1,iB2,iA2,iA1) = overlap;	break;
														case Matrix_Order::B2A1A2B1:	m(iB2,iA1,iA2,iB1) = overlap;	break;
														case Matrix_Order::B2A1B1A2:	m(iB2,iA1,iB1,iA2) = overlap;	break;
														case Matrix_Order::B2A2A1B1:	m(iB2,iA2,iA1,iB1) = overlap;	break;
														case Matrix_Order::B2A2B1A1:	m(iB2,iA2,iB1,iA1) = overlap;	break;
														case Matrix_Order::B2B1A1A2:	m(iB2,iB1,iA1,iA2) = overlap;	break;
														case Matrix_Order::B2B1A2A1:	m(iB2,iB1,iA2,iA1) = overlap;	break;														
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
		}
	}
	return m;
}

template<typename Tdata>
std::array<RI::Tensor<Tdata>,3> Matrix_Orbs22::cal_grad_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B2,
	const Matrix_Order &matrix_order) const
{
	std::array<RI::Tensor<Tdata>,3> m;
	const size_t sizeA1 = index_A1[TA].count_size;
	const size_t sizeA2 = index_A2[TA].count_size;
	const size_t sizeB1 = index_B1[TB].count_size;
	const size_t sizeB2 = index_B2[TB].count_size;
	for(int i=0; i<m.size(); ++i)
	{
		switch(matrix_order)
		{
			case Matrix_Order::A1A2B1B2: m[i] = RI::Tensor<Tdata>({sizeA1, sizeA2, sizeB1, sizeB2});	break;
			case Matrix_Order::A1A2B2B1: m[i] = RI::Tensor<Tdata>({sizeA1, sizeA2, sizeB2, sizeB1});	break;
			case Matrix_Order::A1B1A2B2: m[i] = RI::Tensor<Tdata>({sizeA1, sizeB1, sizeA2, sizeB2});	break;
			case Matrix_Order::A1B1B2A2: m[i] = RI::Tensor<Tdata>({sizeA1, sizeB1, sizeB2, sizeA2});	break;
			case Matrix_Order::A1B2A2B1: m[i] = RI::Tensor<Tdata>({sizeA1, sizeB2, sizeA2, sizeB1});	break;
			case Matrix_Order::A1B2B1A2: m[i] = RI::Tensor<Tdata>({sizeA1, sizeB2, sizeB1, sizeA2});	break;
			case Matrix_Order::A2A1B1B2: m[i] = RI::Tensor<Tdata>({sizeA2, sizeA1, sizeB1, sizeB2});	break;
			case Matrix_Order::A2A1B2B1: m[i] = RI::Tensor<Tdata>({sizeA2, sizeA1, sizeB2, sizeB1});	break;
			case Matrix_Order::A2B1A1B2: m[i] = RI::Tensor<Tdata>({sizeA2, sizeB1, sizeA1, sizeB2});	break;
			case Matrix_Order::A2B1B2A1: m[i] = RI::Tensor<Tdata>({sizeA2, sizeB1, sizeB2, sizeA1});	break;
			case Matrix_Order::A2B2A1B1: m[i] = RI::Tensor<Tdata>({sizeA2, sizeB2, sizeA1, sizeB1});	break;
			case Matrix_Order::A2B2B1A1: m[i] = RI::Tensor<Tdata>({sizeA2, sizeB2, sizeB1, sizeA1});	break;
			case Matrix_Order::B1A1A2B2: m[i] = RI::Tensor<Tdata>({sizeB1, sizeA1, sizeA2, sizeB2});	break;
			case Matrix_Order::B1A1B2A2: m[i] = RI::Tensor<Tdata>({sizeB1, sizeA1, sizeB2, sizeA2});	break;
			case Matrix_Order::B1A2A1B2: m[i] = RI::Tensor<Tdata>({sizeB1, sizeA2, sizeA1, sizeB2});	break;
			case Matrix_Order::B1A2B2A1: m[i] = RI::Tensor<Tdata>({sizeB1, sizeA2, sizeB2, sizeA1});	break;
			case Matrix_Order::B1B2A1A2: m[i] = RI::Tensor<Tdata>({sizeB1, sizeB2, sizeA1, sizeA2});	break;
			case Matrix_Order::B1B2A2A1: m[i] = RI::Tensor<Tdata>({sizeB1, sizeB2, sizeA2, sizeA1});	break;
			case Matrix_Order::B2A1A2B1: m[i] = RI::Tensor<Tdata>({sizeB2, sizeA1, sizeA2, sizeB1});	break;
			case Matrix_Order::B2A1B1A2: m[i] = RI::Tensor<Tdata>({sizeB2, sizeA1, sizeB1, sizeA2});	break;
			case Matrix_Order::B2A2A1B1: m[i] = RI::Tensor<Tdata>({sizeB2, sizeA2, sizeA1, sizeB1});	break;
			case Matrix_Order::B2A2B1A1: m[i] = RI::Tensor<Tdata>({sizeB2, sizeA2, sizeB1, sizeA1});	break;
			case Matrix_Order::B2B1A1A2: m[i] = RI::Tensor<Tdata>({sizeB2, sizeB1, sizeA1, sizeA2});	break;
			case Matrix_Order::B2B1A2A1: m[i] = RI::Tensor<Tdata>({sizeB2, sizeB1, sizeA2, sizeA1});	break;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	for( const auto &co3 : center2_orb22_s.at(TA).at(TB) )
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
								const int LB1 = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB1 = co8.first;
									for( size_t MB1=0; MB1!=2*LB1+1; ++MB1 )
									{
										for( const auto &co9 : co8.second )
										{
											const int LB2 = co9.first;
											for( const auto &co10 : co9.second )
											{
												const size_t NB2 = co10.first;
												for( size_t MB2=0; MB2!=2*LB2+1; ++MB2 )
												{
													const Tdata overlap = co10.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB1, MB2 );
													switch(matrix_order)
													{										
														const std::array<double,3> grad_overlap = RI_Util::Vector3_to_array3(co10.second.cal_grad_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB1, MB2 ));
														const size_t iA1 = index_A1[TA][LA1][NA1][MA1];
														const size_t iA2 = index_A2[TA][LA2][NA2][MA2];
														const size_t iB1 = index_B1[TB][LB1][NB1][MB1];
														const size_t iB2 = index_B2[TB][LB2][NB2][MB2];
														for(size_t i=0; i<m.size(); ++i)
														{
															switch(matrix_order)
															{
																case Matrix_Order::A1A2B1B2:	m[i](iA1,iA2,iB1,iB2) = grad_overlap;	break;
																case Matrix_Order::A1A2B2B1:	m[i](iA1,iA2,iB2,iB1) = grad_overlap;	break;
																case Matrix_Order::A1B1A2B2:	m[i](iA1,iB1,iA2,iB2) = grad_overlap;	break;
																case Matrix_Order::A1B1B2A2:	m[i](iA1,iB1,iB2,iA2) = grad_overlap;	break;
																case Matrix_Order::A1B2A2B1:	m[i](iA1,iB2,iA2,iB1) = grad_overlap;	break;
																case Matrix_Order::A1B2B1A2:	m[i](iA1,iB2,iB1,iA2) = grad_overlap;	break;
																case Matrix_Order::A2A1B1B2:	m[i](iA2,iA1,iB1,iB2) = grad_overlap;	break;
																case Matrix_Order::A2A1B2B1:	m[i](iA2,iA1,iB2,iB1) = grad_overlap;	break;
																case Matrix_Order::A2B1A1B2:	m[i](iA2,iB1,iA1,iB2) = grad_overlap;	break;
																case Matrix_Order::A2B1B2A1:	m[i](iA2,iB1,iB2,iA1) = grad_overlap;	break;
																case Matrix_Order::A2B2A1B1:	m[i](iA2,iB2,iA1,iB1) = grad_overlap;	break;
																case Matrix_Order::A2B2B1A1:	m[i](iA2,iB2,iB1,iA1) = grad_overlap;	break;
																case Matrix_Order::B1A1A2B2:	m[i](iB1,iA1,iA2,iB2) = grad_overlap;	break;
																case Matrix_Order::B1A1B2A2:	m[i](iB1,iA1,iB2,iA2) = grad_overlap;	break;
																case Matrix_Order::B1A2A1B2:	m[i](iB1,iA2,iA1,iB2) = grad_overlap;	break;
																case Matrix_Order::B1A2B2A1:	m[i](iB1,iA2,iB2,iA1) = grad_overlap;	break;
																case Matrix_Order::B1B2A1A2:	m[i](iB1,iB2,iA1,iA2) = grad_overlap;	break;
																case Matrix_Order::B1B2A2A1:	m[i](iB1,iB2,iA2,iA1) = grad_overlap;	break;
																case Matrix_Order::B2A1A2B1:	m[i](iB2,iA1,iA2,iB1) = grad_overlap;	break;
																case Matrix_Order::B2A1B1A2:	m[i](iB2,iA1,iB1,iA2) = grad_overlap;	break;
																case Matrix_Order::B2A2A1B1:	m[i](iB2,iA2,iA1,iB1) = grad_overlap;	break;
																case Matrix_Order::B2A2B1A1:	m[i](iB2,iA2,iB1,iA1) = grad_overlap;	break;
																case Matrix_Order::B2B1A1A2:	m[i](iB2,iB1,iA1,iA2) = grad_overlap;	break;
																case Matrix_Order::B2B1A2A1:	m[i](iB2,iB1,iA2,iA1) = grad_overlap;	break;
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
				}
			}
		}
	}
	return m;
}

// 4-parameter interface, for opt_orb
template<typename Tdata>
std::map < size_t, std::map<size_t, std::map<size_t, std::map<size_t, RI::Tensor<Tdata>>>>> Matrix_Orbs22::cal_overlap_matrix_all(
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B2 ) const
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t, RI::Tensor<Tdata>>>>> matrixes;

	for( const auto &co1 : center2_orb22_s )
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

					matrixes[TA][IA][TB][IB] = cal_overlap_matrix<Tdata>(
						TA,
						TB,
						GlobalC::ucell.atoms[TA].tau[IA],
						GlobalC::ucell.atoms[TB].tau[IB],
						index_A1,
						index_A2,
						index_B1,
						index_B2,
						Matrix_Order::A1B1A2B2);
				}
			}
		}
	}
	return matrixes;
}
#endif

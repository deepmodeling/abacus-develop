//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef MATRIX_ORB11_HPP
#define MATRIX_ORB11_HPP

#include "Matrix_Orbs11.h"
#include "RI_Util.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

template<typename Tdata>
RI::Tensor<Tdata> Matrix_Orbs11::cal_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
	const Matrix_Order &matrix_order) const
{
	RI::Tensor<Tdata> m;
	const size_t sizeA = index_A[TA].count_size;
	const size_t sizeB = index_B[TB].count_size;
	switch(matrix_order)
	{
		case Matrix_Order::AB: m = RI::Tensor<Tdata>({sizeA, sizeB});	break;
		case Matrix_Order::BA: m = RI::Tensor<Tdata>({sizeB, sizeA});	break;
		default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}

	for( const auto &co3 : center2_orb11_s.at(TA).at(TB) )
	{
		const int LA = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA = co4.first;
			for( size_t MA=0; MA!=2*LA+1; ++MA )
			{
				for( const auto &co5 : co4.second )
				{
					const int LB = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NB = co6.first;
						for( size_t MB=0; MB!=2*LB+1; ++MB )
						{
							const Tdata overlap = co6.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA, MB );
							const size_t iA = index_A[TA][LA][NA][MA];
							const size_t iB = index_B[TB][LB][NB][MB];
							switch(matrix_order)
							{
								case Matrix_Order::AB:	m(iA,iB) = overlap;		break;
								case Matrix_Order::BA:	m(iB,iA) = overlap;		break;
								default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
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
std::array<RI::Tensor<Tdata>,3> Matrix_Orbs11::cal_grad_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
	const Matrix_Order &matrix_order) const
{
	std::array<RI::Tensor<Tdata>,3> m;
	const size_t sizeA = index_A[TA].count_size;
	const size_t sizeB = index_B[TB].count_size;
	for(int i=0; i<m.size(); ++i)
	{
		switch(matrix_order)
		{
			case Matrix_Order::AB: m[i] = RI::Tensor<Tdata>({sizeA, sizeB});	break;
			case Matrix_Order::BA: m[i] = RI::Tensor<Tdata>({sizeB, sizeA});	break;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	for( const auto &co3 : center2_orb11_s.at(TA).at(TB) )
	{
		const int LA = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA = co4.first;
			for( size_t MA=0; MA!=2*LA+1; ++MA )
			{
				for( const auto &co5 : co4.second )
				{
					const int LB = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NB = co6.first;
						for( size_t MB=0; MB!=2*LB+1; ++MB )
						{
							const std::array<double,3> grad_overlap = RI_Util::Vector3_to_array3(co6.second.cal_grad_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA, MB ));
							const size_t iA = index_A[TA][LA][NA][MA];
							const size_t iB = index_B[TB][LB][NB][MB];
							for(size_t i=0; i<m.size(); ++i)
							{
								switch(matrix_order)
								{
									case Matrix_Order::AB:	m[i](iA,iB) = grad_overlap[i];		break;
									case Matrix_Order::BA:	m[i](iB,iA) = grad_overlap[i];		break;
									default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
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
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<Tdata>>>>> Matrix_Orbs11::cal_overlap_matrix_all( 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_r, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_c ) const
{
	ModuleBase::TITLE("Matrix_Orbs11","cal_overlap_matrix");
	
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<Tdata>>>>> matrixes;
	
	for( const auto &co1 : center2_orb11_s )
	{
		const size_t TA = co1.first;
		for (size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA)
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );

			for( const auto &co2 : co1.second )
			{
				const size_t TB = co2.first;
				for (size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB)
				{
					const ModuleBase::Vector3<double> &tauB( GlobalC::ucell.atoms[TB].tau[IB] );

					matrixes[TA][IA][TB][IB] = cal_overlap_matrix<Tdata>( TA, TB, tauA, tauB, index_r, index_c, Matrix_Order::AB );
				}
			}
		}
	}
	return matrixes;
}

#endif

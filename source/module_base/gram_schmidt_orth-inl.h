//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-08-03
//==========================================================

#ifndef GRAM_SCHMIDT_ORTH_INL_H
#define GRAM_SCHMIDT_ORTH_INL_H

#include "gram_schmidt_orth.h"

#include "mathzone.h"
#include "blas_connector.h"
#include "math_integral.h" // mohan add 2021-04-03
namespace ModuleBase
{

template<typename Func_Type, typename R_Type>
Gram_Schmidt_Orth<Func_Type,R_Type>::Gram_Schmidt_Orth( const std::vector<R_Type> &rab_in, const Coordinate &coordinate_in )
	:rab(rab_in),
	 coordinate(coordinate_in)
{	
	if( Coordinate::Sphere == coordinate )
	{
		std::vector<R_Type> radial( rab.size() );
		radial[0] = 0;
		for( int ir=1; ir!=radial.size(); ++ir )
			radial[ir] = radial[ir-1] + rab[ir-1];						
		this->radial_2 = Mathzone::Pointwise_Product( radial, radial );
	}
}

template<typename Func_Type, typename R_Type>
std::vector<std::vector<Func_Type>> Gram_Schmidt_Orth<Func_Type,R_Type>::cal_orth( 
	const std::vector<std::vector<Func_Type>> &func,
	const Func_Type norm_threshold )
{
	//	Schmidt: hn to en
	//		e1 = h1 / ||h1||
	//		gn = hn - \sum{i=1 to n-1}(hn,ei)ei
	//		en = gn / ||gn||
	
	std::vector<std::vector<Func_Type>> func_new;
	
	for( size_t if1=0; if1!=func.size(); ++if1 )
	{
		//use CGS2 algorithm to do twice orthogonalization 
		//DOI 10.1007/s00211-005-0615-4
		std::vector<Func_Type> func_try = func[if1];
		for(int niter=0;niter<3;niter++)
		{
			std::vector<Func_Type> func_tmp = func_try;
			for( size_t if_minus=0; if_minus!=func_new.size(); ++if_minus )
			{
				// (hn,ei)
				const std::vector<Func_Type> && mul_func = Mathzone::Pointwise_Product( func_tmp, func_new[if_minus] );
				const Func_Type in_product = cal_norm(mul_func);

				// hn - (hn,ei)ei
				BlasConnector::axpy( mul_func.size(), -in_product, ModuleBase::GlobalFunc::VECTOR_TO_PTR(func_new[if_minus]), 1, ModuleBase::GlobalFunc::VECTOR_TO_PTR(func_try), 1);
			}
		}
		
		// ||gn||
		const std::vector<Func_Type> && func_2 = Mathzone::Pointwise_Product( func_try, func_try );
		const Func_Type norm = sqrt(cal_norm(func_2));
				
		// en = gn / ||gn||
		// if ||gn|| too small, filter out
		if( norm >= norm_threshold )
		{
			BlasConnector::scal( func_try.size(), 1.0/norm, ModuleBase::GlobalFunc::VECTOR_TO_PTR(func_try), 1 );
			func_new.push_back( func_try );
		}
	}
	return func_new;
}

// cal ||f||
template<typename Func_Type, typename R_Type>
Func_Type Gram_Schmidt_Orth<Func_Type,R_Type>::cal_norm( const std::vector<Func_Type> &f )
{
	Func_Type norm = 0.0;
	switch( this->coordinate )
	{
		case Coordinate::Cartesian:
		{
			Integral::Simpson_Integral( f.size(), ModuleBase::GlobalFunc::VECTOR_TO_PTR(f), ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab), norm);		
			break;
		}
		case Coordinate::Sphere:	
		{
			const std::vector<Func_Type> &&tmp_func = Mathzone::Pointwise_Product( f, radial_2 );
			Integral::Simpson_Integral( f.size(), ModuleBase::GlobalFunc::VECTOR_TO_PTR(tmp_func), ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab), norm);	
			break;
		}
		default:
		{
			throw std::invalid_argument("coordinate must be Cartesian or Sphere "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			break;
		}
	}
	return norm;
}

}

#endif	// GRAM_SCHMIDT_ORTH_INL_H

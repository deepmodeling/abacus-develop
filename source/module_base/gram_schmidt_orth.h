//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-08-03
//==========================================================

#ifndef GRAM_SCHMIDT_ORTH_H
#define GRAM_SCHMIDT_ORTH_H

#include<limits>
#include<vector>
namespace ModuleBase
{

template<typename Func_Type, typename R_Type=double>
class Gram_Schmidt_Orth
{
public:

	enum class Coordinate { Cartesian, Sphere };
	
	Gram_Schmidt_Orth( const std::vector<R_Type> &rab, const Coordinate &coordinate );
	
	std::vector<std::vector<Func_Type>> cal_orth( 
		const std::vector<std::vector<Func_Type>> &func,
		const Func_Type norm_threshold = std::numeric_limits<Func_Type>::min() );
	
private:

	const Coordinate coordinate;
	const std::vector<R_Type> &rab;
	std::vector<R_Type> radial_2;	

	Func_Type cal_norm( const std::vector<Func_Type> &f );

};

}
#endif	// GRAM_SCHMIDT_ORTH_H
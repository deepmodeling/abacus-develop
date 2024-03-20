#ifndef RANDOM_H
#define RANDOM_H
#include <cstdlib>
#include <cassert>

namespace ModuleBase
{

class Random
{
	public:
	Random();
	~Random();

	static void between0and1( double *v, const int &num )
	{
		assert( v!= NULL);
		assert( num > 1);
		for(int i=0; i<num; i++)
		{
			v[i] = static_cast<double>( std::rand() ) / RAND_MAX; 
		}
	}

	static double betweenMinus2and2(void)
	{
		return 2.0*betweenMinus1and1();
	}

	static double betweenMinus1and1(void)
	{
		const int a = std::rand() % 2;
		if(a==0) return between0and1();
		else if(a==1) return betweenMinus1and0();
		else throw(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	// Peize Lin add to fix warning 2019-05-01
	}

	static double between0and1(void)
	{
		return static_cast<double>( std::rand() )/RAND_MAX;
	}

	static double betweenMinus1and0(void)
	{
		return -static_cast<double>( std::rand() )/RAND_MAX;
	}

};

}

#endif

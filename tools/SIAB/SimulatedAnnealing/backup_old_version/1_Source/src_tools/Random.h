#ifndef RANDOM_H
#define RANDOM_H
#include <cstdlib>
#include <cassert>

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
		int a = std::rand() % 2;
		if(a==0) return between0and1();
		else if(a==1) return betweenMinus1and0();
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

#endif

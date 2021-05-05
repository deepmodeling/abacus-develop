//==========================================================
// AUTHOR : Peize Lin
// DATE : 2017-04-24
//==========================================================

#ifndef SPH_BESSEL_RECURSIVE_H
#define SPH_BESSEL_RECURSIVE_H

#include<vector>
using namespace std;

class Sph_Bessel_Recursive
{
public:
	class D1;
	class D2;
};



class Sph_Bessel_Recursive::D1
{
public:	
	const vector<vector<double>> & cal_jlx( const int lmax, const size_t ix_size );
	const vector<vector<double>> & get_jlx()const{ return jlx; }
	
	void set_dx(const double dx_in);
	double get_dx()const{ return dx; }

private:
	vector<vector<double>> jlx;		// jlx[l][x]
	double dx;
	bool finish_set_dx = false;
	
	void cal_jlx_0( const int l_size );
	void cal_jlx_smallx( const int l_size, const size_t ix_size );
	void cal_jlx_recursive( const int l_size, const size_t ix_size );
	
	
	double threshold = 1e-8;					// Peize Lin test
};



class Sph_Bessel_Recursive::D2
{
public:	
	const vector<vector<vector<double>>> & cal_jlx( const int lmax, const size_t ix1_size, const size_t ix2_size );
	const vector<vector<vector<double>>> & get_jlx()const{ return jlx; }
	
	void set_dx(const double dx_in);
	double get_dx()const{ return dx; }

private:
	vector<vector<vector<double>>> jlx;		// jlx[l][x1][x2]
	double dx;
	bool finish_set_dx = false;
	
	void cal_jlx_0( const int l_size, const size_t ix1_size, const size_t ix2_size );
	void cal_jlx_smallx( const int l_size, const size_t ix1_size, const size_t ix2_size );
	void cal_jlx_recursive( const int l_size, const size_t ix1_size, const size_t ix2_size );
	
	
	double threshold = 1e-8;					// Peize Lin test
};



class Sph_Bessel_Recursive_Pool
{
public:
	class D1
	{
		public:
		static vector<Sph_Bessel_Recursive::D1> sb_pool;
	};
	class D2
	{
		public:
		static vector<Sph_Bessel_Recursive::D2> sb_pool;
	};
};

#endif	// SPH_BESSEL_RECURSIVE_H
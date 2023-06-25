#ifndef GLOBAL_FUNCTION_H
#define GLOBAL_FUNCTION_H

#include <vector>
#include <valarray>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <complex>
#include <cassert>

#include "tool_title.h" // mohan add 2021-05-05
#include "tool_quit.h" // mohan add 2021-05-07
#include "tool_check.h" // mohan add 2021-05-08
#include "global_variable.h"
#include "global_function-func_each_2.h"		// Peize Lin add 2016-09-07

#include "blas_connector.h"

namespace ModuleBase
{
namespace GlobalFunc
{

void NOTE(const std::string &words);
void NEW_PART(const std::string &words);

//==========================================================
// GLOBAL FUNCTION :
// NAME : OUT( output date for checking )
// NAME : OUT( output date into file "ofs")
// NAME : OUTP( output parameters )
//==========================================================
void OUT(std::ofstream &ofs,const std::string &name);

template <class T>
void OUT(std::ofstream &ofs,const std::string &name,const T &a)
{
	std::stringstream name2;
	name2 << name ;
    ofs<< " " << std::setw(40) << name2.str() << " = " << a <<std::endl;
//	ofs << " " << name << a << std::endl;
    return;
}

template <class T>
void OUT(std::ofstream &ofs,const std::string &name,const T &x, const T&y)
{
    ofs<< " " << std::setw(40) <<name << " = [ " << x << ", " << y <<" ]" << std::endl;
//	ofs << " " << name << a << std::endl;
    return;
}

template <class T>
void OUT(std::ofstream &ofs,const std::string &name,const T &x, const T &y, const T &z)
{
    ofs<< " " << std::setw(40) <<name << " = [ " << x << ", " << y << ", " << z <<" ]" << std::endl;
    return;
}




// output parameters and explanations
template <class T>
void OUTP(std::ofstream &ofs, const std::string &name, const T &a, const std::string &explanation="")
{
	ofs << std::setw(30) << name << " " << a << " #" << explanation << std::endl;
}

template <class T>
void OUT(const std::string &name,const T &a)
{
    std::cout << " " << std::setw(40) << name << " = " << a << std::endl;
//	std::cout << " " << name << a << std::endl;
    return;
}


void OUT_TIME(const std::string &name, time_t &start, time_t &end);

//==========================================================
// GLOBAL FUNCTION :
// NAME : MAKE_DIR( make dir ,using system function)
//==========================================================
void MAKE_DIR( const std::string &file );

//==========================================================
// GLOBAL FUNCTION :
// NAME : AUTO_SET( auto_set variables )
//==========================================================
template <class T>
void AUTO_SET(const std::string &name,const T &a)
{
    GlobalV::ofs_warning <<" AUTO_SET "<<name<<" to "<<a << std::endl;
    return;
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : DONE( ouput information(time) on screen and log)
// 		  we can regard it as a milestone.
//==========================================================
void DONE(std::ofstream &ofs,const std::string &description, bool only_rank0 = false);

//==========================================================
// GLOBAL FUNCTION :
// NAME : ZEROS
// set elements of u as zero which u is 1_d std::complex array
//==========================================================
template<class T, class TI>
inline void ZEROS(std::complex<T> *u,const TI n)		// Peize Lin change int to TI at 2020.03.03
{
    assert(n>=0);
    for (TI i=0;i<n;i++)
    {
        u[i] = std::complex<T>(0.0,0.0);
    }
    return;
}

template<class T, class TI>
inline void ZEROS(T *u,const TI n)		// Peize Lin change int to TI at 2020.03.03
{
    assert(n>=0);
    for (TI i=0;i<n;i++)
    {
        u[i] = 0;
    }
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : TEST_LEVEL
// control the test_level
//==========================================================
void TEST_LEVEL(const std::string &name, bool disable);


//==========================================================
// GLOBAL FUNCTION :
//==========================================================
template <class T>
static void READ_VALUE(std::ifstream &ifs, T &v)
{
    ifs >> v;
    ifs.ignore(150, '\n');
    return;
}

bool SCAN_BEGIN(std::ifstream &ifs, const std::string &TargetName, const bool restart=1, const bool ifwarn=true);
// ifwarn: whether to call GlobalV::ofs_warning when the TargetName is not found, used to avoid invalid warning.
// Mohan warning : the last term can't be written as const bool &restart,
// I don't know why.

void SCAN_END(std::ifstream &ifs, const std::string &TargetName, const bool ifwarn=true);
// ifwarn: whether to call GlobalV::ofs_warning when the TargetName is not found, used to avoid invalid warning.

template<class T>
static inline void DCOPY( const T &a, T &b, const int &dim)
{
    for (int i=0; i<dim; ++i) b[i] = a[i];
}

template<typename T>
inline void COPYARRAY(const T* a, T* b, const long dim);

template<>
inline void COPYARRAY(const std::complex<double>* a, std::complex<double>* b, const long dim)
{
    const int one = 1;
    zcopy_(&dim, a, &one, b, &one);
}

template<>
inline void COPYARRAY(const double* a, double* b, const long dim)
{
    const int one = 1;
    dcopy_(&dim, a, &one, b, &one);
}

void BLOCK_HERE( const std::string &description );

//==========================================================
// GLOBAL FUNCTION :
// NAME : VECTOR_TO_PTR
// change std::vector to pointer
// Peize Lin add 2016-02-25
//==========================================================
template<class T>
static inline T * VECTOR_TO_PTR( std::vector<T> & v )
{
    return &(v[0]);
}
template<class T>
static inline T * VECTOR_TO_PTR( std::valarray<T> & v )
{
    return &(v[0]);
}

template<class T>
static inline const T * VECTOR_TO_PTR( const std::vector<T> & v )
{
    return &(v[0]);
}
template<class T>
static inline const T * VECTOR_TO_PTR( const std::valarray<T> & v )
{
    return &(v[0]);
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : TO_STRING
// change number to std::string
// example: 233 -> "233"
// Peize Lin add 2016-07-18
//==========================================================
template< typename T >
std::string TO_STRING ( const T &n )
{
	std::stringstream newstr;
	newstr<<n;
	return newstr.str();
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : MAP_EXIST
// find whether exists the std::map index
// if exist return the ptr called, else return nullptr
// example: map_exist(ms,i,j,k) -> try to find ms[i][j][k]
// Peize Lin add 2018-07-16
//==========================================================
template< typename T_map, typename T_key1  >
inline void* MAP_EXIST( T_map &ms, const T_key1 &key1 )
{
	auto ms1 = ms.find(key1);
	if( ms1 == ms.end() )	return nullptr;
	return static_cast<void*>(&ms1->second);
}

template< typename T_map, typename T_key1, typename... T_key_tail >
inline void* MAP_EXIST( T_map &ms, const T_key1 &key1, const T_key_tail&... key_tail )
{
	auto ms1 = ms.find(key1);
	if( ms1 == ms.end() )	return nullptr;
	return MAP_EXIST( ms1->second, key_tail... );
}

template< typename T_map, typename T_key1  >
inline const void* MAP_EXIST( const T_map &ms, const T_key1 &key1 )
{
	auto ms1 = ms.find(key1);
	if( ms1 == ms.end() )	return nullptr;
	return static_cast<const void*>(&ms1->second);
}

template< typename T_map, typename T_key1, typename... T_key_tail >
inline const void* MAP_EXIST( const T_map &ms, const T_key1 &key1, const T_key_tail&... key_tail )
{
	auto ms1 = ms.find(key1);
	if( ms1 == ms.end() )	return nullptr;
	return MAP_EXIST( ms1->second, key_tail... );
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : MemAvailable
// read /proc/meminfo
// unit: kB
// Peize Lin add 2019-12-21
//==========================================================
size_t MemAvailable();


//==========================================================
// GLOBAL FUNCTION :
// NAME : DELETE_MUL_PTR
// delete Multi-dimensional array pointer
// example:
//		int*** v;
//		DELETE_MUL_PTR(v,N1,N2);
//	->	for(int i1=0; i1<N1; ++i1){
//			for(int i2=0; i2<N2; ++i2){
//				delete[] v[i1][i2];	v[i1][i2]=nullptr;	}
//			delete[] v[i1];	v[i1]=nullptr;	}
//		delete[] v;	v=nullptr;
// Peize Lin add 2021-05-09
//==========================================================
template <typename T_element>
static inline void DELETE_MUL_PTR(T_element* v)
{
	delete[] v;
	v = nullptr;
}
template <typename T_element, typename T_N_first, typename... T_N_tail>
static inline void DELETE_MUL_PTR(T_element* v, const T_N_first N_first, const T_N_tail... N_tail)
{
	for(T_N_first i=0; i<N_first; ++i)
		DELETE_MUL_PTR(v[i],N_tail...);
	delete[] v;
	v = nullptr;
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : FREE_MUL_PTR
// delete Multi-dimensional array pointer
// example:
//		int*** v;
//		DELETE_MUL_PTR(v,N1,N2);
//	->	for(int i1=0; i1<N1; ++i1){
//			for(int i2=0; i2<N2; ++i2){
//				free(v[i1][i2]);	v[i1][i2]=nullptr;	}
//			free(v[i1]);	v[i1]=nullptr;	}
//		free(v);	v=nullptr;
// Peize Lin add 2021-05-09
//==========================================================
template <typename T_element>
static inline void FREE_MUL_PTR(T_element* v)
{
	free(v);
	v = nullptr;
}
template <typename T_element, typename T_N_first, typename... T_N_tail>
static inline void FREE_MUL_PTR(T_element* v, const T_N_first N_first, const T_N_tail... N_tail)
{
	for(T_N_first i=0; i<N_first; ++i)
		FREE_MUL_PTR(v[i],N_tail...);
	free(v);
	v = nullptr;
}

double ddot_real(
        const int & dim,
        const std::complex<double>* psi_L,
        const std::complex<double>* psi_R,
        const bool reduce = true) ;

//==========================================================
// GLOBAL FUNCTION :
// NAME : IS_COLUMN_MAJOR_KS_SOLVER
// check ks_solver requires column major or not
//==========================================================
static inline bool IS_COLUMN_MAJOR_KS_SOLVER()
{
    return GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx" || GlobalV::KS_SOLVER=="cusolver";
}

}//namespace GlobalFunc
}//namespace ModuleBase

#endif

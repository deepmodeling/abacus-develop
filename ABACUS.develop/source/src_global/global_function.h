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

using namespace std;

void NOTE(const string &words);
void NEW_PART(const string &words);

//==========================================================
// GLOBAL FUNCTION :
// NAME : OUT( output date for checking )
// NAME : OUT( output date into file "ofs")
// NAME : OUTP( output parameters )
//==========================================================
void OUT(ofstream &ofs,const string &name);

template <class T>
void OUT(ofstream &ofs,const string &name,const T &a)
{
	stringstream name2;
	name2 << name ;
    ofs<< " " << setw(40) << name2.str() << " = " << a <<endl;
//	ofs << " " << name << a << endl;
    return;
}

template <class T>
void OUT(ofstream &ofs,const string &name,const T &x, const T&y)
{
	stringstream name2;
	name2 << "[" << name << "]";
    ofs<< " " << setw(40) <<name2.str() << " = " << x << ", " << y << endl;
//	ofs << " " << name << a << endl;
    return;
}

template <class T>
void OUT(ofstream &ofs,const string &name,const T &x, const T &y, const T &z)
{
	stringstream name2;
	name2 << "[" << name << "]";
    ofs<< " " << setw(40) <<name2.str() << " = " << x << ", " << y << ", " << z << endl;
    return;
}




// output parameters and explanations
template <class T>
void OUTP(ofstream &ofs, const string &name, const T &a, const string &explanation="")
{
	ofs << setw(20) << name << a << " #" << explanation << endl;
}

template <class T>
void OUT(const string &name,const T &a)
{
    cout << " " << setw(40) << name << " = " << a;
//	cout << " " << name << a << endl;
    return;
}


void OUT_TIME(const string &name, time_t &start, time_t &end);

//==========================================================
// GLOBAL FUNCTION :
// NAME : MAKE_DIR( make dir ,using system function)
//==========================================================
void MAKE_DIR( const string &file );

//==========================================================
// GLOBAL FUNCTION :
// NAME : AUTO_SET( auto_set variables )
//==========================================================
template <class T>
void AUTO_SET(const string &name,const T &a)
{
    ofs_warning <<" AUTO_SET "<<name<<" to "<<a << endl;
    return;
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : DONE( ouput information(time) on screen and log)
// 		  we can regard it as a milestone.
//==========================================================
void DONE(ofstream &ofs,const string &description, bool only_rank0 = false);

//==========================================================
// GLOBAL FUNCTION :
// NAME : ZEROS
// set elements of u as zero which u is 1_d complex array
//==========================================================
template<class T, class TI>
inline void ZEROS(complex<T> *u,const TI n)		// Peize Lin change int to TI at 2020.03.03
{
    assert(n>=0);
    assert(u!=0);
    for (TI i=0;i<n;i++)
    {
        u[i] = complex<T>(0.0,0.0);
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
void TEST_LEVEL(const string &name);


//==========================================================
// GLOBAL FUNCTION :
//==========================================================
template <class T>
static void READ_VALUE(ifstream &ifs, T &v)
{
    ifs >> v;
    ifs.ignore(150, '\n');
    return;
}

bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart=1);
// Mohan warning : the last term can't be written as const bool &restart,
// I don't know why.

void SCAN_END(ifstream &ifs, const string &TargetName);

template<class T>
static inline void DCOPY( const T &a, T &b, const int &dim)
{
    for (int i=0; i<dim; ++i) b[i] = a[i];
}

void BLOCK_HERE( const string &description );

//==========================================================
// GLOBAL FUNCTION :
// NAME : VECTOR_TO_PTR
// change vector to pointer
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
// change number to string
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
// find whether exists the map index
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

#endif

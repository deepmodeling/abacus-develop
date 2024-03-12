#ifndef CONTAINER_OPERATOR_H
#define CONTAINER_OPERATOR_H

#include <cassert>
#include <vector>
#include <map>

template< typename T>
std::vector<T> operator + ( const std::vector<T> & x1, const std::vector<T> & x2 )
{
	assert(x1.size()==x2.size());
	std::vector<T> x;
	for(std::size_t i=0; i!=x1.size(); ++i )
		x.push_back(x1[i]+x2[i]);
	return x;
}

template< typename T>
std::vector<T> operator - ( const std::vector<T> & x1, const std::vector<T> & x2 )
{
	assert(x1.size()==x2.size());
	std::vector<T> x;
	for(std::size_t i=0; i!=x1.size(); ++i )
		x.push_back(x1[i]-x2[i]);
	return x;
}

template< typename T1, typename T2 >
std::map<T1,T2> operator + ( const std::map<T1,T2> & x1, const std::map<T1,T2> & x2 )
{
	assert(x1.size()==x2.size());
	std::map<T1,T2> x;
	for( const auto &x1i : x1 )
		x.insert(std::make_pair( x1i.first, x1i.second + x2.at(x1i.first) ));
	return x;
}

template< typename T1, typename T2 >
std::map<T1,T2> operator - ( const std::map<T1,T2> & x1, const std::map<T1,T2> & x2 )
{
	assert(x1.size()==x2.size());
	std::map<T1,T2> x;
	for( const auto &x1i : x1 )
		x.insert(std::make_pair( x1i.first, x1i.second - x2.at(x1i.first) ));
	return x;
}

template< typename T1, typename T2 >
std::vector<T2> operator * ( const T1 & x1, const std::vector<T2> & x2 )
{
	std::vector<T2> x;
	for(std::size_t i=0; i!=x2.size(); ++i )
		x.push_back(x1*x2[i]);
	return x;
}

template< typename T1, typename T21, typename T22 >
std::map<T21,T22> operator * ( const T1 & x1, const std::map<T21,T22> & x2 )
{
	std::map<T21,T22> x;
	for( const auto & x2i : x2 )
		x.insert(std::make_pair( x2i.first, x1*x2i.second ));
	return x;
}

#endif	// CONTAINER_OPERATOR_H

// AUTHOR:	Peize Lin
// Date: 	2016-09-07

#ifndef FUNC_EACH_2_H
#define FUNC_EACH_2_H

#include <functional>
#include <map>
#include <vector>

namespace ModuleBase
{
namespace GlobalFunc
{
	
template<typename Ti, typename... T_tail>
void FUNC_EACH_2( 
	Ti & tA, 
	const Ti & tB, 
	std::function< void( Ti&, const Ti&, T_tail... ) > func, 
	const T_tail&... t_tail )
{
	func( tA, tB, t_tail... );
}


template<typename Tv, typename Ti, typename... T_tail>
void FUNC_EACH_2( 
	std::vector<Tv> & tA, 
	const std::vector<Tv> & tB, 
	std::function< void( Ti&, const Ti&, T_tail... ) > func,
	const T_tail&... t_tail )
{
	for( size_t i=0; i!=tA.size(); ++i )
	{
		FUNC_EACH_2( tA[i], tB[i], func, t_tail... );
	}
}


template<typename T1, typename T2, typename Ti, typename... T_tail>
void FUNC_EACH_2( 
	std::map<T1,T2> & tA, 
	const std::map<T1,T2> & tB, 
	std::function< void( Ti&, const Ti&, T_tail... ) > func,
	const T_tail&... t_tail ) 
{
	for( auto & ta : tA )
	{
		FUNC_EACH_2( ta.second, tB.at(ta.first), func, t_tail... );
	}
}

}
}
#endif // FUNC_EACH_2_H

//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#ifndef ELEMENT_BASIS_INDEX_H
#define ELEMENT_BASIS_INDEX_H

#include <cstddef>
#include <vector>

namespace ModuleBase
{

namespace Element_Basis_Index
{
  //private:

	struct NM
	{
		public:
		std::size_t N;
		std::size_t M;
	};

	class Index_TL: public std::vector<std::vector<std::size_t>>
	{
		public:
		std::size_t N;
		std::size_t M;
	};

	class Index_T: public std::vector<Index_TL>
	{
		public:
		std::size_t count_size;
	};

  //public:

	typedef std::vector<std::vector<NM>> Range; 						// range[T][L]
	typedef std::vector<Index_T> IndexLNM;								// index[T][L][N][M]

	extern IndexLNM construct_index( const Range &range );
}

}

#endif
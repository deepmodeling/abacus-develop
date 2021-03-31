#ifndef EXX_ABFS_SCREEN_SCHWARZ
#define EXX_ABFS_SCREEN_SCHWARZ

#include "exx_abfs.h"
#include "abfs.h"
#include "src_global/matrix.h"
#include <map>
#include <atomic>

// |\braket{ik|jl}| \leq \sqrt{ \braket{ik|ik} \braket{jl|jl} }

class Exx_Abfs::Screen::Schwarz
{
public:
	void init( 
		const bool flag_screen_schwarz_in,
		const double threshold_in);
	void cal_max_pair_fock(
		const set<size_t> &atom_centres,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
		const Element_Basis_Index::IndexLNM &index_abfs,
		const Element_Basis_Index::IndexLNM &index_lcaos,
		const Abfs::Vector3_Order<int> &Born_von_Karman_period,
		map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> &Cws,
		map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> &Vws );
	inline bool screen(
		const size_t iat1, const size_t iat2, const size_t iat3, const size_t iat4,
		const Abfs::Vector3_Order<int> & box3, const Abfs::Vector3_Order<int> & box4 ) const;
		
private:
	bool flag_screen_schwarz = false;
	double threshold = 0;
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> max_pair_fock;
	
public:
//	static double num_screen;
//	static double num_cal;
};


bool Exx_Abfs::Screen::Schwarz::screen(
	const size_t iat1, const size_t iat2, const size_t iat3, const size_t iat4,
	const Abfs::Vector3_Order<int> & box3, const Abfs::Vector3_Order<int> & box4 ) const
{
	if(!flag_screen_schwarz)
		return false;
	if( max_pair_fock.at(iat1).at(iat3).at(box3) * max_pair_fock.at(iat2).at(iat4).at(box4) > threshold )
	{
//		num_cal += 1;
		return false;
	}
	else
	{
//		num_screen += 1;
		return true;
	}
}

#endif
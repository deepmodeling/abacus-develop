#ifndef EXX_ABFS_SCREEN_CAUCHY_H
#define EXX_ABFS_SCREEN_CAUCHY_H

#include "exx_abfs.h"
#include "src_ri/abfs-vector3_order.h"
#include <map>
#include <valarray>

class Exx_Abfs::Screen::Cauchy
{
public:
	struct Info_Step
	{
		double C1_norm4_max = 0.0;
		double C3_norm4_max = 0.0;
		double C2_norm4_max = 0.0;
		double C4_norm4_max = 0.0;
		
		double C1_norm2_max = 0.0;
		double C3_norm2_max = 0.0;
		
		double V_norm4 = 0.0;
		
		double D34_norm4_max = 0.0;
		double D32_norm4_max = 0.0;
		double D14_norm4_max = 0.0;
		double D12_norm4_max = 0.0;
	};
	
public:

	void init(
		const bool flag_screen_cauchy_in,
		const double threshold_in,
		const Abfs::Vector3_Order<int> Born_von_Karman_period_in);
	void cal_norm_C_max( 
		const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> & Cs,
		const Element_Basis_Index::IndexLNM & index_lcaos,
		const Element_Basis_Index::IndexLNM & index_abfs);
	void cal_norm_V( const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> & Vs );
	void cal_norm_D_max( const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> & Ds );
		
	Info_Step input_info(
		const size_t iat1, const size_t iat2, const size_t iat3, const size_t iat4,
		const Abfs::Vector3_Order<int> & box2,
		const Abfs::Vector3_Order<int> & box3,
		const Abfs::Vector3_Order<int> & box4 ) const;
	int postcalA( const Info_Step & info_step ) const;
	int postcalB(
		const Info_Step & info_step,
		const matrix & VC_T,			// iw2, \mu1, iw4
		const size_t range_iw2,
		const size_t range_mu1,
		const size_t range_iw4,
		const int postcal_in) const;
	bool postcalC(
		const Info_Step & info_step,
		const matrix & DVC,				// iw1/iw3, \mu1, iw2/iw4
		const size_t range_iw13,
		const size_t range_mu1,
		const size_t range_iw24,
		const size_t iw13H) const;
	bool postcalD( const matrix & H ) const { return (H.absmax()>threshold) ? true : false; }
	
private:
	// max_j \sqrt{ \sum_i m(i,j)^2 }
	double cal_matrix_inner_max( const matrix & m, const size_t ni, const size_t nj ) const;
	// max_i \sqrt{ \sum_j m(i,j)^2 }
	double cal_matrix_outer_max( const matrix & m, const size_t ni, const size_t nj ) const;
	// m m^+
	matrix m_mT( const matrix & m ) const;

private:

	bool flag_screen_cauchy = false;
	double threshold = 0;
	Abfs::Vector3_Order<int> Born_von_Karman_period;
	
private:

	// \sqrt{ || C C^+ || }		for C_{ I i \mu, K k }
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> C_norm4_outer_max;		// max_i \sqrt{ || C_i C_i^+ || }
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> C_norm4_inner_max;		// max_k \sqrt{ || C_k C_k^+ || }
	
	// \sqrt{ || V V^+ || }		for V_{ U \mu, V \nu }
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> V_norm4;
	
	// \sqrt{ || m m^+ || }		for D_{ K k, L l }
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> D_norm4_max;	// max_is \sqrt{ || m m^+ || }
	
	// || C ||					for C_{ I i \mu, K k }
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> C_norm2_outer_max;		// max_i || C_i ||
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> C_norm2_inner_max;		// max_k || C_k ||
		
//public:
//	static double num_screen1;
//	static double num_screen2;
//	static double num_screen3;
//	static double num_cal;
};

#endif
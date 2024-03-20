#ifndef ORB_GAUNT_TABLE_H
#define ORB_GAUNT_TABLE_H

#include <map>
#include "module_base/realarray.h"
#include "module_base/matrix.h"

class ORB_gaunt_table
{
	public:

	ORB_gaunt_table();
	~ORB_gaunt_table();

	/**
	 * Method 2: 
	 * using WIgner 3j symbols
	 * \f$ Y(l1,m1), Y(l2,m2), Y(L,M) \f$
	*/
	
	void init_Gaunt_CH(const int& Lmax);
	double Get_Gaunt_CH(
							const int& l1,
							const int& m1,
							const int& l2,
							const int& m2,
							const int& l3,
							const int& m3	);

	///M defined here are restricted within 0 to 2l+1
	///
	///should be transformed first
	double Get_Gaunt_SH(
							const int& l1,
							const int& mm1,
							const int& l2,
							const int& mm2,
							const int& l3,
							const int& mm3	);

	double Calc_Gaunt_CH(
							const int& l1,
							const int& m1,
							const int& l2,
							const int& m2,
							const int& l3,
							const int& m3	);
					

	/**
	 * MEthod 2
	 *
	 * Directly Calculate integral of
	 * \f$ S(l_1,m_1), S(l_2,m_2), S(L,M) \f$
	 */
	ModuleBase::realArray Gaunt_Coefficients;

	/// (1) Make Ylm_Gaunt Table.
	///----------------
    /*
	void init_Ylm_Gaunt(
		const int &lmax, 
		const double &s1,
		const double &e1,
		const double &s2,
		const double &e2);
    */

	/// (2) Use Ylm_Gaunt to calculate Gaunt Coefficinets element
	///------
    /*
	double Cal_Gaunt_single(
	   	const int &l1, 
		const int &m1, 
	   	const int &l2, 
		const int &m2, 
	   	const int &l, 
		const int &m, 
	   	const double &s1, 
		const double &e1,    
	   	const double &s2, 
		const double &e2);
    */

	/// (3) Make the whole Gaunt Coefficients table
	/// ------------------------------------
	void init_Gaunt(const int &lmax);

	static int get_lm_index(const int l, const int m);

	static int Index_M(const int& m);

	private:
	
	// Index Function
	// Yu's mehtod
	// Peize Lin delete void ModuleBase::GlobalFunc::ZEROS(); 2016-08-26
	
	//int P_EL(const int& L);

	int EP_EL(const int& L);

	int index_func(
			const int& l1,
			const int& l2,
			const int& l3,
			const int& m3	);
	
	double Fact(const int& n);

	void Swap(
			int& l1,
			int& m1,
			int& l2,
			int& m2	);

	//2*Lmax+1
	std::map<int,std::map<int,double>> Gaunt_CH;		// Peize Lin update 2016-08-26
	
	//direct integral
	ModuleBase::matrix Ylm_Gaunt;
};
#endif

//==========================================================
// AUTHOR : mohan
// DATE : 2009-3-29
// Last Modify : 2021-01-04
//==========================================================
#ifndef BESSEL_BASIS_H
#define BESSEL_BASIS_H
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/realarray.h"

#include "../module_cell/unitcell.h"

//==========================================================
// CLASS :
// NAME :  Bessel_Basis
//==========================================================
class Bessel_Basis
{
public:
	Bessel_Basis();
	~Bessel_Basis();

	/// @brief Initialization of Bessel function related matrices.
	/// @details Used for a specific group of C4 coefficients. 2021-01-04, mohan added a new input parameter lmax_in, if we only generate numerical atomic orbitals based on spherical Bessel functions, lmax_in = ucell.lmax. However, if we want to generate Spherical Bessel functions (SBF) for descriptor, then the lmax_in is controlled by user.
	/// @note This function is called in module_io/numerical_basis.cpp and module_io/numerical_descriptor.cpp
	/// @param start_from_file whether read C4 coefficients stored in external files
	/// @param ecutwfc cutoff for numerical atomic orbitals
	/// @param ntype atom types
	/// @param lmax_in maximal angular momentum for numerical orbitals
	/// @param smooth whether smooth SBFs when perform integration to calculate value of matrix element of TableOne. For details, see J. Phys.: Condens. Matter 22 (2010) 445501
	/// @param sigma stddev of Gaussian function for smoothing SBFs
	/// @param rcut_in cutoff radius for SBFs
	/// @param tol_in accurancy control for SBFs
	/// @param dk kspace grid
	/// @param dr realspace grid
	/// @param ucell UnitCell class object, ucell.nmax will be used in this function
	void init( 
		const bool start_from_file,
		const double &ecutwfc,
		const int &ntype,
		const int &lmax_in,
		const bool &smooth,
		const double &sigma,
		const double &rcut_in,
		const double &tol_in,
		const UnitCell& ucell,
		const double &dk = 0.01,
		const double &dr = 0.01
		);
	/// @brief return number of SBFs used for one `chi` (see details for more information)
	/// @details atomic orbital is constructed always with not only one set of SBFs. For different sets, they are marked with different `chi`(s), similar with concept of contracted GTOs. For one `chi`, it is 'q' the summation index, and q is in SBFs like: j_l(q*r), where l is the order of SBF.
	/// @return number of SBFs
	const int& get_ecut_number() const { return Ecut_number;}

	/// @brief Cubic spline interpolation for matrix Faln
	/// @param it atom type index
	/// @param l angular momentum
	/// @param ic chi index
	/// @param gnorm norm of G+k vector
	/// @return interpolated value
	double Polynomial_Interpolation(const int &it, const int &l, const int &ic, const double &gnorm)const;
	/// @brief Cubic spline interpolation for matrix TableOne
	/// @param l angular momentum
	/// @param ie q index (see explanation in note of function BesselBasis::get_ecut_number())
	/// @param gnorm norm of G+k vector
	/// @return interpolated value
	double Polynomial_Interpolation2(const int &l, const int &ie, const double &gnorm)const;

	
	/// @brief  get energy cutoff, which is used to truncate SBF Jlq. 
	/// @param  
	/// @return energy cutoff in Ry
	const double &get_ecut(void) const {return ecut;}
	/// @brief cutoff radius of radial SBF Jlq.
	/// @param  
	/// @return cutoff radius in a.u.
	const double &get_rcut(void) const {return rcut;}

	const double &get_tolerence(void) const {return tolerence;}


	/// @brief check if SBFs are smoothed (mohan add 2009-08-28)
	/// @attention in this case, the Jlq are not the true Jlq.
	/// @param  
	/// @return boolean whether SBFs are smoothed
	const bool &get_smooth(void) const {return smooth;}
	/// @brief get sigma the stddev (standard deviation) used in smooth function (Gaussian function)
	/// @param  
	/// @return stddev of smooth function
	const double &get_sigma(void) const {return sigma;}

private:
	/// @brief the most important array to calculate spillage, has dimension (ntype, lmax+1, max_n, nk)
	ModuleBase::realArray Faln;

	/// @brief Coefficients to be optimized!
	ModuleBase::realArray C4;

	/// @brief matrix whose elements are int{dr r^2 j_l(qr)*j_l(kr)}, has dimension (lmax+1, nq, nk)
	ModuleBase::realArray TableOne;

	/// @brief mesh of k vector, k is in j_l(k*r)
	int kmesh;
	/// @brief grid of k
	double Dk;
	/// @brief number of q vector, q is in j_l(q*r)
	int Ecut_number;
	/// @brief Cutoff radius (in a.u.) of SBFs, for any SBF j_l(qr), r>=rcut, j_l(q*r) = 0 (if not smoothed)
	double rcut;
	/// @brief energy cutoff for determining kmesh and number of SBFs
	double ecut;
	double tolerence;
	/// @brief whether smooth SBFs around cutoff radius, resulting in non-zero values. For importance of smooth of SBFs, see J. Phys.: Condens. Matter 22 (2010) 445501, eqn 6. (mohan add 2009-01-18)
	bool smooth;
	/// @brief stddev of smooth function (Gaussian function, centered at rcut)
	double sigma;

	/// @brief Allocate memory for C4 matrix and initialize all elements to one.
	/// @param ntype number of atom types
	/// @param lmax maximal angular momentum of localized orbitals
	/// @param nmax maximal principal quantum number of localized orbitals
	/// @param ecut_number number of SBFs
	void allocate_C4(
		const int &ntype,
		const int &lmax, 
		const int &nmax,
		const int &ecut_number,
		const UnitCell& ucell
		);

	/// @brief Read C4 from external file. Presently an O(N^2) search algorithm is used. A HTML parser is needed in the future to improve performance.
	/// @param name name of external file where C4-stored file information is contained
	/// @param ntype number of atom types
	/// @param ecut energy cutoff
	/// @param rcut cutoff radius
	/// @param ecut_number number of SBFs
	/// @param tolerence accurancy of SBFs, here only used for consistency check
	void readin_C4(
		const std::string &name,
		const int &ntype,
		const int &ecut,
		const int &rcut,
		const int &ecut_number,
		const double &tolerence,
		const UnitCell& ucell
		);

	void init_TableOne(void);

	/// @brief calculate F_{aln}(it, il, in, ik) = sum_{ie}{C4(it, il, in, ie)*TableOne(il, ie, ik)}, where TableOne is overlap integral between two spherical bessel functions (jle(r) and jlk(r))
	/// @param ntype number of atomtype
	/// @param lmax maximal angular momentum
	/// @param nmax maximal chi
	/// @param ecut_number number of SBFs
	void init_Faln(
		const int &ntype,
		const int &lmax,
		const int &nmax,
		const int &ecut_number,
		const UnitCell& ucell
		);
	
	/// @brief number of localized wave functions
	int nwfc;

	/// @brief calculate element value of TableOne matrix
	/// @details (be called in Bessel_Basis::init(), used for outputing overlap Q matrix) initialize the table whose matrix element is the result of integral int{dr r^2 jle(r)*jlk(r)}, TableOne has three subscript (l, ie, ik), the first runs over orbitals' angular momentum and ie, ik run over ecut_number and kmesh SBFs
	/// @param smooth_in whether jle(r) SBF is smoothed by a Gaussian function
	/// @param sigma_in stddev for controlling smearing of Gaussian function for smoothing jle(r)
	/// @param ecutwfc planewave kinetic energy cutoff for controlling kspace sampling
	/// @param rcut cutoff radius of SBFs
	/// @param dr realspace grid
	/// @param dk kspace grid
	/// @param lmax maximal angular momentum for SBFs
	/// @param ecut_number number of SBFs
	/// @param tolerence accurancy of SBFs
	void init_TableOne(
		const bool smooth_in,
		const double &sigma_in,
		const double &ecut,
		const double &rcut,
		const double &dr,
		const double &dk,
		const int &lmax,
		const int &ecut_number,
		const double &tolerence);
};

#endif

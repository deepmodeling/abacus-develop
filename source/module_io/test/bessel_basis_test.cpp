/// @details unit test for module_io/bessel_basis, not tested functions: readin_C4, allocate_C4

#include <iostream>
#include <fstream>

#include <cmath>

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include <cassert>

#include "../bessel_basis.h"
#include "../../module_cell/unitcell.h"
#include "../../module_elecstate/magnetism.h"

#ifdef __LCAO
#include "../../module_cell/setup_nonlocal.h"
#endif
#include "gtest/gtest.h"

/// @brief Simpson integral
/// @attention this function is a COMPLETE version, but there is no improvement in performance.
/// @param x variable stored in a vector
/// @param y function value stored in a vector
/// @return integral value
double SimpsonIntegral(const std::vector<double> &x, const std::vector<double> &y)
{
    double result = 0.0;
    result += (y[0] + y[y.size() - 1]);
    /* x and y must have the same size and their length must be a same odd number */
    assert(x.size() == y.size());
    assert(x.size() % 2 == 1);
    double h = x[1] - x[0];
    for (int i = 1; i < x.size() - 1; i+=2)
    {
        result += (2*y[i] + 4*y[i+1]);
    }
    result *= h/3;
    return result;
}
/// @brief recursive definition of 1st-kind Spherical Bessel function of the first kind
/// @attention this function is a COMPLETE version, can replace the one in basic mathematical library
/// @param l order of 1st-kind spherical Bessel function
/// @param x variable
/// @return value of 1st-kind spherical Bessel function
double SphericalBessel(int l, double x) {
    /* c++ cannot handle limits well */
    if (x == 0) {
        if (l == 0) {
            return 1;
        }
        else {
            return 0;
        }
    }
    else {
        if (l == 0) {
            return std::sin(x) / x;
        }
        else if (l == 1) {
            return (std::sin(x) / (x * x)) - std::cos(x) / x;
        }
        else {
            return ((2 * l - 1) / x) * SphericalBessel(l - 1, x) - SphericalBessel(l - 2, x);
        }
    }
}
/// @brief 1st-kind Spherical Bessel function of the first kind mapped on realspace grid
/// @attention this function is a COMPLETE version, can replace the one in basic mathematical library
/// @param l order of 1st-kind spherical Bessel function
/// @param q wave vector
/// @param r realspace grid
/// @return function value on realspace grid
std::vector<double> CalculateSphericalBessel(int l, double q, const std::vector<double>& r) {
    std::vector<double> bessel;
    for (const auto& radius : r) {
        double x = q * radius;
        double besselValue = SphericalBessel(l, x);
        bessel.push_back(besselValue);
    }
    return bessel;
}
/// @brief Get the first p zeros of 1st-kind, l-ordered Spherical Bessel function from a constant table
/// @attention this function is a INCOMPLETE version, due to limited support of numerical table.
/// @param order l, angular momentum
/// @param number p, number of zeros needed
/// @return a vector of q, the q is from j_l(qr)
std::vector<double> GetSphericalBesselZeros(int order, int number) {
    std::map<std::pair<int, int>, double> zeros;
    
    zeros[{0, 1}] = 3.14159;  zeros[{0, 2}] = 6.28318;  zeros[{0, 3}] = 9.42477;
    zeros[{0, 4}] = 12.5664;  zeros[{0, 5}] = 15.708;   zeros[{0, 6}] = 18.8495;
    zeros[{0, 7}] = 21.9911;  zeros[{0, 8}] = 25.1327;  zeros[{0, 9}] = 28.2743;
    zeros[{0, 10}] = 31.4159;
    
    zeros[{1, 1}] = 4.49341;  zeros[{1, 2}] = 7.72525;  zeros[{1, 3}] = 10.9041;
    zeros[{1, 4}] = 14.0662;  zeros[{1, 5}] = 17.2208;  zeros[{1, 6}] = 20.3713;
    zeros[{1, 7}] = 23.5181;  zeros[{1, 8}] = 26.6617;  zeros[{1, 9}] = 29.8029;
    zeros[{1, 10}] = 32.9425;
    
    zeros[{2, 1}] = 5.76346;  zeros[{2, 2}] = 9.09501;  zeros[{2, 3}] = 12.3229;
    zeros[{2, 4}] = 15.5146;  zeros[{2, 5}] = 18.6861;  zeros[{2, 6}] = 21.8457;
    zeros[{2, 7}] = 24.9989;  zeros[{2, 8}] = 28.1498;  zeros[{2, 9}] = 31.2997;
    zeros[{2, 10}] = 34.4491;
    
    zeros[{3, 1}] = 7.01559;  zeros[{3, 2}] = 10.4013;  zeros[{3, 3}] = 13.5821;
    zeros[{3, 4}] = 16.7496;  zeros[{3, 5}] = 19.9023;  zeros[{3, 6}] = 23.0446;
    zeros[{3, 7}] = 26.1799;  zeros[{3, 8}] = 29.3105;  zeros[{3, 9}] = 32.4377;
    zeros[{3, 10}] = 35.5629;
    
    zeros[{4, 1}] = 8.26356;  zeros[{4, 2}] = 11.6209;  zeros[{4, 3}] = 14.7965;
    zeros[{4, 4}] = 17.9598;  zeros[{4, 5}] = 21.113;   zeros[{4, 6}] = 24.2583;
    zeros[{4, 7}] = 27.3979;  zeros[{4, 8}] = 30.5325;  zeros[{4, 9}] = 33.6635;
    zeros[{4, 10}] = 36.7914;
    
    zeros[{5, 1}] = 9.51045;  zeros[{5, 2}] = 12.8377;  zeros[{5, 3}] = 16.0106;
    zeros[{5, 4}] = 19.1714;  zeros[{5, 5}] = 22.3224;  zeros[{5, 6}] = 25.4666;
    zeros[{5, 7}] = 28.6055;  zeros[{5, 8}] = 31.7408;  zeros[{5, 9}] = 34.873;
    zeros[{5, 10}] = 38.0025;
    
    std::vector<double> result;
    for (int i = 1; i <= number; ++i) {
        result.push_back(zeros[{order, i}]);
    }
    return result;
}
/// @brief Get mod of q vector of Spherical Bessel functions, all q satisfy when r=`rcut`, j_l(qr)=0. 
/// @details first solve the equation j_l(x) = 0, therefore get the table (l, k) -> x, where l is the order of SBF and k is the k-th zero of j_l(x). Then let x = q*rcut, therefore q = x/rcut, return it.
/// @attention this function itself is a COMPLETE version, while the function it called GetSphericalBesselZeros may be INCOMPLETE, due to limited support of numerical table.
/// @param order the angular momentum of Spherical Bessel functions
/// @param number number of q to be returned
/// @param rcut 'cutoff radius' of Spherical Bessel functions. When r=rcut, Spherical Bessel functions are zero, and for r>rcut, they are zero as required by concept of constructing localized atomic orbital.
/// @return a vector of q, the q is from j_l(qr)
std::vector<double> GetqList(int order, int number, double rcut) {
    std::vector<double> qList;
    std::vector<double> zerosList = GetSphericalBesselZeros(order, number);
    for (const auto& zero : zerosList) {
        qList.push_back(zero / rcut);
    }
    return qList;
}
/// @brief overload operator * for vectors, elements will be multiplied one by one
/// @param x one vector
/// @param y another vector
/// @return vector after multiplication
std::vector<double> operator*(const std::vector<double>& x, const std::vector<double>& y) {
    std::vector<double> result;
    for (int i = 0; i < x.size(); ++i) {
        result.push_back(x[i] * y[i]);
    }
    return result;
}
/// @brief init_TableOne for unit test
/// @attention this is a COMPLETE version of init_TableOne, can replace the one in module_io/bessel_basis.cpp
/// @param smooth whether to smooth the function (gaussian function)
/// @param sigma stddev of gaussian function for smoothing
/// @param ecutwfc control the number of Spheical Bessel functions
/// @param rcut cutoff radius of Spherical Bessel functions, r>rcut, Spherical Bessel functions are zero
/// @param lmax maximal angular momentum of Spherical Bessel functions
/// @param dr grid spacing of r
/// @param dk grid spacing of k
/// @return `std::vector<std::vector<std::vector<double>>>` TableOne[angular momentum][index of Spherical Bessel function][Arbitrary wave vector k]
std::vector<std::vector<std::vector<double>>> GenerateTableOne(const bool smooth, const double sigma, const double ecutwfc, const double rcut, const int lmax, const double dr, const double dk){
    std::vector<double> rGrid;
    std::vector<double> SmoothFactor_rGrid;
    int rGridNum = static_cast<int>(rcut/dr) + 4;
    if (rGridNum % 2 == 0)
    {
        rGridNum += 1;
    }
    for (int indexr=0; indexr<rGridNum; indexr++)
    {
        rGrid.push_back(indexr*dr);
        SmoothFactor_rGrid.push_back(1.0-exp(-std::pow((indexr*dr-rcut),2)/std::pow(sigma,2)));
    }

    std::vector<double> kGrid;
    int kGridNum = static_cast<int>(sqrt(ecutwfc)/dk) + 4 + 1;
    if (kGridNum % 2 == 0)
    {
        kGridNum += 1;
    }
    for (double indexk=0; indexk<kGridNum; indexk++)
    {
        kGrid.push_back(indexk*dk);
    }

    int NumBesselFunction = static_cast<int>(sqrt(ecutwfc)*rcut/M_PI);

    std::vector<std::vector<std::vector<double>>> TableOne;

    for (int l=0; l<lmax+1; l++)
    {
        std::vector<double> qList = GetqList(l, NumBesselFunction, rcut);
        /* should initialize TableOne first */
        if (l == 0)
        {
            for (int l=0; l<lmax+1; l++)
            {
                std::vector<std::vector<double>> Zeros2d;
                for (int indexq=0; indexq<qList.size(); indexq++)
                {
                    std::vector<double> Zeros1d;
                    for (int indexk=0; indexk<kGrid.size(); indexk++)
                    {
                        Zeros1d.push_back(0.0);
                    }
                    Zeros2d.push_back(Zeros1d);
                }
                TableOne.push_back(Zeros2d);
            }
        }
        for (int indexq=0; indexq<qList.size(); indexq++)
        {
            std::vector<double> jle_rGrid = CalculateSphericalBessel(l, qList[indexq], rGrid);
            if (smooth)
            {
                jle_rGrid = jle_rGrid*SmoothFactor_rGrid;
            }
            for (int indexk=0; indexk<kGrid.size(); indexk++)
            {
                std::vector<double> jlk_rGrid = CalculateSphericalBessel(l, kGrid[indexk], rGrid);
                std::vector<double> function_rGrid = jlk_rGrid*jle_rGrid*rGrid*rGrid;
                TableOne[l][indexq][indexk] = SimpsonIntegral(rGrid, function_rGrid);
            }
        }
    }
    return TableOne;
}
/// @brief Improved version of module_io/bessel_basis::readin_C4 and allocate_C4 functions, for generating C4 matrix but now with a higher speed on accessing elements
/// @attention function will read in total number of chi-s both from file and from input, and assert that they are the same. It is also just a INCOMPLETE version, for a complete version, HTML parser library will be included and the other parameter, NumAtomType, will also be used for calibrating data.
/// @param FileName name of external file where C4-stored file information is contained
/// @param NumAtomType number of atom types
/// @param l maximal angular momentum of localized orbitals
/// @param NumChi number of contracted spherical bessel functions to fit one atomic orbital
/// @param NumBesselFunction number of spherical bessel functions in one contracted spherical bessel function
/// @return a map whose key is a string of "atomtype l chi q", and value is the corresponding C4 value
std::unordered_map<std::string, double> ReadinC4(const std::string &FileName, const int &NumAtomType, const int &l, const int &NumChi, const int &NumBesselFunction){
    std::unordered_map<std::string, double> C4Map;
    std::ifstream C4File(FileName);
    /* plan to add an HTML parser library in the future... */
    std::string word;
    bool b_ReadC4Line = false;
    int TotalNumChi = 0;

    assert(!C4File.fail());
    while (C4File.good())
    {
        if (!b_ReadC4Line){
            C4File >> word;
            if (word == "<C4>")
            {
                b_ReadC4Line = true;
                C4File >> word; /* n_chi value read */
                TotalNumChi = std::stoi(word);
                assert(TotalNumChi == NumChi);
                continue;
            }
        }
        else{
            for (int indexchi = 0; indexchi < TotalNumChi; indexchi++){

                C4File >> word; 
                C4File >> word; 
                C4File >> word; /* skip title1, 2 and 3 */

                C4File >> word; std::string key = word; key += " ";
                C4File >> word; key += word; key += " "; 
                C4File >> word; key += word; key += " ";

                for (int indexNumBesselFunction = 0; indexNumBesselFunction < NumBesselFunction; indexNumBesselFunction++)
                {
                    std::string keyForMap = key + std::to_string(indexNumBesselFunction);
                    C4File >> word;
                    C4Map[keyForMap] = std::stod(word);
                }
            }
            break;
        }
    }
    C4File.close();
    return C4Map;
}
/// @brief Generate F_{alpha, l, chi, k} matrix
/// @param FileName name of external file where C4-stored file information is contained
/// @param NumAtomType number of atom types
/// @param l maximal angular momentum of localized orbitals
/// @param NumChi number of contracted spherical bessel functions to fit one atomic orbital
/// @param NumBesselFunction number of spherical bessel functions in one contracted spherical bessel function
/// @param vvv_d_TableOne Integral table, has subscripts (l, q, k), whose element is the result of integral int{dr r^2 jle(r)*jlk(r)}. l runs over angular momentum, q runs over all spherical bessel functions, k runs over k points sampled controlled by ecutwfc and dk by sqrt(ecutwfc)/dk + 1 + 4.
/// @return F_{alpha, l, chi, k} matrix, whose element is the result of sum_{q}{C4(alpha, l, chi, q)*TableOne(l, q, k)}
std::vector<std::vector<std::vector<std::vector<double>>>> GenerateFaln(const std::string &FileName, const int &NumAtomType, const int &lmax, const int &NumChi, const int &NumBesselFunction, std::vector<std::vector<std::vector<double>>> vvv_d_TableOne){
    std::unordered_map<std::string, double> umap_str_d_C4Map = ReadinC4(FileName, NumAtomType, lmax, NumChi, NumBesselFunction);
    std::vector<std::vector<std::vector<std::vector<double>>>> vvvv_d_Faln;
    for (int indexAtomType = 0; indexAtomType < NumAtomType; indexAtomType++)
    {
        std::vector<std::vector<std::vector<double>>> vvv_d_Faln;
        for (int indexl = 0; indexl < lmax+1; indexl++)
        {
            std::vector<std::vector<double>> vv_d_Faln;
            for (int indexChi = 0; indexChi < NumChi; indexChi++)
            {
                std::vector<double> v_d_Faln;
                for (int indexk = 0; indexk < vvv_d_TableOne[0][0].size(); indexk++)
                {
                    double Faln = 0.0;
                    for (int indexBesselFunction = 0; indexBesselFunction < NumBesselFunction; indexBesselFunction++)
                    {
                        std::string key = std::to_string(indexAtomType) + " " + std::to_string(indexl) + " " + std::to_string(indexChi) + " " + std::to_string(indexBesselFunction);
                        Faln += umap_str_d_C4Map[key]*vvv_d_TableOne[indexl][indexBesselFunction][indexk];
                    }
                    v_d_Faln.push_back(Faln);
                }
                vv_d_Faln.push_back(v_d_Faln);
            }
            vvv_d_Faln.push_back(vv_d_Faln);
        }
        vvvv_d_Faln.push_back(vvv_d_Faln);
    }
    return vvvv_d_Faln;
}

/* OVERLOAD (de)constructors... */
pseudo_nc::pseudo_nc()
{
}
pseudo_nc::~pseudo_nc()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}

#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
/* OVERLOAD printM3 function? */
/*
void output::printM3(std::ofstream &ofs, const std::string &description, const ModuleBase::Matrix3 &m)
{
}
*/

/************************************************
 *       unit test of bessel_basis.cpp
 ***********************************************/

/*
 * - Tested Functions:
 *   - InitTest: Bessel_Basis::init(start_from_file, ecutwfc, ntype, lmax_in, smooth, sigma, rcut_in, tol_in, ucell, dk, dr)
 *     - initialize Bessel_Basis class object with following steps:
 *     - call init_TableOne(...): calculate TableOne (l, ie, ik) = int{dr r^2 jle(r)*jlk(r)} (jle(r) and jlk(r) are Spherical Bessel functions)
 *     - return, unless set start_from_file = true, will also do the following:
 *     - call allocate_C4(...): allocate memory for C4 4d-matrix (it, il, in, ie) and set all elements to 1.0
 *     - call readin_C4(...): read C4 value from external file and store in C4 4d-matrix
 *     - call init_Faln(...): calculate F_{aln}(it, il, in, ik) = sum_{ie}{C4(it, il, in, ie)*TableOne(il, ie, ik)}
 *   - PolynomialInterpolation2Test: Bessel_Basis::Polynomial_Interpolation2(const int &l, const int &ie, const double &gnorm)
 *     - return (cubic spline) interpolated element value of TableOne 3d-matrix
 *   - PolynomialInterpolationTest: Bessel_Basis::Polynomial_Interpolation(const int &it, const int &l, const int &ic, const double &gnorm)
 *     - return (cubic spline) interpolated element value of Faln 4d-matrix
 */

#define private public
class TestBesselBasis : public ::testing::Test {
protected:
    UnitCell ucell;
    Bessel_Basis besselBasis;

    void SetUp() override {
        ucell.ntype = 1;
        ucell.lmax = 0;
        ucell.nmax = 1;

        ucell.atoms = new Atom[1];
        ucell.atoms[0].l_nchi = new int[1];

        ucell.atoms[0].nwl = 0;
        ucell.atoms[0].l_nchi[0] = 1;
        /* setup_cell manually */
    }
    void TearDown() override {
        delete[] ucell.atoms->l_nchi;
        delete[] ucell.atoms;
    }
};

 TEST_F(TestBesselBasis, InitTest) {
     /* parameters required by Bessel_Basis::init() */
    bool b_TestFaln = false;
    double d_EnergyCutoff = 4;
    int i_Ntype = 1;
    int i_Lmax = 0;
    bool b_Smooth = false;
    double d_SmoothSigma = 0.1;
    double d_CutoffRadius = 2;
    double d_Tolerance = 0.01;
    double d_dk = 0.01;
    double d_dr = 2;
    /* number of Spherical Bessel functions with given order l,
      used to fit one chi */
    int i_Nq = static_cast<int>(sqrt(d_EnergyCutoff)*d_CutoffRadius/M_PI);
    /* number of SBF is expected to be 1 */
     besselBasis.init(
         b_TestFaln, d_EnergyCutoff, i_Ntype, i_Lmax, b_Smooth, 
         d_SmoothSigma, d_CutoffRadius, d_Tolerance, 
         ucell, d_dk, d_dr
     );
     EXPECT_EQ(besselBasis.get_ecut_number(), i_Nq);
     EXPECT_EQ(besselBasis.get_ecut(), d_EnergyCutoff);
     EXPECT_EQ(besselBasis.get_rcut(), d_CutoffRadius);
     EXPECT_EQ(besselBasis.get_tolerence(), d_Tolerance);
     EXPECT_EQ(besselBasis.get_smooth(), b_Smooth);
     EXPECT_EQ(besselBasis.get_sigma(), d_SmoothSigma);
}

TEST_F(TestBesselBasis, PolynomialInterpolation2Test) {
    /* function Bessel_Basis::Polynomial_Interpolation2 is to do
    cubic interpolation on TableOne, this matrix has, dimension l*nq*nk */
    /* parameters required by Bessel_Basis::init() */
    bool b_TestFaln = false;
    double d_EnergyCutoff = 4;
    int i_Ntype = 1;
    int i_Lmax = 0;
    bool b_Smooth = false;
    double d_SmoothSigma = 0.1;
    double d_CutoffRadius = 2.0;
    double d_Tolerance = 0.01;
    double d_dk = 2.0;
    double d_dr = 0.01; /* for d_dr will largely impact accurancy of numerical integration, will not give a toy-value */

    /* number of angular momentum considered should be 1 (lmax = 0)->The first dimension */
    int i_Nq = static_cast<int>(sqrt(d_EnergyCutoff)*d_CutoffRadius/M_PI);
    /* number of SBF is expected to be 1->The second dimension */
    int i_kMesh = static_cast<int>(sqrt(d_EnergyCutoff)/d_dk) + 4 + 1;
    /* the third dimension is at least to be 6 */
    /* therefore the expected dimension of TableOne is 1*1*6 */

    besselBasis.init(
        b_TestFaln, d_EnergyCutoff, i_Ntype, i_Lmax, b_Smooth, 
        d_SmoothSigma, d_CutoffRadius, d_Tolerance, 
        ucell, d_dk, d_dr
    );
    /* gnorm for interpolation */
    double d_Gnorm = 1.0;
    double d_position = d_Gnorm/d_dk;
    int i_position = static_cast<int>(d_position);
    assert(i_position < i_kMesh-4);
    double d_x0 = d_position - static_cast<double>(i_position);
    double d_x1 = 1.0 - d_x0;
    double d_x2 = 2.0 - d_x0;
    double d_x3 = 3.0 - d_x0;

   std::vector<std::vector<std::vector<double>>> vvv_d_TableOne = GenerateTableOne(
        b_Smooth, d_SmoothSigma, d_EnergyCutoff, d_CutoffRadius, 
        i_Lmax, d_dr, d_dk
    );
    double d_yExpected = vvv_d_TableOne[0][0][i_position]*d_x1*d_x2*d_x3/6.0+
                         vvv_d_TableOne[0][0][i_position]*d_x0*d_x2*d_x3/2.0-
                         vvv_d_TableOne[0][0][i_position]*d_x1*d_x0*d_x3/2.0+
                         vvv_d_TableOne[0][0][i_position]*d_x1*d_x2*d_x0/6.0;
    double d_yTested = besselBasis.Polynomial_Interpolation2(0, 0, d_Gnorm);
    EXPECT_NEAR(d_yExpected, d_yTested, 0.01);
}
TEST_F(TestBesselBasis, PolynomialInterpolationTest) {
    /* function Bessel_Basis::Polynomial_Interpolation is to do
    cubic interpolation on Faln, this matrix has, dimension atomtype*l*nchi*nk.
    To obtain Faln matrix, it is needed to first get C4 from external file.
    The C4 is coefficent of SBF, has dimension atomtype*l*nchi*nq.
    Therefore the Faln, equals the contraction of nq between C4 and TableOne.
    C4: atomtype*l*nchi*nq
    TableOne: l*nq*nk -> Faln: atomtype*l*nchi*nk */

    /* parameters required by Bessel_Basis::init() */
    bool b_TestFaln = true;
    double d_EnergyCutoff = 4;
    int i_Ntype = 1;
    int i_Lmax = 0;
    bool b_Smooth = false;
    double d_SmoothSigma = 0.1;
    double d_CutoffRadius = 2.0;
    double d_Tolerance = 0.01;
    double d_dk = 2.0;
    double d_dr = 0.01; /* for d_dr will largely impact accurancy of numerical integration, will not give a toy-value */

    int i_Nq = static_cast<int>(sqrt(d_EnergyCutoff)*d_CutoffRadius/M_PI);
    int i_kMesh = static_cast<int>(sqrt(d_EnergyCutoff)/d_dk) + 4 + 1;
     /*
      manipulate Bessel_Basis::init_Faln function 
      because for(int it=0; it<ntype; it++)
      then        for(int il=0; il<GlobalC::ucell.atoms[it].nwl+1; il++)
      and             for(int in=0; in<GlobalC::ucell.atoms[it].l_nchi[il]; in++)
      Parameter ntype is controlled by input, but
      Parameter nwl, l_nchi are controlled by GlobalC, and will determine exact dimension of Faln
    */
    /* 
      Therefore Faln is expected to have dimension 1*1*1*6
    */
    
    besselBasis.init(
        b_TestFaln, d_EnergyCutoff, i_Ntype, i_Lmax, b_Smooth, 
        d_SmoothSigma, d_CutoffRadius, d_Tolerance, 
        ucell, d_dk, d_dr
    );
    /* Gnorm for interpolation */
    double d_Gnorm = 1.0;
    double d_position = d_Gnorm/d_dk;
    int i_position = static_cast<int>(d_position);
    assert(i_position < i_kMesh-4);
    double d_x0 = d_position - static_cast<double>(i_position);
    double d_x1 = 1.0 - d_x0;
    double d_x2 = 2.0 - d_x0;
    double d_x3 = 3.0 - d_x0;

    std::vector<std::vector<std::vector<double>>> vvv_d_TableOne = GenerateTableOne(
        b_Smooth, d_SmoothSigma, d_EnergyCutoff, d_CutoffRadius, 
        i_Lmax, d_dr, d_dk
    );
    std::vector<std::vector<std::vector<std::vector<double>>>> vvvv_d_Faln = GenerateFaln(
        "./support/BesselBasis_UnitTest_C4_AtomType0.html", i_Ntype, i_Lmax, 1, i_Nq, vvv_d_TableOne
    );
    double d_yExpected = vvvv_d_Faln[0][0][0][i_position]*d_x1*d_x2*d_x3/6.0+
                         vvvv_d_Faln[0][0][0][i_position]*d_x0*d_x2*d_x3/2.0-
                         vvvv_d_Faln[0][0][0][i_position]*d_x1*d_x0*d_x3/2.0+
                         vvvv_d_Faln[0][0][0][i_position]*d_x1*d_x2*d_x0/6.0;
    double d_yTested = besselBasis.Polynomial_Interpolation(0, 0, 0, d_Gnorm);
    EXPECT_NEAR(d_yExpected, d_yTested, 0.01);
}
#undef private

#include "gtest/gtest.h"
#include "source_base/math_integral.h"
#include "source_base/global_variable.h"
#include <algorithm>
#include <string>
#include <vector>
#include <numeric>
#include <iomanip>

#ifdef __MPI
#include <mpi.h>
#endif

#include "source_basis/module_ao/orb_atomic_lm.h"

/***********************************************************
 *      unit test of class "Numerical_Orbital_Lm"
 ***********************************************************/

/* Tested functions:
 *
 * - set_orbital_info
 *   copies the input into class members,
 *   computes psi(k) from psi(r) (or psi(r) from psi(k) if input is psi(k)),
 *   interpolates psi(r) on a uniform grid,
 *   and optionally saves all psi data to file
 *
 * - cal_kradial
 *   applies a radial Fourier transform to psi(r) to get psi(k)
 *
 * - cal_kradial_sbpool
 *   same as cal_kradial (with a more efficient algorithm)
 *
 * - cal_rradial_sbpool
 *   applies an inverse radial Fourier transform to psi(k) to get psi(r)
 *
 * - plot
 *   save psi, psik, psi_uniform & dpsi_uniform to file
 *
 * - extra_uniform
 *   interpolates psi on a uniform grid to get psi_uniform, and calculates
 *   its first order derivative dpsi_uniform
 *   FIXME what is zty supposed to be?
 *
 * - all "getters"
 *   get access to class members
 *
 ***********************************************************/

class NumericalOrbitalLmTest : public ::testing::Test
{
protected:
    void SetUp();
    void TearDown();

    // objects under unit test
    std::vector<Numerical_Orbital_Lm> nolm_;

    // helper functions
    void init();
    void pour_data();
    void init_with_different_k(double const& ecut, double const& dk);
    size_t calc_nk(double const& ecutwfc, double const& dk);
    size_t calc_nr_uniform(double const& rcut, double const& dr_uniform);
    bool check_file_match(size_t const& nline, double const* col1, double const* col2, double const& tol, std::string const& fname);

    // Numerical_Orbital_Lm declares this fixture a friend, but a TEST_F body
    // lives in a class derived from it and friendship is not inherited, so the
    // calls into the private radial transforms are routed through these.
    void cal_kradial(size_t i) { nolm_[i].cal_kradial(); }
    void cal_kradial_sbpool(size_t i) { nolm_[i].cal_kradial_sbpool(); }
    void cal_rradial_sbpool(size_t i) { nolm_[i].cal_rradial_sbpool(); }
    void plot(size_t i) const { nolm_[i].plot(); }

    // radial real-space mesh spacing
    double dr_;

    // data to pour into nolm_ (used by set_orbital_info)
    std::string elem_label_;
    int index_atom_type_;
    std::vector<int> l_;
    std::vector<int> index_chi_l_; // index of orbital within its L
    int nr_;
    double* r_radial_;
    double* rab_;
    Numerical_Orbital_Lm::Psi_Type psi_type_;
    std::vector<double*> chi_;
    size_t nk_;
    double dk_;
    double dr_uniform_;
    bool flag_plot_;
    bool flag_sbpool_;
    bool force_flag_;
};


size_t NumericalOrbitalLmTest::calc_nk(double const& ecutwfc, double const& dk) {

    // current formula for calculating nk from ecutwfc & dk
    // see source_basis/module_ao/orb_read.cpp, function "Read_Orbitals"

    size_t nk = 0;

    if(ecutwfc < 20) {
        nk = static_cast<int>( 2 * sqrt(ecutwfc) / dk )  + 4;
    } else {
        nk = static_cast<int>( sqrt(ecutwfc) / dk )  + 4;
    }

    if (nk%2 == 0) {
        ++nk;
    }

    return nk;
}


size_t NumericalOrbitalLmTest::calc_nr_uniform(double const& rcut, double const& dr_uniform) {
    return static_cast<int>(rcut/dr_uniform) + 10;
}


void NumericalOrbitalLmTest::SetUp() {

    ///////////////////////////////////////////////////
    //                  Parameters
    ///////////////////////////////////////////////////

    // numerical atomic orbital file to read data from
    std::string orb_file = "./lcao_H2O/O_gga_7au_60Ry_2s2p1d.orb";

    // energy cutoff for 1D integration (radial Fourier transform)
    // this number is used in most tests except r2k2r_consistency
    // where a much larger value is used.
    double ecutwfc = 100.0;

    // In normal ABACUS calculation, dk, nk & dr_uniform
    // are retrieved from the LCAO_Orbitals object,
    // see setupNonlocal() in source_cell/setup_nonlocal.cpp.
    // Here we just provide some reasonable values.
    dk_ = 0.01;
    dr_uniform_ = 0.001;

    nk_ = calc_nk(ecutwfc, dk_);

    // not really meaningful in this unit test
    index_atom_type_ = 42;

    // plot() is tested alone rather than called in set_orbital_file
    flag_plot_ = false;

    // in agreement with the current code
    flag_sbpool_ = true;

    // PARAM.input.cal_force
    // if true, extra_uniform will compute zty
    force_flag_ = true;


    ///////////////////////////////////////////////////
    //              Read orbital file
    ///////////////////////////////////////////////////

    std::ifstream ifs(orb_file);

    // variables read from orb_file
    double ecut = 0.0; // energy cutoff in orbital file, not the one in INPUT
    double rcut = 0.0;
    unsigned lmax = 0;
    std::vector<unsigned> nchi_l; // number of chi of each angular momentum
    int nr_read = 0; // number of mesh points

    // nr_ has to be odd, probably due to the way Simpson_Integral works
    // nr_ equals nr_read if nr_read is odd,
    // and equals nr_read+1 if nr_read is even

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "Element");
    ifs >> elem_label_;

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "Cutoff(Ry)");
    ifs >> ecut;

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "Cutoff(a.u.)");
    ifs >> rcut;

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "Lmax");
    ifs >> lmax;

    // number of chi for each angular momentum
    nchi_l.resize(lmax+1);
    std::vector<std::string> symbol = {"S", "P", "D", "F", "G", "H", "I", "K"};
    for (size_t l = 0; l <= lmax; ++l) {
        std::string key = symbol[l] + "orbital-->";
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, key);
        ifs >> nchi_l[l];
    }

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "Mesh");
    ifs >> nr_read;

    nr_ = (nr_read%2) ? nr_read : nr_read+1;

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "dr");
    ifs >> dr_;

    size_t nchi_tot = std::accumulate(nchi_l.begin(), nchi_l.end(), 0);

    l_.resize(nchi_tot);
    index_chi_l_.resize(nchi_tot);
    chi_.resize(nchi_tot);
    nolm_.resize(nchi_tot);

    r_radial_ = new double[nr_];
    rab_ = new double[nr_];
    for (int i = 0; i != nr_; ++i) {
        r_radial_[i] = i*dr_;
        rab_[i] = dr_;
    }

    // used in the normalization of input chi
    double* integrand= new double[nr_];
    for (int ir = 0; ir != nr_; ++ir) {
        integrand[ir] = 0.0;
    }

    // the orbital file contains psi(r)
    psi_type_ = Numerical_Orbital_Lm::Psi_Type::Psi;

    size_t ichi_tot = 0;

    for (size_t l = 0; l <= lmax; ++l) {
        for (size_t ichi_l = 0; ichi_l < nchi_l[l]; ++ichi_l) {

            // the next block of the orbital file is like
            // Type    L    N
            //    0    0    0
            std::string dummy;
            unsigned type = 0, L = 0, N = 0;
            ifs >> dummy >> dummy >> dummy;
            ifs >> type >> L >> N;
            assert(L == l);
            assert(N == ichi_l);
            //std::cout << "l = " << l << "   index_chi_l = " << N << std::endl;

            l_[ichi_tot] = l;
            index_chi_l_[ichi_tot] = ichi_l;

            // read & normalize chi
            chi_[ichi_tot] = new double[nr_];

            for (int ir = 0; ir != nr_; ++ir) {
                chi_[ichi_tot][ir] = 0.0;
            }

            for (int ir = 0; ir != nr_read; ++ir) {
                ifs >> chi_[ichi_tot][ir];
                integrand[ir] = std::pow(chi_[ichi_tot][ir]*r_radial_[ir], 2);
            }

            // radint = \int_0^{\infty} dr r^2 [chi(r)]^2
            double radint = 0.0;
            ModuleBase::Integral::Simpson_Integral(nr_, integrand, rab_, radint);

            for (int ir = 0; ir != nr_; ++ir) {
                chi_[ichi_tot][ir] /= std::sqrt(radint);

            }

            ++ichi_tot;

        }
    }

    delete[] integrand;
}


void NumericalOrbitalLmTest::TearDown() {
    delete[] r_radial_;
    delete[] rab_;
    for (size_t ichi = 0; ichi != chi_.size(); ++ichi) {
        delete[] chi_[ichi];
    }
}


void NumericalOrbitalLmTest::init() {

    // initialized the tested objects by pouring
    // the data collected in SetUp() into nolm_

    for (size_t ichi_tot = 0; ichi_tot != chi_.size(); ++ichi_tot) {
        nolm_[ichi_tot].set_orbital_info(elem_label_, index_atom_type_,
                l_[ichi_tot], index_chi_l_[ichi_tot], nr_, rab_,
                r_radial_, psi_type_, chi_[ichi_tot], nk_, dk_,
                dr_uniform_, flag_plot_, flag_sbpool_, force_flag_);
    }
}


void NumericalOrbitalLmTest::init_with_different_k(double const& ecut, double const& dk) {

    // initialize nolm_ with a k mesh specified by given ecut & dk

    this->dk_ = dk;
    this->nk_ = calc_nk(ecut, dk);
    this->init();

}


bool NumericalOrbitalLmTest::check_file_match(size_t const& nline, double const* col1, double const* col2, double const& tol, std::string const& fname) {

    /* This function checks whether the content of file named "fname"
     * contains certain data with certain format.
     *
     * The file should contain only two colums of floating-point numbers,
     * each column should match col1 or col2 within a tolerance of tol,
     * and the number of lines should be nline. No empty line is allowed.
     * No space is allowed except for the delimiter between the two columns.
     *
     * Return true if all the criteria are satisfied; false otherwise.
     */

    std::ifstream ifs;
    std::string tmp1, tmp2, tmp_line;

    ifs.open(fname);

    for (size_t i = 0; i != nline; ++i) {
        std::stringstream ss;
        std::getline(ifs, tmp_line);
        ss << tmp_line;
        ss >> tmp1 >> tmp2;

        if( std::abs( col1[i] - std::stod(tmp1.c_str()) ) > tol ||
            std::abs( col2[i] - std::stod(tmp2.c_str()) ) > tol ||
            ss.tellg() != -1 ) {
            return false;
        }
    }

    std::getline(ifs, tmp_line);
    if ( ifs.tellg() != -1 ) {
        return false;
    }

    ifs.close();

    return true;
}


TEST_F(NumericalOrbitalLmTest, Init) {

    // this test checks whether set_orbital_info works as expected

    // before
    // a brief check of the default constructor
    for (size_t ichi_tot = 0; ichi_tot != chi_.size(); ++ichi_tot) {
        EXPECT_EQ(nolm_[ichi_tot].getLabel(), "");
        EXPECT_EQ(nolm_[ichi_tot].getType(), 0);
        EXPECT_EQ(nolm_[ichi_tot].getL(), 0);
        EXPECT_EQ(nolm_[ichi_tot].getChi(), 0);
        EXPECT_EQ(nolm_[ichi_tot].getNr(), 1);
        EXPECT_EQ(nolm_[ichi_tot].getNk(), 1);
        EXPECT_EQ(nolm_[ichi_tot].nr_uniform, 1);

        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].getRcut(), 0.0);
        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].getKcut(), 0.0);
        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].getDk(), 0.0);
        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].dr_uniform, -1.0);
        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].zty, 0.0);

        EXPECT_TRUE(nolm_[ichi_tot].get_r_radial().empty());
        EXPECT_TRUE(nolm_[ichi_tot].get_k_radial().empty());
        EXPECT_TRUE(nolm_[ichi_tot].get_psi().empty());
        EXPECT_TRUE(nolm_[ichi_tot].get_psir().empty());
        EXPECT_TRUE(nolm_[ichi_tot].get_psif().empty());
        EXPECT_TRUE(nolm_[ichi_tot].get_psi_k().empty());
        EXPECT_TRUE(nolm_[ichi_tot].get_psi_k2().empty());
    }

    this->init();

    // after
    for (size_t ichi_tot = 0; ichi_tot != chi_.size(); ++ichi_tot) {
        EXPECT_EQ(nolm_[ichi_tot].getType(), index_atom_type_);
        EXPECT_EQ(nolm_[ichi_tot].getNk(), nk_);
        EXPECT_EQ(nolm_[ichi_tot].nr_uniform,
                this->calc_nr_uniform(nolm_[ichi_tot].getRcut(), dr_uniform_));

        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].getKcut(), (nk_-1)*dk_);
        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].getDk(), dk_);
        EXPECT_DOUBLE_EQ(nolm_[ichi_tot].dr_uniform, dr_uniform_);
        // TODO zty yet to be understood
        //EXPECT_DOUBLE_EQ(nolm_[ichi_tot].zty, 0.0);

        EXPECT_EQ(nolm_[ichi_tot].get_r_radial().size(), nr_);
        EXPECT_EQ(nolm_[ichi_tot].get_k_radial().size(), nk_);
        EXPECT_EQ(nolm_[ichi_tot].get_psi().size(), nr_);
        EXPECT_EQ(nolm_[ichi_tot].get_psir().size(), nr_);
        EXPECT_EQ(nolm_[ichi_tot].get_psif().size(), nk_);
        EXPECT_EQ(nolm_[ichi_tot].get_psi_k().size(), nk_);
        EXPECT_EQ(nolm_[ichi_tot].get_psi_k2().size(), nk_);

        for (int ir = 0; ir != nr_; ++ir) {
            EXPECT_DOUBLE_EQ(nolm_[ichi_tot].get_r_radial()[ir], ir*0.01);
            EXPECT_DOUBLE_EQ(nolm_[ichi_tot].get_psi()[ir], chi_[ichi_tot][ir]);
            EXPECT_DOUBLE_EQ(nolm_[ichi_tot].get_psir()[ir], ir*0.01*chi_[ichi_tot][ir]);
        }

        // whether psif makes sense or not is checked in r2k2r_consistency

        for (size_t ik = 0; ik != nk_; ++ik) {
            EXPECT_DOUBLE_EQ(nolm_[ichi_tot].get_k_radial()[ik], ik*dk_);
            EXPECT_DOUBLE_EQ(nolm_[ichi_tot].get_psi_k()[ik], ik*dk_*nolm_[ichi_tot].get_psif()[ik]);
            EXPECT_DOUBLE_EQ(nolm_[ichi_tot].get_psi_k2()[ik], ik*dk_*nolm_[ichi_tot].get_psi_k()[ik]);
        }
    }

    // see orb file for details

    double max_tol = 1e-12;

    EXPECT_EQ(nolm_[0].getLabel(), "O");
    EXPECT_EQ(nolm_[0].getL(), 0);
    EXPECT_EQ(nolm_[0].getChi(), 0);
    EXPECT_EQ(nolm_[0].getNr(), 701);
    EXPECT_DOUBLE_EQ(nolm_[0].getRcut(), 7.0);
    EXPECT_NEAR(nolm_[0].get_psi()[0], 1.208504975904e+00, max_tol);
    EXPECT_NEAR(nolm_[0].get_psi()[1], 1.208605373194e+00, max_tol);
    EXPECT_NEAR(nolm_[0].get_psi()[4], 1.210103935461e+00, max_tol);
    EXPECT_NEAR(nolm_[0].get_psi()[699], 4.465396560257e-08, max_tol);
    EXPECT_NEAR(nolm_[0].get_psi()[700], 0.0, max_tol);

    EXPECT_EQ(nolm_[1].getLabel(), "O");
    EXPECT_EQ(nolm_[1].getL(), 0);
    EXPECT_EQ(nolm_[1].getChi(), 1);
    EXPECT_EQ(nolm_[1].getNr(), 701);
    EXPECT_DOUBLE_EQ(nolm_[1].getRcut(), 7.0);
    EXPECT_NEAR(nolm_[1].get_psi()[0], 7.254873428942e-01, max_tol);
    EXPECT_NEAR(nolm_[1].get_psi()[1], 7.256666701836e-01, max_tol);
    EXPECT_NEAR(nolm_[1].get_psi()[4], 7.283448557011e-01, max_tol);
    EXPECT_NEAR(nolm_[1].get_psi()[699], -1.916246212603e-06, max_tol);
    EXPECT_NEAR(nolm_[1].get_psi()[700], 0.0, max_tol);

    EXPECT_EQ(nolm_[2].getLabel(), "O");
    EXPECT_EQ(nolm_[2].getL(), 1);
    EXPECT_EQ(nolm_[2].getChi(), 0);
    EXPECT_EQ(nolm_[2].getNr(), 701);
    EXPECT_DOUBLE_EQ(nolm_[2].getRcut(), 7.0);
    EXPECT_NEAR(nolm_[2].get_psi()[0], 0.0, max_tol);
    EXPECT_NEAR(nolm_[2].get_psi()[1], 4.626669306440e-02, max_tol);
    EXPECT_NEAR(nolm_[2].get_psi()[4], 1.845014292772e-01, max_tol);
    EXPECT_NEAR(nolm_[2].get_psi()[699], 2.870401658966e-07, max_tol);
    EXPECT_NEAR(nolm_[2].get_psi()[700], 0.0, max_tol);

    EXPECT_EQ(nolm_[3].getLabel(), "O");
    EXPECT_EQ(nolm_[3].getL(), 1);
    EXPECT_EQ(nolm_[3].getChi(), 1);
    EXPECT_EQ(nolm_[3].getNr(), 701);
    EXPECT_DOUBLE_EQ(nolm_[3].getRcut(), 7.0);
    EXPECT_NEAR(nolm_[3].get_psi()[0], 0.0, max_tol);
    EXPECT_NEAR(nolm_[3].get_psi()[1], 3.375340101333e-02, max_tol);
    EXPECT_NEAR(nolm_[3].get_psi()[4], 1.346256082234e-01, max_tol);
    EXPECT_NEAR(nolm_[3].get_psi()[699], -2.771091616120e-06, max_tol);
    EXPECT_NEAR(nolm_[3].get_psi()[700], 0.0, max_tol);

    EXPECT_EQ(nolm_[4].getLabel(), "O");
    EXPECT_EQ(nolm_[4].getL(), 2);
    EXPECT_EQ(nolm_[4].getChi(), 0);
    EXPECT_EQ(nolm_[4].getNr(), 701);
    EXPECT_DOUBLE_EQ(nolm_[4].getRcut(), 7.0);
    EXPECT_NEAR(nolm_[4].get_psi()[0], 0.0, max_tol);
    EXPECT_NEAR(nolm_[4].get_psi()[1], -3.343626342662e-04, max_tol);
    EXPECT_NEAR(nolm_[4].get_psi()[4], -5.337546547975e-03, max_tol);
    EXPECT_NEAR(nolm_[4].get_psi()[699], 1.396308876444e-06, max_tol);
    EXPECT_NEAR(nolm_[4].get_psi()[700], 0.0, max_tol);
}


TEST_F(NumericalOrbitalLmTest, Getters) {

    // this test checks whether all the getters work as expected
    // whether the values make sense or not is tested in "initialize"

    this->init();

    for (size_t i = 0; i != chi_.size(); ++i) {
        EXPECT_EQ(nolm_[i].getLabel(), nolm_[i].getLabel());
        EXPECT_EQ(nolm_[i].getType(), nolm_[i].getType());
        EXPECT_EQ(nolm_[i].getL(), nolm_[i].getL());
        EXPECT_EQ(nolm_[i].getChi(), nolm_[i].getChi());

        EXPECT_DOUBLE_EQ(nolm_[i].getDk(), nolm_[i].getDk());
        EXPECT_DOUBLE_EQ(nolm_[i].getDruniform(), dr_uniform_);
        EXPECT_EQ(nolm_[i].getPsiuniform(), &nolm_[i].psi_uniform[0]);
        EXPECT_EQ(nolm_[i].getDpsiuniform(), &nolm_[i].dpsi_uniform[0]);
        EXPECT_EQ(nolm_[i].getNruniform(), nolm_[i].nr_uniform);
        EXPECT_EQ(nolm_[i].getDruniform(), nolm_[i].dr_uniform);

        EXPECT_EQ(nolm_[i].getNr(), nolm_[i].getNr());
        EXPECT_EQ(nolm_[i].getNk(), nolm_[i].getNk());

        EXPECT_EQ(nolm_[i].getRcut(), nolm_[i].getRcut());
        EXPECT_EQ(nolm_[i].getKcut(), nolm_[i].getKcut());

        EXPECT_EQ(nolm_[i].getRadial(), &nolm_[i].get_r_radial()[0]);
        EXPECT_EQ(nolm_[i].get_r_radial(), nolm_[i].get_r_radial());

        EXPECT_EQ(nolm_[i].getRab(), &nolm_[i].get_rab()[0]);
        EXPECT_EQ(nolm_[i].get_rab(), nolm_[i].get_rab());

        EXPECT_EQ(nolm_[i].getDk(), nolm_[i].getDk());
        EXPECT_EQ(nolm_[i].getKpoint(), &nolm_[i].get_k_radial()[0]);
        EXPECT_EQ(nolm_[i].get_k_radial(), nolm_[i].get_k_radial());

        EXPECT_EQ(nolm_[i].getPsi(), &nolm_[i].get_psi()[0]);
        EXPECT_EQ(nolm_[i].getPsi_r(), &nolm_[i].get_psir()[0]);
        EXPECT_EQ(nolm_[i].getPsif(), &nolm_[i].get_psif()[0]);
        EXPECT_EQ(nolm_[i].getPsi_k(), &nolm_[i].get_psi_k()[0]);
        EXPECT_EQ(nolm_[i].getPsi_k2(), &nolm_[i].get_psi_k2()[0]);

        EXPECT_EQ(nolm_[i].get_psi(), nolm_[i].get_psi());
        EXPECT_EQ(nolm_[i].get_psif(), nolm_[i].get_psif());
        EXPECT_EQ(nolm_[i].get_psi_k(), nolm_[i].get_psi_k());
        EXPECT_EQ(nolm_[i].get_psi_k2(), nolm_[i].get_psi_k2());

        for (size_t ir = 0; ir != nolm_[i].get_r_radial().size(); ++ir) {
            EXPECT_EQ(nolm_[i].getRadial(ir), nolm_[i].get_r_radial()[ir]);
            EXPECT_EQ(nolm_[i].getRab(ir), nolm_[i].get_rab()[ir]);
            EXPECT_EQ(nolm_[i].getPsi(ir), nolm_[i].get_psi()[ir]);
            EXPECT_EQ(nolm_[i].getPsi_r(ir), nolm_[i].get_psir()[ir]);
        }

        for (size_t ik = 0; ik != nolm_[i].get_k_radial().size(); ++ik) {
            EXPECT_EQ(nolm_[i].getKpoint(ik), nolm_[i].get_k_radial()[ik]);
            EXPECT_EQ(nolm_[i].getPsif(ik), nolm_[i].get_psif()[ik]);
            EXPECT_EQ(nolm_[i].getPsi_k(ik), nolm_[i].get_psi_k()[ik]);
            EXPECT_EQ(nolm_[i].getPsi_k2(ik), nolm_[i].get_psi_k2()[ik]);
        }
    }
}


TEST_F(NumericalOrbitalLmTest, PsiNormalization) {

    // This test checks the normalization of
    // 1. the radial function in the real space
    // 2. the radial function in the k space
    // 3. the interpolated radial function (psi_uniform) in the real space

    // default ecutwfc might be inadequate, use a larger one instead
    this->init_with_different_k(1600.0, dk_);

    double radint = 0.0;
    double* rintegrand = new double[nr_];
    double* kintegrand = new double[nk_];

    for (size_t i = 0; i != chi_.size(); ++i) {

        // normalization check of chi(r)
        for (int ir = 0; ir != nr_; ++ir) {
            rintegrand[ir] = std::pow(nolm_[i].getPsi_r(ir), 2);
        }

        ModuleBase::Integral::Simpson_Integral(nr_, rintegrand, rab_, radint);
        EXPECT_NEAR(radint, 1.0, 1e-10);

        // normalization check of chi(k)
        for (size_t ik = 0; ik != nk_; ++ik) {
            kintegrand[ik] = std::pow(nolm_[i].getPsi_k(ik), 2);
        }

        ModuleBase::Integral::Simpson_Integral(nk_, kintegrand, dk_, radint);
        EXPECT_NEAR(radint, 1.0, 1e-6); // what is a reasonable error?

    }

    delete[] rintegrand;
    delete[] kintegrand;


    // check the normalization of psi_uniform
    for (size_t i = 0; i != chi_.size(); ++i) {

        int nr = nolm_[i].nr_uniform;

        // note that Simpson_Integral only accepts an odd number of mesh points
        if (nr%2 == 0) {
            ++nr;
        }

        rintegrand = new double[nr];
        for (int ir = 0; ir != nr; ++ir) {
            rintegrand[ir] = 0.0;
        }

        // normalization check of psi_uniform
        for (int ir = 0; ir != nolm_[i].nr_uniform; ++ir) {
            rintegrand[ir] = std::pow(ir*dr_uniform_*nolm_[i].psi_uniform[ir], 2);
        }

        ModuleBase::Integral::Simpson_Integral(nr, rintegrand, dr_uniform_, radint);
        EXPECT_NEAR(radint, 1.0, 1e-6);

        delete[] rintegrand;
    }
}


TEST_F(NumericalOrbitalLmTest, K2RConsistency) {

    // This test checks whether the results of
    // cal_kradial & cal_kradial_sbpool agree.

    this->init();

    double* chik_ = new double[nk_];
    double* kchik_ = new double[nk_];
    double* k2chik_ = new double[nk_];

    for (size_t i = 0; i != chi_.size(); ++i) {
        // save previous chi(k) results obtained by init()
        for (size_t ik = 0; ik != nk_; ++ik) {
            chik_[ik] = nolm_[i].getPsif(ik);
            kchik_[ik] = nolm_[i].getPsi_k(ik);
            k2chik_[ik] = nolm_[i].getPsi_k2(ik);
        }

        // use a different method than which used in init()
        if (flag_sbpool_) {
            cal_kradial(i);
        } else {
            cal_kradial_sbpool(i);
        }

        double max_tol = 1e-6;
        for (size_t ik = 0; ik != nk_; ++ik) {
            EXPECT_NEAR(chik_[ik], nolm_[i].getPsif(ik), max_tol);
            EXPECT_NEAR(kchik_[ik], nolm_[i].getPsi_k(ik), max_tol);
            EXPECT_NEAR(k2chik_[ik], nolm_[i].getPsi_k2(ik), max_tol);
        }
    }
}


TEST_F(NumericalOrbitalLmTest, R2K2RConsistency) {

    // This test checks whether cal_rradial_sbpool brings chi(k)
    // back to the original chi(r) by looking at the error
    //      \int dr r^2 (chi_in(r)-chi_out(r))^2

    this->init_with_different_k(1600.0, dk_);

    // original r*psi(r)
    double* rchi_ = new double[nr_];
    double* err_integrand = new double[nr_];
    double err = 0.0;

    // maximum tolerance of err
    double max_tol = 1e-6;

    for (size_t i = 0; i != chi_.size(); ++i) {
        for (int ir = 0; ir != nr_; ++ir) {
            rchi_[ir] = nolm_[i].getPsi_r(ir);
        }

        cal_rradial_sbpool(i);

        for (int ir = 0; ir != nr_; ++ir) {
            err_integrand[ir] = std::pow(rchi_[ir]-nolm_[i].getPsi_r(ir), 2);
        }

        ModuleBase::Integral::Simpson_Integral(nr_, err_integrand, rab_, err);
        EXPECT_LT(err, max_tol);
    }
}


TEST_F(NumericalOrbitalLmTest, FiniteDiffPsiUniform) {

    // this test checks whether dpsi_uniform agrees with the
    // finite difference of psi_uniform

    this->init();

    double max_tol = 1e-3;

    for (size_t i = 0; i != nolm_.size(); ++i) {
        std::vector<double>& f = nolm_[i].psi_uniform;
        for (int ir = 4; ir != nolm_[i].nr_uniform-4; ++ir) {
            double fd =
            ( +1.0/280*f[ir-4] - 4.0/105*f[ir-3] + 1.0/5*f[ir-2] - 4.0/5*f[ir-1]
              -1.0/280*f[ir+4] + 4.0/105*f[ir+3] - 1.0/5*f[ir+2] + 4.0/5*f[ir+1]
            ) / nolm_[i].dr_uniform;
            EXPECT_NEAR(fd, nolm_[i].dpsi_uniform[ir], max_tol);
        }
    }

}


TEST_F(NumericalOrbitalLmTest, PsiSave) {

    // This test checks whether plot() works as expected.
    //
    // set flag_plot to true will have set_orbital_info call plot(),
    // which saves psi, psik, psi_uniform & dpsi_uniform to files.

    // the directory that holds the data files is not created yet,
    // so files cannot be opened and plot() should fail at this stage.
    // but even if file cannot be opened, plot should not throw!
    flag_plot_ = true;
    ASSERT_NO_THROW(this->init());

    std::vector<std::string> orb{"s", "s", "p", "p", "d"};
    std::ifstream ifs;
    std::string tmp_r, tmp_k, tmp_beta, tmp;
    double tol = 1e-5;

    // should be PARAM.sys.global_out_dir+label
    // but in this unit test global_out_dir is empty string
    // see Numerical_Nonlocal_Lm::plot() for details
    std::string dir = "./O/";

    // we now creat the directory to hold data files
    mkdir(dir.c_str(), 0777);

    for (size_t i = 0; i != nolm_.size(); ++i) {

        // this call should successfully write data to files
        ASSERT_NO_THROW(plot(i));

        auto get_fname = [&] (std::string const& suffix) -> std::string {
            return dir+"/O-" + orb[i] + std::to_string(nolm_[i].getChi()+1)
                + "-orbital-" + suffix + ".dat";
        };

        std::string psi_fname = get_fname("r");
        std::string psik_fname = get_fname("k");
        std::string psiru_fname = get_fname("ru");
        std::string psidru_fname = get_fname("dru");

        EXPECT_TRUE(this->check_file_match(nolm_[i].getNr(),
                    nolm_[i].getRadial(), nolm_[i].getPsi(), tol, psi_fname));
        EXPECT_TRUE(this->check_file_match(nolm_[i].getNk(),
                    nolm_[i].getKpoint(), nolm_[i].getPsi_k(), tol, psik_fname));

        double* ru_mesh = new double[nolm_[i].nr_uniform];
        for (int ir = 0; ir != nolm_[i].nr_uniform; ++ir) {
            ru_mesh[ir] = ir*nolm_[i].dr_uniform;
        }

        EXPECT_TRUE(this->check_file_match(nolm_[i].nr_uniform,
                    ru_mesh, nolm_[i].getPsiuniform(), tol, psiru_fname));
        EXPECT_TRUE(this->check_file_match(nolm_[i].nr_uniform,
                    ru_mesh, nolm_[i].getDpsiuniform(), tol, psidru_fname));

        remove(psi_fname.c_str());
        remove(psik_fname.c_str());
        remove(psiru_fname.c_str());
        remove(psidru_fname.c_str());
    }

    remove(dir.c_str());

}


TEST_F(NumericalOrbitalLmTest, VariousPsiType) {

    // this test checks the behavior of set_orbital_info
    // under various input Psi_Type

    this->init_with_different_k(1600.0, dk_);

    // this tolerance is used for element-wise comparison
    // rather than an integration (as in r2k2r_consistency)
    double max_tol = 1e-3;

    for (size_t i = 0; i != nolm_.size(); ++i) {
        std::vector<double> psi_ref, psif_ref, psik_ref, psik2_ref;

        psi_ref = nolm_[i].get_psi();
        psif_ref = nolm_[i].get_psif();
        psik_ref = nolm_[i].get_psi_k();
        psik2_ref = nolm_[i].get_psi_k2();


        // alternative Psi_Type input

        // Psi_Type == Psif
        nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psif,
                &psif_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, flag_sbpool_, force_flag_);

        for (int ir = 0; ir != nolm_[i].getNr(); ++ir) {
            EXPECT_NEAR(nolm_[i].get_psi()[ir], psi_ref[ir], max_tol);
        }

        // Psi_Type == Psik
        nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psik,
                &psik_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, flag_sbpool_, force_flag_);

        for (int ir = 0; ir != nolm_[i].getNr(); ++ir) {
            EXPECT_NEAR(nolm_[i].get_psi()[ir], psi_ref[ir], max_tol);
        }

        // Psi_Type == Psik2
        nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psik2,
                &psik2_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, flag_sbpool_, force_flag_);

        for (int ir = 0; ir != nolm_[i].getNr(); ++ir) {
            EXPECT_NEAR(nolm_[i].get_psi()[ir], psi_ref[ir], max_tol);
        }

    }
}


TEST_F(NumericalOrbitalLmTest, TurnOffSphBesPool) {

    // checks the behavior of set_orbital_info when sbpool is turned off
    //
    // expected behavior:
    // if Psi_Type == Psi, set_orbital_info should work
    // if Psi_Type is any Fourier space type, set_orbital_info should throw

    flag_sbpool_ = false;
    EXPECT_NO_THROW(this->init());

    for (size_t i = 0; i != nolm_.size(); ++i) {
        std::vector<double> psi_ref, psif_ref, psik_ref, psik2_ref;

        psi_ref = nolm_[i].get_psi();
        psif_ref = nolm_[i].get_psif();
        psik_ref = nolm_[i].get_psi_k();
        psik2_ref = nolm_[i].get_psi_k2();

        EXPECT_NO_THROW(nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psi,
                &psi_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, false, force_flag_));
        EXPECT_THROW(nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psif,
                &psif_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, false, force_flag_), std::domain_error);
        EXPECT_THROW(nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psik,
                &psik_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, false, force_flag_), std::domain_error);
        EXPECT_THROW(nolm_[i].set_orbital_info(elem_label_, index_atom_type_,
                l_[i], index_chi_l_[i], nr_, rab_,
                r_radial_, Numerical_Orbital_Lm::Psi_Type::Psik2,
                &psik2_ref[0], nk_, dk_,
                dr_uniform_, flag_plot_, false, force_flag_), std::domain_error);
    }
}


int main(int argc, char **argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}



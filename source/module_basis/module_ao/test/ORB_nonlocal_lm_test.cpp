#include "gtest/gtest.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/math_integral.h"
#include <fstream>
#include <iomanip>

#define private public
#include "module_basis/module_ao/ORB_nonlocal_lm.h"
#undef private


#ifdef __MPI
#include <mpi.h>
#endif


/***********************************************************
 *      unit test of class "Numerical_Nonlocal_Lm"
 ***********************************************************/

/**
 * Tested functions
 *
 * - set_NL_proj
 *   copies the input into class members
 *
 * - get_kradial
 *   applies a radial Fourier transform to beta_r to obtain beta_k
 *
 * - freemem
 *   deallocates the allocated memory of r_radial, rab, beta_r, 
 *   beta_uniform, dbeta_uniform, k_radial & beta_k, and set them 
 *   to nullptr
 *
 * - renew
 *   allocates memory for r_radial, rab, beta_r, beta_uniform,
 *   dbeta_uniform, k_radial & beta_k, and set their content to all zero.
 *
 * - plot
 *   saves beta_r, beta_k & beta_uniform to file
 *
 * - all "getters"
 *   get access to class members
 *
 * - operator=
 *   enables deep copy
 */

class NumericalNonlocalLmTest : public ::testing::Test
{
protected:

    void SetUp();
    void TearDown();

    // a vector of objects under unit test
    // corresponds to all the nonlocal projectors in the above upf file
    std::vector<Numerical_Nonlocal_Lm> nnl;

    // helper functions
    void init();
    std::string trim(std::string const&);
    double err_r2k2r(Numerical_Nonlocal_Lm&);
    size_t calc_nk(double const& ecutwfc, double const& dk);
    void change_k(Numerical_Nonlocal_Lm&, double const& ecut, double const& dk);
    bool check_file_match(size_t const& nline, double const* col1, double const* col2, double const& tol, std::string const& fname);

    // number of beta projectors
    size_t nproj_;

    // data used to initialized Numerical_Nonlocal_Lm objects
    std::string elem_label_;
    int index_atom_type_;
    double* rab_;
    double* r_radial_;
    std::vector<int> l_;
    std::vector<int> nr_; // number of non-zero points of beta
    std::vector<double*> beta_r_;
    size_t nk_;
    double dk_;
    double dr_uniform_;
};


void NumericalNonlocalLmTest::SetUp() {

    ///////////////////////////////////////////////////
    //                  Parameters
    ///////////////////////////////////////////////////

    // upf file to read data from
    std::string upf_file = "./lcao_H2O/O_ONCV_PBE-1.0.upf";

    // energy cutoff for 1D integration (radial Fourier transform)
    // this number is used in most tests except r2k2r_consistency
    // where a much larger value is used.
    double ecutwfc = 100.0;

    // In normal ABACUS calculation, dk, nk & dr_uniform 
    // are retrieved from the LCAO_Orbitals object, 
    // see setupNonlocal() in module_cell/setup_nonlocal.cpp.
    // Here we just provide some reasonable values.
    dk_ = 0.01;
    dr_uniform_ = 0.001;

    nk_ = calc_nk(ecutwfc, dk_);

    // not really meaningful in this unit test
    index_atom_type_ = 42;

    ///////////////////////////////////////////////////
    //                  Read UPF file
    ///////////////////////////////////////////////////

    // variables below will be read from upf_file
    size_t nmesh_upf = 0;

    std::ifstream ifs(upf_file);

    // see read_pseudo_upf201 in module_cell/read_pp_upf201.cpp
    // the following code only works for UPF files of version 2.0.1

    //----------- read header ------------
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_HEADER");

    // each line looks like 
    //    label="   value"
    // (note the whitespaces at the beginning & within the quotation marks)
    // the last line has a "/>" marking the end of PP_HEADER
    std::string linebuf;
    std::string label;
    std::string val;

    while(std::getline(ifs, linebuf)) {

        std::string::size_type pos = linebuf.find('=');
        if (pos != std::string::npos) {
            // extract the label
            // skip leading whitespaces
            auto label_start = linebuf.find_first_not_of(" \t");
            label = linebuf.substr(label_start, pos-label_start);
        } else if (linebuf.find("/>") != std::string::npos) {
            // reach the end of PP_HEADER
            break;
        } else {
            // skip lines without '=' or "/>"
            continue;
        }

        val = trim(linebuf);

        // only the entries below are relevant to the current unit test
        if (label == "element") {
            elem_label_ = val;
        } else if (label == "number_of_proj") {
            nproj_ = std::atoi(val.c_str());
        } else if (label == "mesh_size") {
            nmesh_upf = std::atoi(val.c_str());
        }

        if (linebuf.find("/>") != std::string::npos) {
            // reach the end of PP_HEADER
            break;
        }
    }

    //----------- read mesh ------------
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R");
    // skip the rest of the line whose info is unimportant or duplicate
    std::getline(ifs, linebuf);

    r_radial_ = new double[nmesh_upf];
    for (size_t ir = 0; ir != nmesh_upf; ++ir) {
        ifs >> r_radial_[ir];
    }

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB");
    // skip the rest of the line whose info is unimportant or duplicate
    std::getline(ifs, linebuf);

    rab_ = new double[nmesh_upf];
    for (size_t ir = 0; ir != nmesh_upf; ++ir) {
        ifs >> rab_[ir];
    }
    
    //----------- read nonlocal projectors & initialize objects ------------

    l_.resize(nproj_);
    nr_.resize(nproj_);
    beta_r_.resize(nproj_);
    nnl.resize(nproj_);

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
    for (size_t iproj = 0; iproj != nproj_; ++iproj) {

        beta_r_[iproj] = new double[nmesh_upf];

        // process beta headers
        while(std::getline(ifs, linebuf)) {
            if (linebuf.find('>') != std::string::npos) {
                // end of PP_BETA header, beta mesh starts
                break;
            }

            if (linebuf.find("angular_momentum") != std::string::npos) {
                l_[iproj] = std::atoi(trim(linebuf).c_str());
            }
        }

        // read beta mesh (UPF <PP_BETA> mesh is r*beta(r)!)
        for (size_t ir = 0; ir != nmesh_upf; ++ir) {
            ifs >> beta_r_[iproj][ir];
        }

        // determine the actual mesh size by ignoring the trailing 
        // zeros (or small numbers) of beta

        // the code could fail if the original upf mesh size is even
        // and there's no trailing zero. Currently we assume that 
        // there's always at least one trailing zero.
        for (nr_[iproj] = nmesh_upf; nr_[iproj] > 1; --nr_[iproj]) {
            if (std::abs(beta_r_[iproj][nr_[iproj]-1])>1e-14) {
                break;
            }
        }

        if (nr_[iproj]%2 == 0) {
            ++nr_[iproj];
        }

        while(std::getline(ifs, linebuf)) {
            if (linebuf.find("/PP_BETA") != std::string::npos) {
                // reach beta mesh ending symbol </PP_BETA.x>
                break;
            }
        }

    }
}


void NumericalNonlocalLmTest::TearDown() {
    delete[] r_radial_;
    delete[] rab_;
    for (size_t ip = 0; ip != beta_r_.size(); ++ip) {
        delete[] beta_r_[ip];
    }
}


void NumericalNonlocalLmTest::init() {
    // initialized the tested objects by pouring the data
    // collected in SetUp() to nnl;
    for (size_t iproj = 0; iproj != nproj_; ++iproj) {
        nnl[iproj].set_NL_proj(elem_label_, index_atom_type_, l_[iproj], 
                nr_[iproj], rab_, r_radial_, beta_r_[iproj], 
                nk_, dk_, dr_uniform_);

        /*
        // normalization check 
        double* tmp = new double[nr_[iproj]];
        for (int ir = 0; ir != nr_[iproj]; ++ir) {
            tmp[ir] = beta_r_[iproj][ir]*beta_r_[iproj][ir];
        }
        double radint = 0.0;
        ModuleBase::Integral::Simpson_Integral(nr_[iproj], tmp, rab_, radint);
        std::cout << "proj " << iproj 
            << std::fixed << std::setprecision(12)
            << "    radial integral = " << radint << std::endl;
        */
    }
}


std::string NumericalNonlocalLmTest::trim(std::string const& str) {

    // extract the substring between quotation marks (with whitespace trimmed)
    // str MUST contain a pair of quotation marks

    auto start = str.find('"');
    auto end = str.find_last_of('"');
    std::string tmp = str.substr(start+1, end-start-1);

    if (tmp.length() == 0) {
        return tmp;
    }

    start = tmp.find_first_not_of(" \t");
    end = tmp.find_last_not_of(" \t");
    return tmp.substr(start, end-start+1);
}


double NumericalNonlocalLmTest::err_r2k2r(Numerical_Nonlocal_Lm& nnl_tmp) {

    // err_r2k2r makes use of Numerical_Nonlocal_Lm::get_kradial()
    // to transform beta_k back to beta_r. 
    // The error is computed as the difference between the origional 
    // and the transformed beta_r.

    double* beta_r_old = new double[nnl_tmp.nr];
    double* rab_old = new double[nnl_tmp.nr];
    double* err = new double[nnl_tmp.nr];

    for (int ir = 0; ir != nnl_tmp.nr; ++ir) {
        beta_r_old[ir] = nnl_tmp.beta_r[ir];
        rab_old[ir] = nnl_tmp.rab[ir];
    }

    std::swap(nnl_tmp.nr, nnl_tmp.nk);
    std::swap(nnl_tmp.r_radial, nnl_tmp.k_radial);
    std::swap(nnl_tmp.beta_r, nnl_tmp.beta_k);

    delete[] nnl_tmp.rab;
    nnl_tmp.rab = new double[nnl_tmp.nr];
    for (int ik = 0; ik != nnl_tmp.nr; ++ik) {
        nnl_tmp.rab[ik] = nnl_tmp.dk;
    }

    nnl_tmp.get_kradial();

    for (int ir = 0; ir != nnl_tmp.nk; ++ir) {
        err[ir] = std::pow(nnl_tmp.getBeta_k(ir) - beta_r_old[ir], 2);
    }

    double errint = 0.0;
    ModuleBase::Integral::Simpson_Integral(nnl_tmp.nk, err, rab_old, errint);

    return errint;
}


size_t NumericalNonlocalLmTest::calc_nk(double const& ecutwfc, double const& dk) {

    // current formula for nk
    // see module_basis/module_ao/ORB_read.cpp, function "Read_Orbitals"

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


void NumericalNonlocalLmTest::change_k(Numerical_Nonlocal_Lm& nnl_, double const& ecut, double const& dk) {

    // recalculates k mesh & beta_k with given ecut & dk
    // used in r2k2r_consistency

    Numerical_Nonlocal_Lm tmp;
    tmp = nnl_;

    nnl_.set_NL_proj(
            tmp.label,
            tmp.index_atom_type,
            tmp.angular_momentum_l,
            tmp.nr,
            tmp.rab,
            tmp.r_radial,
            tmp.beta_r,
            this->calc_nk(ecut, dk),
            dk,
            tmp.dr_uniform
            );
}


TEST_F(NumericalNonlocalLmTest, Init) {
    
    this->init();

    for (size_t ip = 0; ip != nproj_; ++ip) {
        EXPECT_EQ(elem_label_, nnl[ip].label);
        EXPECT_EQ(index_atom_type_, nnl[ip].index_atom_type);
        EXPECT_EQ(l_[ip], nnl[ip].angular_momentum_l);
        EXPECT_EQ(dr_uniform_, nnl[ip].dr_uniform);
        EXPECT_EQ(nr_[ip], nnl[ip].nr);
        EXPECT_EQ(r_radial_[nr_[ip]-1], nnl[ip].rcut);
        EXPECT_EQ(nk_, nnl[ip].nk);
        EXPECT_EQ(dk_, nnl[ip].dk);
        
        // freemem() & renew() will be tested elsewhere

        for (int ir = 0; ir != nr_[ip]; ++ir) {
            EXPECT_EQ(r_radial_[ir], nnl[ip].r_radial[ir]);
            EXPECT_EQ(rab_[ir], nnl[ip].rab[ir]);
            EXPECT_EQ(beta_r_[ip][ir], nnl[ip].beta_r[ir]);
        }

        for (size_t ik = 0; ik != nk_; ++ik) {
            EXPECT_EQ(ik*dk_, nnl[ip].k_radial[ik]);
        }
        EXPECT_EQ((nk_-1)*dk_, nnl[ip].kcut);

        // get_kradial() will be tested elsewhere
    }
}


TEST_F(NumericalNonlocalLmTest, Getters) {

    // this test checks all the getter functions

    this->init();

    for (size_t iproj = 0; iproj != nnl.size(); ++iproj) {
        EXPECT_EQ(nnl[iproj].getType(), 42); // index_atom_type
        EXPECT_DOUBLE_EQ(nnl[iproj].getDk(), 0.01);

        ASSERT_NE(nnl[iproj].getRadial(), nullptr);
        ASSERT_NE(nnl[iproj].getKpoint(), nullptr);
        ASSERT_NE(nnl[iproj].getBeta_r(), nullptr);
        ASSERT_NE(nnl[iproj].getBeta_k(), nullptr);

        for (int ir = 0; ir != nnl[iproj].nr; ++ir) {
            EXPECT_DOUBLE_EQ(nnl[iproj].getRadial(ir), 0.01*ir);
            EXPECT_DOUBLE_EQ(nnl[iproj].getRadial()[ir], 0.01*ir);
            EXPECT_DOUBLE_EQ(nnl[iproj].getBeta_r()[ir], nnl[iproj].getBeta_r(ir));
        }

        for (int ik = 0; ik != nnl[iproj].nk; ++ik) {
            EXPECT_DOUBLE_EQ(nnl[iproj].getKpoint(ik), ik*0.01);
            EXPECT_DOUBLE_EQ(nnl[iproj].getKpoint()[ik], ik*0.01);
            EXPECT_DOUBLE_EQ(nnl[iproj].getBeta_k()[ik], nnl[iproj].getBeta_k(ik));
        }
    }

    // see upf file for details
    EXPECT_EQ(nnl[0].getL(), 0);
    EXPECT_DOUBLE_EQ(nnl[0].getRcut(), 1.3200);
    EXPECT_DOUBLE_EQ(nnl[0].getBeta_r(0), 0.0);
    EXPECT_DOUBLE_EQ(nnl[0].getBeta_r(1), -8.2277987587e-02);
    EXPECT_DOUBLE_EQ(nnl[0].getBeta_r(4), -3.2850076507e-01);
    EXPECT_DOUBLE_EQ(nnl[0].getBeta_r(132), 1.6220072646e-06);


    EXPECT_EQ(nnl[1].getL(), 0);
    EXPECT_DOUBLE_EQ(nnl[1].getRcut(), 1.3200);
    EXPECT_DOUBLE_EQ(nnl[1].getBeta_r(0), 0.0);
    EXPECT_DOUBLE_EQ(nnl[1].getBeta_r(1), -1.1723087215e-02);
    EXPECT_DOUBLE_EQ(nnl[1].getBeta_r(4), -4.2201369170e-02);
    EXPECT_DOUBLE_EQ(nnl[1].getBeta_r(132), -2.7996494097e-06);

    EXPECT_EQ(nnl[2].getL(), 1);
    EXPECT_DOUBLE_EQ(nnl[2].getRcut(), 1.5000);
    EXPECT_DOUBLE_EQ(nnl[2].getBeta_r(0), 0.0);
    EXPECT_DOUBLE_EQ(nnl[2].getBeta_r(1), 3.5860269827e-03);
    EXPECT_DOUBLE_EQ(nnl[2].getBeta_r(4), 5.6837367274e-02);
    EXPECT_DOUBLE_EQ(nnl[2].getBeta_r(150), 9.6147639048e-06);

    EXPECT_EQ(nnl[3].getL(), 1);
    EXPECT_DOUBLE_EQ(nnl[3].getRcut(), 1.5000);
    EXPECT_DOUBLE_EQ(nnl[3].getBeta_r(0), 0.0);
    EXPECT_DOUBLE_EQ(nnl[3].getBeta_r(1), 9.2893242255e-04);
    EXPECT_DOUBLE_EQ(nnl[3].getBeta_r(4), 1.4588689808e-02);
    EXPECT_DOUBLE_EQ(nnl[3].getBeta_r(150), 5.9608986436e-06);
}


TEST_F(NumericalNonlocalLmTest, DeepCopy) {

    // this test checks whether the operator= overload properly
    // performs a deep copy

    this->init();

    Numerical_Nonlocal_Lm tmp;

    size_t iproj = 3;
    tmp = nnl[iproj];

    EXPECT_EQ(tmp.label, nnl[iproj].label);
    EXPECT_EQ(tmp.index_atom_type, nnl[iproj].index_atom_type);
    EXPECT_EQ(tmp.angular_momentum_l, nnl[iproj].angular_momentum_l);
    EXPECT_EQ(tmp.nr, nnl[iproj].nr);
    EXPECT_EQ(tmp.nk, nnl[iproj].nk);
    EXPECT_EQ(tmp.index_proj, nnl[iproj].index_proj);

    EXPECT_DOUBLE_EQ(tmp.rcut, nnl[iproj].rcut);
    EXPECT_DOUBLE_EQ(tmp.kcut, nnl[iproj].kcut);
    EXPECT_DOUBLE_EQ(tmp.dk, nnl[iproj].dk);

    ASSERT_NE(tmp.getRadial(), nullptr);
    ASSERT_NE(tmp.getKpoint(), nullptr);
    ASSERT_NE(tmp.getBeta_k(), nullptr);
    ASSERT_NE(tmp.getBeta_r(), nullptr);

    for (int ir = 0; ir != nnl[iproj].nr; ++ir) {
        EXPECT_DOUBLE_EQ(tmp.r_radial[ir], nnl[iproj].r_radial[ir]);
        EXPECT_DOUBLE_EQ(tmp.rab[ir], nnl[iproj].rab[ir]);
        EXPECT_DOUBLE_EQ(tmp.beta_r[ir], nnl[iproj].beta_r[ir]);
    }

    for (int ik = 0; ik != nnl[iproj].nk; ++ik) {
        EXPECT_DOUBLE_EQ(tmp.k_radial[ik], nnl[iproj].k_radial[ik]);
        EXPECT_DOUBLE_EQ(tmp.beta_k[ik], nnl[iproj].beta_k[ik]);
    }
}


TEST_F(NumericalNonlocalLmTest, R2K2RConsistency) {

    /*
     * This test checks whether get_kradial() works as expected 
     * by applying it to k*beta(k) and see if it gives back the 
     * original r*beta(r).
     *
     * This is NOT a convergence test, so a large ecutwfc is used.
     *
     * get_kradial() transforms r*beta(r) to k*beta(k) according to
     * 
     * beta(k) = sqrt(2/pi) * \int_0^{\infty} beta(r) * jl(kr) * r^2 dr
     *
     * whose inverse transform
     *
     * beta(r) = sqrt(2/pi) * \int_0^{\infty} beta(k) * jl(kr) * k^2 dk
     *
     * is just itself.
     *
     */

    this->init();

    double ecut = 1000.0;
    for (size_t iproj = 0; iproj != nnl.size(); ++iproj) {
        Numerical_Nonlocal_Lm tmp;
        tmp = nnl[iproj];
        this->change_k(tmp, ecut, tmp.dk);
        EXPECT_LT(err_r2k2r(tmp), 1e-6);
    }
}


/*
TEST_F(NumericalNonlocalLmTest, R2K2RConsistencyMany) {

    this->init();

    std::vector<double> ecut_list{100.0, 200.0, 400.0, 800.0, 1600.0};

    std::cout << "ecut convergence test" << std::endl;

    for (size_t iproj = 0; iproj != nnl.size(); ++iproj) {
        for (size_t ie = 0; ie != ecut_list.size(); ++ie) {
            Numerical_Nonlocal_Lm tmp;
            tmp = nnl[iproj];
            this->change_k(tmp, ecut_list[ie], tmp.dk);
            double err = err_r2k2r(tmp);
            std::cout << "proj = " << iproj
                << "    ecut = "  << std::setw(8) << ecut_list[ie]
                << "    dk = "    << std::setw(6) << tmp.dk
                << "    error = " << std::setw(10) << err 
                << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;


    std::vector<double> dk_list{0.001, 0.01, 0.1, 
        0.3, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0};

    std::cout << "dk convergence test" << std::endl;

    double ecut = 400.0;
    for (size_t iproj = 0; iproj != nnl.size(); ++iproj) {
        for (size_t idk = 0; idk != dk_list.size(); ++idk) {
            Numerical_Nonlocal_Lm tmp;
            tmp = nnl[iproj];
            this->change_k(tmp, ecut, dk_list[idk]);
            double err = err_r2k2r(tmp);
            std::cout << "proj = " << iproj
                << "    ecut = "  << std::setw(6)  << ecut
                << "    dk = "    << std::setw(6)  << dk_list[idk]
                << "    error = " << std::setw(12) << err
                << "    nk = "    << std::setw(5)  << tmp.nr
                << std::endl;
        }
        std::cout << std::endl;
    }
}
*/


TEST_F(NumericalNonlocalLmTest, FreeAndRenew) {

    this->init();

    EXPECT_NE(nnl[0].r_radial, nullptr);
    EXPECT_NE(nnl[0].rab, nullptr);
    EXPECT_NE(nnl[0].beta_r, nullptr);
    EXPECT_NE(nnl[0].beta_uniform, nullptr);
    EXPECT_NE(nnl[0].dbeta_uniform, nullptr);
    EXPECT_NE(nnl[0].k_radial, nullptr);
    EXPECT_NE(nnl[0].beta_k, nullptr);

    nnl[0].freemem();

    EXPECT_EQ(nnl[0].r_radial, nullptr);
    EXPECT_EQ(nnl[0].rab, nullptr);
    EXPECT_EQ(nnl[0].beta_r, nullptr);
    EXPECT_EQ(nnl[0].beta_uniform, nullptr);
    EXPECT_EQ(nnl[0].dbeta_uniform, nullptr);
    EXPECT_EQ(nnl[0].k_radial, nullptr);
    EXPECT_EQ(nnl[0].beta_k, nullptr);

    nnl[0].renew();

    ASSERT_NE(nnl[0].r_radial, nullptr);
    ASSERT_NE(nnl[0].rab, nullptr);
    ASSERT_NE(nnl[0].beta_r, nullptr);
    ASSERT_NE(nnl[0].beta_uniform, nullptr);
    ASSERT_NE(nnl[0].dbeta_uniform, nullptr);
    ASSERT_NE(nnl[0].k_radial, nullptr);
    ASSERT_NE(nnl[0].beta_k, nullptr);

    for (int ir = 0; ir != nnl[0].nr; ++ir) {
        EXPECT_DOUBLE_EQ(nnl[0].r_radial[ir], 0.0);
        EXPECT_DOUBLE_EQ(nnl[0].rab[ir], 0.0);
        EXPECT_DOUBLE_EQ(nnl[0].beta_r[ir], 0.0);
    }

    for (int ir = 0; ir != nnl[0].nr_uniform; ++ir) {
        EXPECT_DOUBLE_EQ(nnl[0].beta_uniform[ir], 0.0);
        EXPECT_DOUBLE_EQ(nnl[0].dbeta_uniform[ir], 0.0);
    }
        
    for (int ik = 0; ik != nnl[0].nk; ++ik) {
        EXPECT_DOUBLE_EQ(nnl[0].k_radial[ik], 0.0);
        EXPECT_DOUBLE_EQ(nnl[0].beta_k[ik], 0.0);
    }
}


bool NumericalNonlocalLmTest::check_file_match(size_t const& nline, double const* col1, double const* col2, double const& tol, std::string const& fname) {

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


TEST_F(NumericalNonlocalLmTest, BetaSave) {

    this->init();

    std::vector<std::string> orb{"s", "s", "p", "p"};
    std::ifstream ifs;
    std::string tmp_r, tmp_k, tmp_beta, tmp;
    double tol = 1e-5;

    // should be GlobalV::global_out_dir+label
    // but in this unit test global_out_dir is empty string
    // see Numerical_Nonlocal_Lm::plot() for details
    std::string dir = "./O/";

    mkdir(dir.c_str(), 0777);

    for (size_t i : {0, 2}) {

        ASSERT_NO_THROW(nnl[i].plot(0));

        std::string betar_fname = dir+"/O-" + orb[i] + "-proj-r.dat";
        std::string betak_fname = dir+"/O-" + orb[i] + "-proj-k.dat";
        std::string betaru_fname = dir+"/O-" + orb[i] + "-proj-ru.dat";

        EXPECT_EQ(true, this->check_file_match(nnl[i].nr, 
                    nnl[i].r_radial, nnl[i].beta_r, tol, betar_fname));
        EXPECT_EQ(true, this->check_file_match(nnl[i].nk, 
                    nnl[i].k_radial, nnl[i].beta_k, tol, betak_fname));

        double* r_uniform_mesh = new double[nnl[i].nr_uniform];
        for (int ir = 0; ir != nnl[i].nr_uniform; ++ir) {
            r_uniform_mesh[ir] = ir*nnl[i].dr_uniform;
        }
        EXPECT_EQ(true, this->check_file_match(nnl[i].nr_uniform, 
                    r_uniform_mesh, nnl[i].beta_uniform, tol, betaru_fname));

        remove(betar_fname.c_str());
        remove(betak_fname.c_str());
        remove(betaru_fname.c_str());
    }

    remove(dir.c_str());

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






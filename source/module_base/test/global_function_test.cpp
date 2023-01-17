#include "../global_function.h"
#include "../global_variable.h"
#include "../vector3.h"
#include "../blas_connector.h"
#include "../tool_quit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <time.h>

/************************************************
 *  unit test of functions in global_function
 ***********************************************/

/**
 * - Tested Function
 * - NewPart
 *   - note the start of new calculation part
 * - OutV1
 *   - output a string name with format
 * - OutV2
 *   - output a string name and its value with format
 * - OutV3
 *   - output a string name and three value with format
 * - OutP
 *   - output a string paramter, its value and explanation with format
 * - MakeDir
 *   - make a directory
 * - OutTime
 *   - print out calculation time > 0.1 minites
 * - AutoSet
 *   - note the setting of values in code
 * - Done
 *   - output info. and time on screen and log
 * - Zero
 *   - zero out complex array
 * - Scan
 *   - SCAN_BEGIN and SCAN_END used to read xml files
 * - MapExist
 *   - search the existence of an index in a map
 * - ReadValue
 *   - read a value and delete "\n"
 * - Dcopy
 *   - copy vector
 * - VectorToPointer
 *   - get the first pointer of (const) vector or valarray
 * - Note
 *   - print out warning info in running.log file
 * - COPYARRAY
 *   - copy complex or double arrays
 */

inline void EXPECT_COMPLEX_EQ(const std::complex<double>& a, const std::complex<double>& b)
{
    EXPECT_DOUBLE_EQ(a.real(), b.real());
    EXPECT_DOUBLE_EQ(a.imag(), b.imag());
}

class GlobalFunctionTest : public testing::Test
{
  protected:
    std::ofstream ofs;
    std::ifstream ifs;
    time_t start, end;
    // for capturing output in files and on screen
    std::string output;
    void SetUp()
    {
        GlobalV::ofs_warning.open("warning.log");
        GlobalV::ofs_running.open("running.log");
    }
    void TearDown()
    {
        GlobalV::ofs_warning.close();
        GlobalV::ofs_running.close();
        remove("warning.log");
        remove("running.log");
        remove("tmp");
    }
};

TEST_F(GlobalFunctionTest, NewPart)
{
    ModuleBase::GlobalFunc::NEW_PART("New Part Starts ...");
    GlobalV::ofs_running.close();
    ifs.open("running.log");
    getline(ifs, output);
    getline(ifs, output);
    getline(ifs, output);
    getline(ifs, output);
    // output in running.log file
    EXPECT_THAT(output, testing::HasSubstr("New Part Starts ..."));
    ifs.close();
}

TEST_F(GlobalFunctionTest, OutScreen)
{
	testing::internal::CaptureStdout();
	int nbx = 100;
	double rcut = 10.5;
	ModuleBase::GlobalFunc::OUT("nbx", nbx);
	ModuleBase::GlobalFunc::OUT("rcut", rcut);
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nbx = 100"));
	EXPECT_THAT(output,testing::HasSubstr("rcut = 10.5"));
}

TEST_F(GlobalFunctionTest, OutV1)
{
    ofs.open("tmp");
    ModuleBase::GlobalFunc::OUT(ofs, "abacus");
    ofs.close();
    ifs.open("tmp");
    getline(ifs, output);
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("abacus"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, OutV2)
{
    ofs.open("tmp");
    double ecut = 400;
    long nkpt = 1000;
    unsigned long nbands = 1000;
    char para[5] = "abcd";
    char key[6] = "abcde";
    char alph[29] = "abcdefghijklmnopqrstuvwxyzxy";
    ModuleBase::GlobalFunc::OUT(ofs, "ecut", ecut);
    ModuleBase::GlobalFunc::OUT(ofs, "nkpt", nkpt);
    ModuleBase::GlobalFunc::OUT(ofs, "nbands", nbands);
    ModuleBase::GlobalFunc::OUT(ofs, "para", para);
    ModuleBase::GlobalFunc::OUT(ofs, "key", key);
    ModuleBase::GlobalFunc::OUT(ofs, "alph", alph);
    ofs.close();
    ifs.open("tmp");
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("ecut = 400"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("nkpt = 1000"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("nbands = 1000"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("para = abcd"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("key = abcde"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("alph = abcdefghijklmnopqrstuvwxyzxy"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, OutV3)
{
    ofs.open("tmp");
    int nx = 100;
    int ny = 125;
    int nz = 375;
    double ax = 1.1;
    double ay = 2.2;
    double az = 3.3;
    ModuleBase::GlobalFunc::OUT(ofs, "grid", nx, ny, nz);
    ModuleBase::GlobalFunc::OUT(ofs, "direct", ax, ay, az);
    ofs.close();
    ifs.open("tmp");
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("[grid] = 100, 125, 375"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("[direct] = 1.1, 2.2, 3.3"));
    ifs.close();
}
// P for parameters
TEST_F(GlobalFunctionTest, OutP)
{
    ofs.open("tmp");
    double ecut = 400;
    std::string explanation = "energy cutoff";
    ofs << std::setiosflags(std::ios::left);
    ModuleBase::GlobalFunc::OUTP(ofs, "ecut", ecut, explanation);
    ofs.close();
    ifs.open("tmp");
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("ecut                           400 #energy cutoff"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, MakeDir)
{
    GlobalV::MY_RANK = 0;
    ModuleBase::GlobalFunc::MAKE_DIR("scf");
    std::system("test -d ");
    std::system("rm -r scf ");
    SUCCEED();
}

TEST_F(GlobalFunctionTest, OutTime)
{
    std::string name = "scf";
    start = time(NULL);
    end = time(NULL) + 200;
    ModuleBase::GlobalFunc::OUT_TIME(name, start, end);
    GlobalV::ofs_warning.close();
    ifs.open("warning.log");
    getline(ifs, output);
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("NAME < scf >"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, AutoSet)
{
    bool fftwan = true;
    ModuleBase::GlobalFunc::AUTO_SET("fftwan", fftwan);
    int NBANDS = 100;
    ModuleBase::GlobalFunc::AUTO_SET("NBANDS", NBANDS);
    GlobalV::ofs_warning.close();
    ifs.open("warning.log");
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("AUTO_SET fftwan to 1"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("AUTO_SET NBANDS to 100"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, Done)
{
    ofs.open("tmp");
    testing::internal::CaptureStdout();
    ModuleBase::GlobalFunc::DONE(ofs, "SETUP UNITCELL");
    // output on screen
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("DONE"));
    EXPECT_THAT(output, testing::HasSubstr("SETUP UNITCELL"));
    ofs.close();
    // output in file
    ifs.open("tmp");
    std::string outputf;
    getline(ifs, outputf);
    EXPECT_THAT(outputf, testing::HasSubstr("DONE"));
    EXPECT_THAT(outputf, testing::HasSubstr("SETUP UNITCELL"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, Zero)
{
    // int size
    std::complex<double>* porter = nullptr;
    double* pt = nullptr;
    bool* poc = nullptr;
    ModuleBase::Vector3<double>* u = nullptr;
    int size = 1000;
    porter = new std::complex<double>[size];
    pt = new double[size];
    poc = new bool[size];
    u = new ModuleBase::Vector3<double>[size];
    for (int i = 0; i < size; ++i)
    {
	    u[i].set(1.1,2.2,3.3);
    }
    double vu = 4.8;
    double zo = 0.0;
    bool r = true;
    std::complex<double> value{1.1, 2.2};
    std::complex<double> zero{0.0, 0.0};
    std::fill(&porter[0], &porter[size], value);
    std::fill(&pt[0],&pt[size],vu);
    std::fill(&poc[0],&poc[size],r);
    ModuleBase::GlobalFunc::ZEROS(porter, size);
    ModuleBase::GlobalFunc::ZEROS(pt, size);
    ModuleBase::GlobalFunc::ZEROS(poc, size);
    ModuleBase::GlobalFunc::ZEROS(u, size);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_COMPLEX_EQ(porter[i], zero);
        EXPECT_DOUBLE_EQ(pt[i],zo);
        EXPECT_FALSE(poc[i]);
        EXPECT_DOUBLE_EQ(u[i].x,zo);
        EXPECT_DOUBLE_EQ(u[i].y,zo);
        EXPECT_DOUBLE_EQ(u[i].z,zo);
    }
    delete[] porter;
    delete[] pt;
    delete[] poc;
    delete[] u;
    // long size
    long size1= 100;
    porter = new std::complex<double>[size1];
    pt = new double[size1];
    std::fill(&porter[0], &porter[size1], value);
    std::fill(&pt[0],&pt[size1],vu);
    ModuleBase::GlobalFunc::ZEROS(porter, size1);
    ModuleBase::GlobalFunc::ZEROS(pt, size1);
    for (int i = 0; i < size1; ++i)
    {
        EXPECT_COMPLEX_EQ(porter[i], zero);
        EXPECT_DOUBLE_EQ(pt[i],zo);
    }
    delete[] porter;
    delete[] pt;
    // unsigned long size
    unsigned long size2= 500;
    porter = new std::complex<double>[size2];
    pt = new double[size2];
    std::fill(&porter[0], &porter[size2], value);
    std::fill(&pt[0],&pt[size2],vu);
    ModuleBase::GlobalFunc::ZEROS(porter, size2);
    ModuleBase::GlobalFunc::ZEROS(pt, size2);
    for (int i = 0; i < size2; ++i)
    {
        EXPECT_COMPLEX_EQ(porter[i], zero);
        EXPECT_DOUBLE_EQ(pt[i],zo);
    }
    delete[] porter;
    delete[] pt;
}

TEST_F(GlobalFunctionTest, Scan)
{
    ofs.open("tmp");
    ofs << "<PP_MESH>" << std::endl;
    ofs << "100 100 100" << std::endl;
    ofs << "</PP_MESH>" << std::endl;
    ofs.close();
    ifs.open("tmp");
    EXPECT_FALSE(ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP>"));
    getline(ifs, output);
    getline(ifs, output);
    // std::cout << output << std::endl;
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP>");
    EXPECT_TRUE(ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH>"));
    getline(ifs, output);
    getline(ifs, output);
    // std::cout << output << std::endl;
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_MESH>");
    ifs.close();
    ifs.open("warning.log");
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("In SCAN_BEGIN, can't find: <PP> block."));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("In SCAN_END, can't find: </PP> block."));
    ifs.close();
}

TEST_F(GlobalFunctionTest, MapExist)
{
    std::map<int, double> SPIN = {{1, 2}, {3, 4}, {5, 6}};
    EXPECT_EQ(ModuleBase::GlobalFunc::MAP_EXIST(SPIN, 1), &SPIN[1]);
    EXPECT_EQ(ModuleBase::GlobalFunc::MAP_EXIST(SPIN, 3), &SPIN[3]);
    EXPECT_EQ(ModuleBase::GlobalFunc::MAP_EXIST(SPIN, 5), &SPIN[5]);
}

TEST_F(GlobalFunctionTest, ReadValue)
{
    ofs.open("tmp");
    ofs << "50 60 70" << std::endl;
    ofs.close();
    ifs.open("tmp");
    int nx, ny, nz;
    // source/module_cell/read_atoms.cpp line 153:154
    ifs >> nx >> ny;
    ModuleBase::GlobalFunc::READ_VALUE(ifs, nz);
    ifs.close();
    EXPECT_EQ(nz, 70);
}

TEST_F(GlobalFunctionTest, Dcopy)
{
    int size = 100;
    std::vector<std::complex<double>> aa(size, std::complex<double>(1.0, 2.0));
    std::vector<std::complex<double>> bb(size);
    std::vector<double> daa(size,1.1);
    std::vector<double> dbb(size);
    ModuleBase::GlobalFunc::DCOPY(aa, bb, size);
    ModuleBase::GlobalFunc::DCOPY(daa, dbb, size);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_COMPLEX_EQ(bb[i], aa[i]);
        EXPECT_DOUBLE_EQ(dbb[i], daa[i]);
    }
}

TEST_F(GlobalFunctionTest, VectorToPointer)
{
    int size = 100;
    std::vector<double> aa(size, 1.0);
    EXPECT_EQ(ModuleBase::GlobalFunc::VECTOR_TO_PTR(aa), aa.data());
    std::valarray<double> bb(1.0, size);
    EXPECT_EQ(ModuleBase::GlobalFunc::VECTOR_TO_PTR(bb), &bb[0]);
    const std::vector<double> cc(size, 1.0);
    EXPECT_EQ(ModuleBase::GlobalFunc::VECTOR_TO_PTR(cc), cc.data());
    const std::valarray<double> dd(1.0, size);
    EXPECT_EQ(ModuleBase::GlobalFunc::VECTOR_TO_PTR(dd), &dd[0]);
}

TEST_F(GlobalFunctionTest, COPYARRAY)
{
    long size = 100;
    std::complex<double>* aa = nullptr;
    std::complex<double>* bb = nullptr;
    aa = new std::complex<double>[size];
    bb = new std::complex<double>[size];
    std::complex<double> value{1.1, 2.2};
    std::fill(&aa[0], &aa[size], value);
    ModuleBase::GlobalFunc::COPYARRAY(aa,bb,size);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_COMPLEX_EQ(bb[i], value);
    }
    double* daa = nullptr;
    double* dbb = nullptr;
    daa = new double[size];
    dbb = new double[size];
    std::fill(&daa[0],&daa[size],3.3);
    ModuleBase::GlobalFunc::COPYARRAY(daa,dbb,size);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_DOUBLE_EQ(dbb[i], 3.3);
    }
    delete[] aa;
    delete[] bb;
    delete[] daa;
    delete[] dbb;
}
/*
TEST_F(GlobalFunctionTest, Note)
{
    ModuleBase::GlobalFunc::NOTE("Wrong Settings!");
    GlobalV::ofs_running.close();
    ifs.open("running.log");
    getline(ifs, output);
    getline(ifs, output);
    // output in runnint.log file
    EXPECT_THAT(output, testing::HasSubstr("Wrong Settings!"));
    ifs.close();
}
*/

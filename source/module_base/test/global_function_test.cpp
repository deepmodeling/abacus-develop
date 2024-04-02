#include "../global_function.h"
#include "../global_variable.h"
#include "../vector3.h"
#include "../blas_connector.h"
#include "../tool_quit.h"
#include <string>
#include <cstring>
#include <fstream>
#include <streambuf>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
 * - TO_STRING
 *   - change an arbitrary type of data to string
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
 * - IS_COLUMN_MAJOR_KS_SOLVER()
 *   - judge whether the KS_SOLVER is column major
 * - VECTOR_TO_PTR
 *   - get a copy of the ptr of a vector
 * - VECTOR_TO_PTR_v3double
 *    - get a copy of the ptr of a vector whose elements' type belongs to Vector3<double>
 * - MemAvailable
 *    - get the current memory valus
 * - TETS_LEVEL
 *    - set the test level
 * - BLOCK_HERE
 * 	  - add the block
 */

inline void EXPECT_COMPLEX_FLOAT_EQ(const std::complex<float>& a, const std::complex<float>& b)
{
    EXPECT_FLOAT_EQ(a.real(), b.real());
    EXPECT_FLOAT_EQ(a.imag(), b.imag());
}

inline void EXPECT_COMPLEX_DOUBLE_EQ(const std::complex<double>& a, const std::complex<double>& b)
{
    EXPECT_DOUBLE_EQ(a.real(), b.real());
    EXPECT_DOUBLE_EQ(a.imag(), b.imag());
}


template<typename T>
inline void CHECK_ZEROS(T &size)
{
    bool* pt_b = nullptr;
    int* pt_i = nullptr;
    float* pt_f = nullptr;
    double* pt_d = nullptr;
    std::complex<float>* pt_cf = nullptr;
    std::complex<double>* pt_cd = nullptr;
    ModuleBase::Vector3<double>* pt_v3 = nullptr;
    pt_b = new bool[size];
    pt_i = new int[size];
    pt_f = new float[size];
    pt_d = new double[size];
    pt_cf = new std::complex<float>[size];
    pt_cd = new std::complex<double>[size];
    pt_v3 = new ModuleBase::Vector3<double>[size];
    // long long size
    long long size_ll = 100;
    bool value_b = true;
    int value_i = 2;
    float value_f = 3.0;
    double value_d = 4.8;
    std::complex<float> value_cf{1.3, 2.2};
    std::complex<double> value_cd{1.1, 2.2};
    std::fill(&pt_b[0], &pt_b[size], value_b);
    std::fill(&pt_i[0], &pt_i[size], value_i);
    std::fill(&pt_f[0], &pt_f[size], value_f);
    std::fill(&pt_d[0], &pt_d[size], value_d);
    std::fill(&pt_cf[0], &pt_cf[size], value_cf);
    std::fill(&pt_cd[0], &pt_cd[size], value_cd);
    for (int i = 0; i < size; ++i)
    {
	    pt_v3[i].set(1.1,2.2,3.3);
    }
    ModuleBase::GlobalFunc::ZEROS(pt_b, size);
    ModuleBase::GlobalFunc::ZEROS(pt_i, size);
    ModuleBase::GlobalFunc::ZEROS(pt_f, size);
    ModuleBase::GlobalFunc::ZEROS(pt_d, size);
    ModuleBase::GlobalFunc::ZEROS(pt_cf, size);
    ModuleBase::GlobalFunc::ZEROS(pt_cd, size);
    ModuleBase::GlobalFunc::ZEROS(pt_v3, size);
    int zero_i = 0;
    float zero_f = 0.0;
    double zero_d = 0.0;
    std::complex<float> zero_cf{0.0, 0.0};
    std::complex<double> zero_cd{0.0, 0.0};
    for (int i = 0; i < size; ++i)
    {
        EXPECT_FALSE(pt_b[i]);
        EXPECT_EQ(pt_i[i],zero_i);
        EXPECT_FLOAT_EQ(pt_f[i],zero_f);
        EXPECT_DOUBLE_EQ(pt_d[i],zero_d);
        EXPECT_COMPLEX_FLOAT_EQ(pt_cf[i], zero_cf);
        EXPECT_COMPLEX_DOUBLE_EQ(pt_cd[i], zero_cd);
        EXPECT_DOUBLE_EQ(pt_v3[i].x,zero_d);
        EXPECT_DOUBLE_EQ(pt_v3[i].y,zero_d);
        EXPECT_DOUBLE_EQ(pt_v3[i].z,zero_d);
    }
    delete[] pt_b;
    delete[] pt_i;
    delete[] pt_f;
    delete[] pt_d;
    delete[] pt_cf;
    delete[] pt_cd;
    delete[] pt_v3;
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
    bool tmp_bool = true;
    int tmp_int = 1;
    long tmp_long = 2;
    unsigned long tmp_unsigned_long = 3;
    float tmp_float = 4.0;
    double tmp_double = 5.0;
    std::string tmp_string = "string";
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_bool", tmp_bool);
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_int", tmp_int);
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_long", tmp_long);
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_unsigned_long", tmp_unsigned_long);
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_float", tmp_float);
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_double", tmp_double);
    ModuleBase::GlobalFunc::OUT(ofs, "tmp_string", tmp_string);
    std::string para = "";
    int length = 0;
    for (int i=0;i<50;i++)
    {
    	para += "a";
	length = para.size()+1;
	if(length == 5){
		char tmp_char[5];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 6){
		char tmp_char[6];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 13){
		char tmp_char[13];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 15){
		char tmp_char[15];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 20){
		char tmp_char[20];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 22){
		char tmp_char[22];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 23){
		char tmp_char[23];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 25){
		char tmp_char[25];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 28){
		char tmp_char[28];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 29){
		char tmp_char[29];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 30){
		char tmp_char[30];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
	else if (length == 32){
		char tmp_char[32];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::OUT(ofs, "para", tmp_char);}
    }
    ofs.close();
    ifs.open("tmp");
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("tmp_bool = 1"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_int = 1"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_long = 2"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_unsigned_long = 3"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_float = 4"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_double = 5"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_string = string"));
    std::string tmp_para = "a";
    for (int i=0;i<50;i++)
    {
	tmp_para += "a";
	length = tmp_para.size()+1;
	if (length == 32) EXPECT_THAT(str, testing::HasSubstr(tmp_para));
    }
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
    EXPECT_THAT(output, testing::HasSubstr("grid = [ 100, 125, 375 ]"));
    getline(ifs, output);
    EXPECT_THAT(output, testing::HasSubstr("direct = [ 1.1, 2.2, 3.3 ]"));
    ifs.close();
}
// P for parameters
TEST_F(GlobalFunctionTest, OutP)
{
    ofs.open("tmp");
    bool tmp_bool = true;
    int tmp_int = 1;
    double tmp_double = 2.0;
    std::string tmp_string = "string";
    std::string tmp_bool_ex = "tmp_bool_ex";
    std::string tmp_int_ex = "tmp_int_ex";
    std::string tmp_double_ex = "tmp_double_ex";
    std::string tmp_string_ex = "tmp_string_ex";
    ofs << std::setiosflags(std::ios::left);
    ModuleBase::GlobalFunc::OUTP(ofs, "tmp_bool", tmp_bool, tmp_bool_ex);
    ModuleBase::GlobalFunc::OUTP(ofs, "tmp_int", tmp_int, tmp_int_ex);
    ModuleBase::GlobalFunc::OUTP(ofs, "tmp_double", tmp_double, tmp_double_ex);
    ModuleBase::GlobalFunc::OUTP(ofs, "tmp_string", tmp_string, tmp_string_ex);
    ofs.close();
    ifs.open("tmp");
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("tmp_bool                       1 #tmp_bool_ex"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_int                        1 #tmp_int_ex"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_double                     2 #tmp_double_ex"));
    EXPECT_THAT(str, testing::HasSubstr("tmp_string                     string #tmp_string_ex"));
    ifs.close();
}

TEST_F(GlobalFunctionTest, ToString)
{
    bool tmp_bool = true;
    int tmp_int = 1;
    long tmp_long = 2;
    unsigned long tmp_unsigned_long = 3;
    float tmp_float = 4.0;
    double tmp_double = 5.0;
    std::string tmp_string = "string";
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_bool),"1");
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_int),"1");
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_long),"2");
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_unsigned_long),"3");
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_float),"4");
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_double),"5");
    EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_string),"string");
    std::string para = "";
    int length = 0;
    for (int i=0;i<100;i++)
    {
    	para += "a";
	length = para.size()+1;
	if(length == 42){
		char tmp_char[42];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 47){
		char tmp_char[47];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 50){
		char tmp_char[50];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 52){
		char tmp_char[52];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 53){
		char tmp_char[53];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 63){
		char tmp_char[63];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 64){
		char tmp_char[64];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 74){
		char tmp_char[74];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 81){
		char tmp_char[81];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
	else if (length == 83){
		char tmp_char[83];
		strcpy(tmp_char,para.c_str());
    		EXPECT_EQ(ModuleBase::GlobalFunc::TO_STRING(tmp_char),para);}
    }
}

TEST_F(GlobalFunctionTest, MakeDir)
{
    GlobalV::MY_RANK = 0;
    ModuleBase::GlobalFunc::MAKE_DIR("scf");
    auto error1 = std::system("test -d ");
    EXPECT_EQ(error1, 0);
    auto error2 = std::system("rm -r scf ");
    EXPECT_EQ(error2, 0);
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
    bool tmp_b = false;
    int tmp_i = 1;
    float tmp_f = 2.0;
    double tmp_d = 3.0;
    std::string tmp_string = "string";
    ModuleBase::GlobalFunc::AUTO_SET("tmp_b", tmp_b);
    ModuleBase::GlobalFunc::AUTO_SET("tmp_i", tmp_i);
    ModuleBase::GlobalFunc::AUTO_SET("tmp_f", tmp_f);
    ModuleBase::GlobalFunc::AUTO_SET("tmp_d", tmp_d);
    ModuleBase::GlobalFunc::AUTO_SET("tmp_string", tmp_string);
    std::string para = "";
    int length = 0;
    for (int i=0;i<10;i++)
    {
    	para += "a";
	length = para.size()+1;
	if(length == 2){
		char tmp_char[2];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::AUTO_SET("tmp_char",tmp_char);}
	else if (length == 3){
		char tmp_char[3];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::AUTO_SET("tmp_char",tmp_char);}
	else if (length == 6){
		char tmp_char[6];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::AUTO_SET("tmp_char",tmp_char);}
	else if (length == 8){
		char tmp_char[8];
		strcpy(tmp_char,para.c_str());
    		ModuleBase::GlobalFunc::AUTO_SET("tmp_char",tmp_char);}
    }
    GlobalV::ofs_warning.close();
    ifs.open("warning.log");
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("AUTO_SET tmp_b to 0"));
    EXPECT_THAT(str, testing::HasSubstr("AUTO_SET tmp_i to 1"));
    EXPECT_THAT(str, testing::HasSubstr("AUTO_SET tmp_f to 2"));
    EXPECT_THAT(str, testing::HasSubstr("AUTO_SET tmp_d to 3"));
    EXPECT_THAT(str, testing::HasSubstr("AUTO_SET tmp_string to string"));
    EXPECT_THAT(str, testing::HasSubstr("AUTO_SET tmp_char to aaaaaaa"));
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

TEST_F(GlobalFunctionTest, Zeros)
{
    int size_i = 1000;
    CHECK_ZEROS(size_i);
    long size_l = 1000;
    CHECK_ZEROS(size_l);
    unsigned long size_ul = 1000;
    CHECK_ZEROS(size_ul);
    long long size_ll = 1000;
    CHECK_ZEROS(size_ll);
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
    ofs << "100" << std::endl;
    ofs << "3.0" << std::endl;
    ofs << "string" << std::endl;
    ofs.close();
    ifs.open("tmp");
    int tmp_int = 0;
    double tmp_double = 0.0;
    std::string tmp_string;
    // source/module_cell/read_atoms.cpp line 153:154
    ModuleBase::GlobalFunc::READ_VALUE(ifs, tmp_int);
    ModuleBase::GlobalFunc::READ_VALUE(ifs, tmp_double);
    ModuleBase::GlobalFunc::READ_VALUE(ifs, tmp_string);
    ifs.close();
    EXPECT_EQ(tmp_int, 100);
    EXPECT_DOUBLE_EQ(tmp_double, 3.0);
    EXPECT_EQ(tmp_string, "string");
}

TEST_F(GlobalFunctionTest, Dcopy)
{
    int size = 100;
    std::vector<std::complex<double>> aa(size, std::complex<double>(1.0, 2.0));
    std::vector<std::complex<double>> bb(size);
    std::vector<double> daa(size,1.1);
    std::vector<double> dbb(size);
    std::complex<double>* aalist = new std::complex<double>[size];
    std::complex<double>* bblist = new std::complex<double>[size];
    for (int i=0;i<size;i++)
    {
	    aalist[i] = std::complex<double>(1.0,2.0);
	    bblist[i] = std::complex<double>(0.0,0.0);
    }
    double* daalist = new double[size];
    double* dbblist = new double[size];
    for (int i=0;i<size;i++)
    {
	    daalist[i] = 3.0;
	    dbblist[i] = 0.0;
    }
    ModuleBase::GlobalFunc::DCOPY(aa, bb, size);
    ModuleBase::GlobalFunc::DCOPY(daa, dbb, size);
    ModuleBase::GlobalFunc::DCOPY(aalist, bblist, size);
    ModuleBase::GlobalFunc::DCOPY(daalist, dbblist, size);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_COMPLEX_DOUBLE_EQ(bb[i], aa[i]);
        EXPECT_COMPLEX_DOUBLE_EQ(bblist[i], aalist[i]);
        EXPECT_DOUBLE_EQ(dbb[i], daa[i]);
        EXPECT_DOUBLE_EQ(dbblist[i], daalist[i]);
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
        EXPECT_COMPLEX_DOUBLE_EQ(bb[i], value);
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

TEST_F(GlobalFunctionTest,IsColumnMajor)
{
	GlobalV::KS_SOLVER = "genelpa";
	EXPECT_TRUE(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER());
}

TEST_F(GlobalFunctionTest,Vector2Ptr)
{
    int size = 100;
    std::vector<std::complex<double>> aa(size, std::complex<double>(1.0, 2.0));
    std::complex<double>* ptr_d = nullptr;
    ptr_d=ModuleBase::GlobalFunc::VECTOR_TO_PTR(aa);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_COMPLEX_DOUBLE_EQ(ptr_d[i],std::complex<double>(1.0,2.0));
    }
}

TEST_F(GlobalFunctionTest,MemAvailable)
{
    for(int i=0;i<5;i++)
    {
        std::ifstream ifs("/proc/meminfo");
        while (ifs.good())
        {
            std::string label, size, kB;
            ifs >> label >> size >> kB;
            if (label == "MemAvailable:")
            {
                EXPECT_LE(std::stol(size)-1000,ModuleBase::GlobalFunc::MemAvailable());
                EXPECT_GE(std::stol(size)+1000,ModuleBase::GlobalFunc::MemAvailable());
            }
        }
    }
}

TEST_F(GlobalFunctionTest,TEST_LEVEL)
{
    std::string name;
    bool test_bool=false;
    name="none";
    ModuleBase::GlobalFunc::TEST_LEVEL(name,test_bool);
    EXPECT_EQ(GlobalV::test_wf,0);
    EXPECT_EQ(GlobalV::test_potential,0);
    EXPECT_EQ(GlobalV::test_charge,0);
    name="init_potential";
    ModuleBase::GlobalFunc::TEST_LEVEL(name,test_bool);
    EXPECT_EQ(GlobalV::test_wf,1);
    EXPECT_EQ(GlobalV::test_potential,1);
    EXPECT_EQ(GlobalV::test_charge,1);
    name="init_read";
    ModuleBase::GlobalFunc::TEST_LEVEL(name,test_bool);
    EXPECT_EQ(GlobalV::test_input,1);
    EXPECT_EQ(GlobalV::test_winput,1);
    EXPECT_EQ(GlobalV::test_kpoint,1);
    EXPECT_EQ(GlobalV::test_atom,1);
    EXPECT_EQ(GlobalV::test_unitcell,1);
#ifndef __EPM
        EXPECT_EQ(GlobalV::test_pseudo_cell,1);
#else
        EXPECT_EQ(test_epm_unitcell,1);
#endif
    name="pw_init";
    ModuleBase::GlobalFunc::TEST_LEVEL(name,test_bool);
    EXPECT_EQ(GlobalV::test_pw,1);
}

TEST_F(GlobalFunctionTest,BlockHere)
{
	std::string output2;
	std::string block_in="111";
	GlobalV::MY_RANK=1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::GlobalFunc::BLOCK_HERE(block_in), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr("\n********************************************"
		"\n Here is a Block, 1: go on 0: quit"
		"\n 111"
		"\n********************************************"));
}

TEST_F(GlobalFunctionTest,BlockHere2)
{
	std::string output2;
	std::string block_in="111";
	GlobalV::MY_RANK=0;
	std::string fake_input = "1";
	std::istringstream iss{fake_input};
	std::cin.rdbuf(iss.rdbuf());
	testing::internal::CaptureStdout();
//	EXPECT_EXIT(ModuleBase::GlobalFunc::BLOCK_HERE(block_in), ::testing::ExitedWithCode(0),"");
	ModuleBase::GlobalFunc::BLOCK_HERE(block_in);
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr("\n********************************************"
		"\n Here is a Block, 1: go on 0: quit"
		"\n 111"
		"\n********************************************"));
}

TEST_F(GlobalFunctionTest,BlockHere3)
{
	std::string output2;
	std::string block_in="111";
	GlobalV::MY_RANK=0;
	testing::internal::CaptureStdout();
	std::string fake_input = "0";
	std::istringstream iss{fake_input};
	std::cin.rdbuf(iss.rdbuf());
	EXPECT_EXIT(ModuleBase::GlobalFunc::BLOCK_HERE(block_in), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr("\n********************************************"
		"\n Here is a Block, 1: go on 0: quit"
		"\n 111"
		"\n********************************************"));
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

TEST_F(GlobalFunctionTest,Vector2Ptr_v3double)
{
    int size = 100;
    std::vector<ModuleBase::Vector3<double>> abcd(size, ModuleBase::Vector3<double>(1.1,2.2,3.3));
    ModuleBase::Vector3<double>* ptr_v3d = nullptr;
    ptr_v3d=ModuleBase::GlobalFunc::VECTOR_TO_PTR(abcd);
    for (int i = 0; i < size; ++i)
    {
        EXPECT_EQ(ptr_v3d[i],ModuleBase::Vector3<double>(1.1,2.2,3.3));
    }
}



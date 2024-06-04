#include "gtest/gtest.h"
#include "../write_wfc_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "../write_wfc_nao.h"
#include "module_base/scalapack_connector.h"

TEST(GenWfcLcaoFnameTest, OutType1GammaOnlyOutAppFlagTrue)
{
    int out_type = 1;
    bool gamma_only = true;
    bool out_app_flag = true;
    int ik = 0;
    int istep = 0;
    GlobalV::global_out_dir = "/path/to/global_out_dir/";

    std::string expected_output = "/path/to/global_out_dir/WFC_LCAO_GAMMA1.txt";
    std::string result = ModuleIO::wfc_lcao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);

    EXPECT_EQ(result, expected_output);
}

TEST(GenWfcLcaoFnameTest, OutType2GammaOnlyOutAppFlagFalse)
{
    int out_type = 2;
    bool gamma_only = true;
    bool out_app_flag = false;
    int ik = 1;
    int istep = 2;
    GlobalV::global_out_dir = "/path/to/global_out_dir/";

    std::string expected_output = "/path/to/global_out_dir/WFC_LCAO_GAMMA2_ION3.dat";
    std::string result = ModuleIO::wfc_lcao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);

    EXPECT_EQ(result, expected_output);
}

TEST(GenWfcLcaoFnameTest, OutTypeInvalid)
{
    int out_type = 3;
    bool gamma_only = false;
    bool out_app_flag = true;
    int ik = 2;
    int istep = 3;
    GlobalV::global_out_dir = "/path/to/global_out_dir/";

    std::string expected_output = "/path/to/global_out_dir/WFC_LCAO_K3.txt";
    // catch the screen output
    testing::internal::CaptureStdout();
    std::string result = ModuleIO::wfc_lcao_gen_fname(out_type, gamma_only, out_app_flag, ik, istep);
    std::string output = testing::internal::GetCapturedStdout();


    EXPECT_EQ(result, expected_output);
}


// mock write_wfc_nao as this function is implemented in another file
void ModuleIO::write_wfc_nao(const std::string &name, const double* ctot, const int nlocal, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary)
{
    int myrank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

    if (myrank == 0)
    {
        std::ofstream ofs(name, std::ios::binary);
        ofs.write((char*)ctot, nlocal * ekb.nc * sizeof(double));
        ofs.close();
    }
}

void ModuleIO::write_wfc_nao_complex(const std::string &name, const std::complex<double>* ctot, const int nlocal,const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary)
{
    int myrank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    if (myrank == 0)
    {
        std::ofstream ofs(name, std::ios::binary);
        ofs.write((char*)ctot, nlocal * ekb.nc * sizeof(std::complex<double>));
        ofs.close();
    }
}

template <typename T>
void read_bin(const std::string &name, std::vector<T> &data)
{
    std::ifstream ifs(name, std::ios::binary);
    ifs.seekg(0, std::ios::end);
    std::streamsize size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    size_t count = size / sizeof(T);
    data.resize(count);
    ifs.read(reinterpret_cast<char*>(data.data()), size);
    ifs.close();   
}

class WriteWfcLcaoTest : public testing::Test
{
protected:
    ModuleBase::matrix ekb;
    ModuleBase::matrix wg;
    std::vector<ModuleBase::Vector3<double>> kvec_c;
    std::vector<double> psi_init_double;
    std::vector<std::complex<double>> psi_init_complex;
    std::vector<double> psi_local_double;
    std::vector<std::complex<double>> psi_local_complex;
    Parallel_Orbitals pv;
    int nk = 2;
    int nbands = 3;
    int nbasis = 4;
    int nbands_local;
    int nbasis_local;
    int my_rank = 0;

    void SetUp() override
    {
        GlobalV::out_app_flag = 1;
        ekb.create(1, nbands); // in this test the value of ekb and wg is not important and not used.
        wg.create(1, nbands);
        psi_init_double.resize(nk * nbands * nbasis);
        psi_init_complex.resize(nk * nbands * nbasis);
        for (int i = 0; i < nk * nbands * nbasis; i++)
        {
            psi_init_double[i] = i * 0.1;
            psi_init_complex[i] = std::complex<double>(i * 0.2, i * 0.3);
        }
#ifdef __MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        pv.init(nbasis,nbands,1,MPI_COMM_WORLD);
        for (int i = 0; i < 9; i++)
        {
            pv.desc_wfc[i] = pv.desc[i];
        }

        // create a global PV, and distribute the psi_init to local
        Parallel_2D pv_glb;
        pv_glb.init(nbasis, nbands, std::max(nbasis,nbands), MPI_COMM_WORLD);
        psi_local_double.resize(nk * pv.get_row_size() * pv.get_col_size());
        psi_local_complex.resize(nk * pv.get_row_size() * pv.get_col_size());
        for(int ik=0;ik<nk;ik++)
        {
            Cpxgemr2d(nbasis, nbands, psi_init_double.data() + ik * nbands * nbasis, 1,1, pv_glb.desc, psi_local_double.data() + ik * pv.get_row_size() * pv.get_col_size(), 1,1, pv.desc, pv.blacs_ctxt);
            Cpxgemr2d(nbasis, nbands, psi_init_complex.data() + ik * nbands * nbasis, 1,1, pv_glb.desc, psi_local_complex.data() + ik * pv.get_row_size() * pv.get_col_size(), 1,1, pv.desc, pv.blacs_ctxt);
        }
        nbands_local = pv.get_col_size();
        nbasis_local = pv.get_row_size();
#else
        psi_local_double = psi_init_double;
        psi_local_complex = psi_init_complex;
        nbands_local = nbands;
        nbasis_local = nbasis;
#endif
    }
};

TEST_F(WriteWfcLcaoTest, WriteWfcLcao)
{
    // create a psi object
    psi::Psi<double> my_psi(psi_local_double.data(),nk,nbands_local,nbasis_local);
    GlobalV::global_out_dir = "./";
    ModuleIO::write_wfc_lcao(2,my_psi,ekb,wg,kvec_c,pv,-1) ;

    // check the output
    if (my_rank == 0)
    {
        for(int ik=0; ik < nk;ik++)
        {
            std::string fname = ModuleIO::wfc_lcao_gen_fname(2, true, GlobalV::out_app_flag, ik, -1);
            std::ifstream file1(fname);
            EXPECT_TRUE(file1.good());
            std::vector<double> data;
            read_bin(fname, data);

            EXPECT_EQ(data.size(), nbands*nbasis);

            for (int i = 0; i < nbands*nbasis; i++)
            {
                EXPECT_DOUBLE_EQ(data[i], psi_init_double[nbands*nbasis*ik + i]);
            }
            // remove the output files
            std::remove(fname.c_str());

        }
    }
}

TEST_F(WriteWfcLcaoTest, WriteWfcLcaoComplex)
{
    psi::Psi<std::complex<double>> my_psi(psi_local_complex.data(),nk,nbands_local,nbasis_local);
    GlobalV::global_out_dir = "./";
    ModuleIO::write_wfc_lcao(2,my_psi,ekb,wg,kvec_c,pv,-1) ;

    // check the output file
    if (my_rank == 0)
    {
        for(int ik=0; ik < nk;ik++)
        {
            std::string fname = ModuleIO::wfc_lcao_gen_fname(2, false, GlobalV::out_app_flag, ik, -1);
            std::ifstream file1(fname);
            EXPECT_TRUE(file1.good());
            std::vector<std::complex<double>> data;
            read_bin(fname, data);

            EXPECT_EQ(data.size(), nbands*nbasis);
            for (int i = 0; i < nbands*nbasis; i++)
            {
                EXPECT_DOUBLE_EQ(data[i].real(), psi_init_complex[nbands*nbasis*ik + i].real());
                EXPECT_DOUBLE_EQ(data[i].imag(), psi_init_complex[nbands*nbasis*ik + i].imag());
            }
            // remove the output files
            std::remove(fname.c_str());
        }
    }
}

#ifdef __MPI
#include "mpi.h"
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif


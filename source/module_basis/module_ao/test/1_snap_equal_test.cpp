#include<gtest/gtest.h>
#include"ORB_unittest.h"
#include "module_base/global_variable.h"

//Test whether the 2-center-int results
// and its derivative from two clases are equal. 
//	- ORB_gen_table::snap_psipsi(job=0) and Center2_Orb::Orb11::cal_overlap
// - ORB_gen_table::snap_psipsi(job=1) and Center2_Orb::Orb11::cal_grad_overlap
TEST_F(test_orb, equal_test)
{
	
	this->set_center2orbs();
	//equal test
	//orb
	double olm_0[1] = { 0 };
	double olm_1[3] = { 0,0,0 };
	//center2orb
    double clm_0 = 0;
	ModuleBase::Vector3<double> clm_1;

	//test parameters
	const double rmax = 5; //Ry
	srand((unsigned)time(NULL));
    ModuleBase::Vector3<double> R1(0, 0, 0);
	ModuleBase::Vector3<double> R2(randr(rmax), randr(rmax), randr(rmax));
	std::cout << "random R2=(" << R2.x << "," << R2.y << "," << R2.z << ")" << std::endl;
	ModuleBase::Vector3<double> dR = ModuleBase::Vector3<double>(0.001, 0.001, 0.001);
    //4. calculate overlap and grad_overlap by both methods
	int T1 = 0;
	
	for (int T2 = 0;T2 < ORB.get_ntype();++T2)
	{
		for (int L1 = 0;L1 < ORB.Phi[T1].getLmax();++L1)
		{
			for (int N1 = 0;N1 < ORB.Phi[T1].getNchi(L1);++N1)
			{
				for (int L2 = 0;L2 < ORB.Phi[T2].getLmax();++L2)
				{
					for (int N2 = 0;N2 < ORB.Phi[T2].getNchi(L2);++N2)
					{
						for (int m1 = 0;m1 < 2 * L1 + 1;++m1)
						{
							for (int m2 = 0;m2 < 2 * L2 + 1;++m2)
							{
								OGT.snap_psipsi(
									ORB, olm_0, 0, 'S',
									R1, T1, L1, m1, N1,
									R2, T2, L2, m2, N2
									);
								OGT.snap_psipsi(
									ORB, olm_1, 1, 'S',
									R1, T1, L1, m1, N1,
									R2, T2, L2, m2, N2
									);
								//std::cout << this->mock_center2_orb11[T1][T2][L1][N1][L2][N2]->cal_overlap(R1, R2, m1, m2);
								clm_0 =
									test_center2_orb11[T1][T2][L1][N1][L2][N2]->cal_overlap(R1, R2, m1, m2);
								clm_1 =
									test_center2_orb11[T1][T2][L1][N1][L2][N2]->cal_grad_overlap(R1, R2, m1, m2);
								EXPECT_NEAR(olm_0[0], clm_0, 1e-10);
								EXPECT_NEAR(olm_1[0], clm_1.x, 1e-10);
								EXPECT_NEAR(olm_1[1], clm_1.y, 1e-10);
								EXPECT_NEAR(olm_1[2], clm_1.z, 1e-10);
								ModuleBase::GlobalFunc::ZEROS(olm_1, 3);
							}
						}
					}

				}
			}
		}
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


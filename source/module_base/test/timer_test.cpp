#include "../timer.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <cstdio>
#include <chrono>
#include <thread>
#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of class timer
 ***********************************************/

/**
 * - Tested Functions:
 *   - Tick
 *     - tick 1st time, set start_flag to false
 *     - tick 2nd time, calculate time duration
 *   - Start
 *     - start total time calculation
 *   - PrintAll
 *     - print computational processes with time > 0.1 s 
 *   - Finish
 *     - finish total time calculation 
 *     - print computational processes with time > 0.1 s 
 *   - PrintUntilNow
 *     - stop total time calculation
 *     - print total time until now
 *     - then start total time calculation again
 */

class TimerTest : public testing::Test
{
protected:
	//
	// for capturing stdout
	std::string output;
	// for output in file
	std::ofstream ofs;
	std::ifstream ifs;
	int T_Elapse = 100; // microseconds = 0.1 milliseconds
	void TearDown()
	{
		remove("tmp");
	}
};


TEST_F(TimerTest, Tick)
{
	ModuleBase::timer::tick("wavefunc","evc");
	// after 1st call of tick, start_flag becomes false
	EXPECT_FALSE(ModuleBase::timer::timer_pool["wavefunc"]["evc"].start_flag);
	std::this_thread::sleep_for(std::chrono::microseconds(T_Elapse)); // 0.1 ms
	// then we can have time elapsed in cpu_second
	ModuleBase::timer::tick("wavefunc","evc");
	EXPECT_GT(ModuleBase::timer::timer_pool["wavefunc"]["evc"].cpu_second,0.0001);
}

TEST_F(TimerTest, Start)
{
	ModuleBase::timer::start();
	// start() called tick() once
	EXPECT_FALSE(ModuleBase::timer::timer_pool[""]["total"].start_flag);
}

TEST_F(TimerTest, PrintAll)
{
	ModuleBase::timer::tick("wavefunc","evc");
	std::this_thread::sleep_for(std::chrono::microseconds(T_Elapse)); // 0.1 ms
	ModuleBase::timer::tick("wavefunc","evc");
	// call print_all
	ofs.open("tmp");
	testing::internal::CaptureStdout();
	ModuleBase::timer::print_all(ofs);
	output = testing::internal::GetCapturedStdout();
	ofs.close();
	// checout output on screen
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	// check output in file
	ifs.open("tmp");
	getline(ifs,output);
	getline(ifs,output);
	getline(ifs,output);
	getline(ifs,output);
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	ifs.close();
}

TEST_F(TimerTest, PrintUntilNow)
{
	long double time = ModuleBase::timer::print_until_now();
	EXPECT_TRUE(time>0.0);
}

TEST_F(TimerTest, Finish)
{
	ModuleBase::timer::tick("wavefunc","evc");
	std::this_thread::sleep_for(std::chrono::microseconds(T_Elapse)); // 0.1 ms
	ModuleBase::timer::tick("wavefunc","evc");
	// call print_all
	ofs.open("tmp");
	testing::internal::CaptureStdout();
	ModuleBase::timer::finish(ofs);
	output = testing::internal::GetCapturedStdout();
	ofs.close();
	// checout output on screen
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	// check output in file
	ifs.open("tmp");
	getline(ifs,output);
	getline(ifs,output);
	getline(ifs,output);
	getline(ifs,output);
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	ifs.close();
}

// use __MPI to activate parallel environment
#ifdef __MPI
int main(int argc, char **argv)
{

	MPI_Init(&argc,&argv);

	testing::InitGoogleTest(&argc,argv);
	int result = RUN_ALL_TESTS();

	MPI_Finalize();

	return result;
}
#endif


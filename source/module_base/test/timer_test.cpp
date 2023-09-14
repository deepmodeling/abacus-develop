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

TEST_F(TimerTest, write_to_json)
{
	ModuleBase::timer::tick("wavefunc","evc");
	std::this_thread::sleep_for(std::chrono::microseconds(T_Elapse)); // 0.1 ms
	ModuleBase::timer::tick("wavefunc","evc");
	ModuleBase::timer::write_to_json("tmp.json");
	// check if tmp.json exists
	ifs.open("tmp.json");
	EXPECT_TRUE(ifs.good());

	// read all lines and remove all spaces and tabs and newlines
	std::string line;
	std::string tmp;
	std::string content;
	while(getline(ifs,line))
	{
		tmp = line;
		tmp.erase(std::remove(tmp.begin(),tmp.end(),' '),tmp.end());
		tmp.erase(std::remove(tmp.begin(),tmp.end(),'\t'),tmp.end());
		tmp.erase(std::remove(tmp.begin(),tmp.end(),'\n'),tmp.end());
		content += tmp;
	}

	// check if the content is correct
	// shuold be like this:
	// {"total": 0, "sub":[{"class_name": "wavefunc","sub":[{"name":"evc","cpu_second": 0.000318,"calls":2,"cpu_second_per_call":0.000159,"cpu_second_per_total": null}]}]}
	EXPECT_THAT(content,testing::HasSubstr("\"total\":"));
	EXPECT_THAT(content,testing::HasSubstr("\"sub\":[{\"class_name\":\"wavefunc\",\"sub\":[{\"name\":\"evc\",\"cpu_second\":"));
	EXPECT_THAT(content,testing::HasSubstr("\"calls\":2,\"cpu_second_per_call\":"));
	EXPECT_THAT(content,testing::HasSubstr("\"cpu_second_per_total\":"));
	EXPECT_THAT(content,testing::HasSubstr("}]}]}"));
	ifs.close();
	remove("tmp.json");
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
	std::cout << "Get captured stdout: \n" << std::endl;
	std::cout << output << std::endl;
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
	EXPECT_THAT(output,testing::HasSubstr("CLASS_NAME"));
	EXPECT_THAT(output,testing::HasSubstr("NAME"));
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("CALLS"));
	EXPECT_THAT(output,testing::HasSubstr("AVG(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("PER(%)"));
	// check output in file
	ifs.open("tmp");
	std::cout << "Capture contents line by line from output file: \n" << std::endl;
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
	getline(ifs,output);
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("CLASS_NAME"));
	EXPECT_THAT(output,testing::HasSubstr("NAME"));
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("CALLS"));
	EXPECT_THAT(output,testing::HasSubstr("AVG(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("PER(%)"));
	ifs.close();
	remove("time.json");
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
	std::cout << "Get captured stdout: \n" << std::endl;
	std::cout << output << std::endl;
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
	EXPECT_THAT(output,testing::HasSubstr("CLASS_NAME"));
	EXPECT_THAT(output,testing::HasSubstr("NAME"));
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("CALLS"));
	EXPECT_THAT(output,testing::HasSubstr("AVG(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("PER(%)"));
	// check output in file
	ifs.open("tmp");
	std::cout << "Capture contents line by line from output file: \n" << std::endl;
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
	getline(ifs,output);
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("CLASS_NAME"));
	EXPECT_THAT(output,testing::HasSubstr("NAME"));
	EXPECT_THAT(output,testing::HasSubstr("TIME(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("CALLS"));
	EXPECT_THAT(output,testing::HasSubstr("AVG(Sec)"));
	EXPECT_THAT(output,testing::HasSubstr("PER(%)"));
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


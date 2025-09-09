#include "timer.h"

#include <cmath>

#ifdef __MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

#include "chrono"
#include "source_base/formatter.h"

#if defined(__CUDA) && defined(__USE_NVTX)
#include <nvToolsExt.h>
#include "source_io/module_parameter/parameter.h"
#endif

namespace ModuleBase
{

//----------------------------------------------------------
// EXPLAIN :
//----------------------------------------------------------
bool timer::disabled = false;
size_t timer::n_now = 0;
std::map<std::string,std::map<std::string,timer::Timer_One>> timer::timer_pool;

void timer::finish(std::ofstream &ofs,const bool print_flag)
{
	timer::tick("","total");
	if(print_flag) {
		print_all( ofs );
}
}

//----------------------------------------------------------
//
//----------------------------------------------------------
void timer::start()
{
	// first init ,then we can use tick
	timer::tick("","total");
}

double timer::cpu_time()
{
//----------------------------------------------------------
// EXPLAIN : here static is important !!
// only first call can let t0 = 0,clock begin
// when enter this function second time , t0 > 0
//----------------------------------------------------------
	static auto t1 = std::chrono::system_clock::now();
	const auto t2 = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
	return double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
}

void timer::tick(const std::string &class_name,const std::string &name)
{
//----------------------------------------------------------
// EXPLAIN : if timer is disabled , return
//----------------------------------------------------------
	if (disabled) {
		return;
}

#ifdef _OPENMP
	if(!omp_get_thread_num())
#endif
	{
		Timer_One &timer_one = timer_pool[class_name][name];

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : cpu_time
//
// EXPLAIN :
// if start_flag == true,means a new clock counting begin,
// hence we record the start time of this clock counting.
// if start_flag == false, means it's the end of this counting,
// so we add the time during this two 'time point'  to the clock time storage.
//----------------------------------------------------------
		if(timer_one.start_flag)
		{
#ifdef __MPI
			int is_initialized = 0;
    		MPI_Initialized(&is_initialized);
			if(is_initialized)
			{
				timer_one.cpu_start = MPI_Wtime();
			}
#else
			timer_one.cpu_start = cpu_time();
#endif
			++timer_one.calls;
			timer_one.start_flag = false;
#if defined(__CUDA) && defined(__USE_NVTX)
            if (PARAM.inp.timer_enable_nvtx){
                std::string label = class_name + ":" + name;
                nvtxRangePushA(label.data());
            }
#endif
		}
		else
		{
#ifdef __MPI
			int is_initialized = 0;
    		MPI_Initialized(&is_initialized);
			if(is_initialized)
			{
				timer_one.cpu_second += MPI_Wtime() - timer_one.cpu_start;
			}
#else
			timer_one.cpu_second += (cpu_time() - timer_one.cpu_start);
#endif
			timer_one.start_flag = true;
#if defined(__CUDA) && defined(__USE_NVTX)
            if (PARAM.inp.timer_enable_nvtx){
                nvtxRangePop();
            }
#endif
		}
	} // end if(!omp_get_thread_num())
}

long double timer::print_until_now()
{
	// stop the clock
	timer::tick("","total");
	// start again
	timer::tick("","total");
	return timer_pool[""]["total"].cpu_second;
}

void timer::write_to_json(std::string file_name)
{
#ifdef __MPI
    // in some unit test, the mpi is not initialized, so we need to check it
	// if mpi is not initialized, we do not run this function
	int is_initialized = 0;
    MPI_Initialized(&is_initialized);
	if (!is_initialized) {
		return;
}
	int my_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if (my_rank != 0) {
		return;
}
#endif

    // check if a double is inf, if so, return "null", else return a string of the input double
	auto double_to_string = [](double d) -> std::string
	{
		if(std::isinf(d))
        {
			return "Infinity";
        }
		else
        {
			return FmtCore::format("%.15f", d);
        }
	};

	// The output json file format is like this:
	// {
	//     "total": 1.0,
	//     "sub": [
	//         {
	//             "class_name": "wavefunc",
	//             "sub": [
	//                 {
	//                     "name": "evc",
	//                     "cpu_second": 0.000318,
	//                     "calls": 2,
	//                     "cpu_second_per_call": 0.000159,
	//                     "cpu_second_per_total": 0.000318
	//                 }
	//             ]
	//         }
	//     ]
	// }

	std::ofstream ofs(file_name);
	std::string indent = "    ";
	int order_a = 0;
	ofs << "{\n";
	ofs << indent << "\"total\": " << timer_pool[""]["total"].cpu_second << ",\n";
	ofs << indent << "\"sub\": [\n";
	for(auto &timer_pool_A : timer_pool)
	{
		order_a ++;
		// if calss_name == "", it means total time, so we skip it
		if(timer_pool_A.first == "")
		{
			continue;
		}
		int order_b = 0;
		const std::string class_name = timer_pool_A.first;
		ofs << indent << indent << "{\n";
		ofs << indent << indent << indent << "\"class_name\": \"" << class_name << "\",\n";
		ofs << indent << indent << indent << "\"sub\": [\n";
		for(auto &timer_pool_B : timer_pool_A.second)
		{
			order_b ++;
			const std::string name = timer_pool_B.first;
			const Timer_One timer_one = timer_pool_B.second;
			ofs << indent << indent << indent << indent << "{\n";
			ofs << indent << indent << indent << indent << "\"name\": \"" << name << "\",\n";
			ofs << indent << indent << indent << indent << "\"cpu_second\": "
				<< std::setprecision(15) << timer_one.cpu_second << ",\n";
			ofs << indent << indent << indent << indent << "\"calls\": " << timer_one.calls << ",\n";
			ofs << indent << indent << indent << indent << "\"cpu_second_per_call\": "
				<< double_to_string(timer_one.cpu_second/timer_one.calls) << ",\n";
			ofs << indent << indent << indent << indent << "\"cpu_second_per_total\": "
				<< double_to_string(timer_one.cpu_second/timer_pool[""]["total"].cpu_second) << "\n";

			if (order_b == timer_pool_A.second.size())
			{
				ofs << indent << indent << indent << indent << "}\n";
			}
			else
			{
				ofs << indent << indent << indent << indent << "},\n";
			}
		}
		ofs << indent << indent << indent << "]\n";
		if (order_a == timer_pool.size())
		{
			ofs << indent << indent << "}\n";
		}
		else
		{
			ofs << indent << indent << "},\n";
		}
	}
	ofs << indent << "]\n";
	ofs << "}\n";
	ofs.close();
}

void timer::print_all(std::ofstream &ofs)
{
	constexpr double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6

	std::vector<std::pair<std::pair<std::string,std::string>,Timer_One>> timer_pool_order;
	for(auto &timer_pool_A : timer_pool)
	{
		const std::string class_name = timer_pool_A.first;
		for(auto &timer_pool_B : timer_pool_A.second)
		{
			const std::string name = timer_pool_B.first;
			const Timer_One timer_one = timer_pool_B.second;
			if(timer_pool_order.size() < timer_one.order+1)
			{
				timer_pool_order.resize(timer_one.order+1);
			}
			//timer_pool_order[timer_one.order] = {{class_name, name}, timer_one}; //qianrui change it to make it compatible with old compiler version
			timer_pool_order[timer_one.order] = std::pair<std::pair<std::string,std::string>, Timer_One> {
            std::pair<std::string,std::string >{class_name,name}, timer_one};
		}
	}
	std::vector<std::string> class_names;
	std::vector<std::string> names;
	std::vector<double> times;
	std::vector<int> calls;
	std::vector<double> avgs;
	std::vector<double> pers;
	for(auto &timer_pool_order_A : timer_pool_order)
	{
		const std::string &class_name = timer_pool_order_A.first.first;
		const std::string &name = timer_pool_order_A.first.second;
		const Timer_One &timer_one = timer_pool_order_A.second;

		if(timer_one.cpu_second < 0)
		{
			continue;
		}

		// only print out timers that are larger than 1%
		// mohan add 2025-03-09
		const double percentage_thr = 1.0;
        const double percentage = timer_one.cpu_second / timer_pool_order[0].second.cpu_second * 100;
		if(percentage<percentage_thr)
		{
			continue;
		}

		class_names.push_back(class_name);
		names.push_back(name);
		times.push_back(timer_one.cpu_second);
		calls.push_back(timer_one.calls);
		avgs.push_back(timer_one.cpu_second/timer_one.calls);


		// if the total time is too small, we do not calculate the percentage
		if (timer_pool_order[0].second.cpu_second < 1e-9)
		{
			pers.push_back(0);
		}
		else
		{
			pers.push_back(percentage);
		}
	}
	assert(class_names.size() == names.size());
	assert(class_names.size() == times.size());
	assert(class_names.size() == calls.size());
	assert(class_names.size() == avgs.size());
	assert(class_names.size() == pers.size());

	std::vector<std::string> titles = {"CLASS_NAME", "NAME", "TIME/s", "CALLS", "AVG/s", "PER/%"};
	std::vector<std::string> formats = {"%-10s", "%-10s", "%6.2f", "%8d", "%6.2f", "%6.2f"};
	FmtTable time_statistics(/*titles=*/titles,
							 /*nrows=*/pers.size(),
							 /*formats=*/formats,
							 /*indent=*/0,
							 /*align=*/{/*value*/FmtTable::Align::LEFT, /*title*/FmtTable::Align::CENTER});
	time_statistics << class_names << names << times << calls << avgs << pers;
	const std::string table = "\nTIME STATISTICS\n" + time_statistics.str();
	std::cout<<table<<std::endl;
	ofs<<table<<std::endl;
	write_to_json("time.json");

}
}

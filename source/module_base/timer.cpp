//==========================================================
// AUTHOR : fangwei , mohan
// DATE : 2008-11-06,
// UPDATE : Peize Lin at 2019-11-21
//==========================================================
#include "timer.h"
#include "chrono"
#include <vector>
#include <math.h>
#include "module_base/formatter.h"

#ifdef __MPI
#include "mpi.h"
#endif

#ifdef _OPENMP
#include <omp.h>
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
	if(print_flag)
		print_all( ofs );
}

//----------------------------------------------------------
//
//----------------------------------------------------------
void timer::start(void)
{
	// first init ,then we can use tick
	timer::tick("","total");
}

double timer::cpu_time(void)
{
//----------------------------------------------------------
// EXPLAIN : here static is important !!
// only first call can let t0 = 0,clock begin
// when enter this function second time , t0 > 0
//----------------------------------------------------------
	// static clock_t t0 = clock();
	// const clock_t t1 = clock() - t0;
	// return (t1<0) ? 0 : (double)t1/CLOCKS_PER_SEC;

	// static time_t t0 = time(NULL);
	// const time_t t1 = time(NULL);
	// double res = difftime(t1, t0);
	// return (res<0) ? 0 : res;
	static auto t1 = std::chrono::system_clock::now();
	const auto t2 = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
	return double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
		// mohan add, abandon the cross point time 2^32 ~ -2^32 .
}

void timer::tick(const std::string &class_name,const std::string &name)
{
//----------------------------------------------------------
// EXPLAIN : if timer is disabled , return
//----------------------------------------------------------
	if (disabled)
		return;

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
			int is_initialized;
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
		}
		else
		{
#ifdef __MPI
			int is_initialized;
    		MPI_Initialized(&is_initialized);
			if(is_initialized)
			{
				timer_one.cpu_second += MPI_Wtime() - timer_one.cpu_start;
			}
#else
			// if(class_name=="electrons"&&name=="c_bands")
			// {
			// 	cout<<"call times"<<timer_one.calls<<endl;
			// 	cout<<"electrons c_bands cost time:"<<endl;
			// 	cout<<cpu_time()<<"-"<<timer_one.cpu_start<<endl;
			// }

			timer_one.cpu_second += (cpu_time() - timer_one.cpu_start);
#endif
			timer_one.start_flag = true;
		}
	} // end if(!omp_get_thread_num())
}

long double timer::print_until_now(void)
{
	// stop the clock
	timer::tick("","total");
	// start again
	timer::tick("","total");
	return timer_pool[""]["total"].cpu_second;
}

void timer::write_to_json(std::string file_name)
{
	// check if a double is inf, if so, return "null", else return a string of the input double
#ifdef __MPI	
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if (my_rank != 0)
		return;
#endif
	auto double_to_string = [](double d) -> std::string
	{
		formatter::Fmt fmt(0, 15, ' ', false, false, false);
		if(isinf(d))
			return "Infinity";
		else
			return fmt.format(d);
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
			continue;
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
			ofs << indent << indent << indent << indent << "\"cpu_second\": " << std::setprecision(15) << timer_one.cpu_second << ",\n";
			ofs << indent << indent << indent << indent << "\"calls\": " << timer_one.calls << ",\n";
			ofs << indent << indent << indent << indent << "\"cpu_second_per_call\": " << double_to_string(timer_one.cpu_second/timer_one.calls) << ",\n";
			ofs << indent << indent << indent << indent << "\"cpu_second_per_total\": " << double_to_string(timer_one.cpu_second/timer_pool[""]["total"].cpu_second) << "\n";
			if (order_b == timer_pool_A.second.size())
				ofs << indent << indent << indent << indent << "}\n";
			else
				ofs << indent << indent << indent << indent << "},\n";
		}
		ofs << indent << indent << indent << "]\n";
		if (order_a == timer_pool.size())
			ofs << indent << indent << "}\n";
		else
			ofs << indent << indent << "},\n";
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
				timer_pool_order.resize(timer_one.order+1);
			//timer_pool_order[timer_one.order] = {{class_name, name}, timer_one}; //qianrui change it to make it compatible with old compiler version
			timer_pool_order[timer_one.order] = std::pair<std::pair<std::string,std::string>, Timer_One> {std::pair<std::string,std::string >{class_name,name}, timer_one};
		}
	}
	std::vector<std::string> class_names;
	std::vector<std::string> names;
	std::vector<double> times;
	std::vector<int> calls;
	std::vector<double> avgs;
	std::vector<double> pers;
	std::string table;
	context.set_context({"mid_title", "mid_title", "time", "int_w8", "time", "percentage"});
	for(auto &timer_pool_order_A : timer_pool_order)
	{
		const std::string &class_name = timer_pool_order_A.first.first;
		const std::string &name = timer_pool_order_A.first.second;
		const Timer_One &timer_one = timer_pool_order_A.second;

		if(timer_one.cpu_second < small)
			continue;
		class_names.push_back(class_name);
		names.push_back(name);
		times.push_back(timer_one.cpu_second);
		calls.push_back(timer_one.calls);
		avgs.push_back(timer_one.cpu_second/timer_one.calls);
		pers.push_back(timer_one.cpu_second / timer_pool_order[0].second.cpu_second * 100);
	}
	context.enable_title();
	context<<"CLASS_NAME"<<class_names<<"NAME"<<names<<"TIME(Sec)"<<times<<"CALLS"<<calls<<"AVG(Sec)"<<avgs<<"PER(%)"<<pers;
	context.center_title();
	context.set_overall_title("TIME STATISTICS");
	table = context.str();
	std::cout<<table<<std::endl;
	ofs<<table<<std::endl;
	write_to_json("time.json");
}
}

/*
void timer::print_all(std::ofstream &ofs)
{
//	std::cout<<"\n timer::print_all()"<<std::endl;
	const double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6

	std::cout << std::setprecision(2);

	// prepare
	bool *print_flag = new bool[n_clock];
	for(int i=0; i<n_clock; i++)
	{
		print_flag[i] = false;
	}

	int type = 1; // 2:calls 1:total_time
	bool non_reorder = 1;

	std::cout<<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << std::endl;
	ofs <<"\n\n\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << std::endl;
	ofs << std::setprecision(3);
	for (int i=0; i<n_clock; i++)
	{
		int k = 0;
		double tmp = -1.0;

		if(non_reorder)
		{
			k = i;
		}
		else
		{
			// search in all clocks
			for(int j=0; j<n_clock; j++)
			{
				if(print_flag[j])
				{
					continue;
				}
				if(type==1)
				{
					if(tmp < cpu_second[j])
					{
						k = j;
						tmp = cpu_second[j];
					}
				}
				else if(type==2)
				{
					if(tmp < calls[j])
					{
						k = j;
						tmp = calls[j];
					}
				}
			}
		}
		print_flag[k]=true;

		if ((cpu_second[k] >= 0 && cpu_second[k] < small) ||
		        (cpu_second[k] <= 0 && cpu_second[k] > -small))
		{
			continue;
		}

		if( level[k] > 'X' ) continue;


		const long double spend_time = cpu_second[k];
		const double average_spend_time = spend_time/calls[k];


		ofs  << " "
			 << std::setw(2) << level[k]
			 << std::setw(20) << class_name[k]
			 << std::setw(20) << name[k]
			 << std::setw(15) << spend_time
			 << std::setw(10) << calls[k]
			 << std::setw(10) << std::setprecision(2) << average_spend_time
			 << std::setw(10) << spend_time / cpu_second[0] * 100 << "%" << std::endl;


		std::cout << std::resetiosflags(ios::scientific);

		std::cout  << " "
		     << std::setw(2) << level[k]
			 << std::setw(20) << class_name[k]
			 << std::setw(20) << name[k]
			 << std::setw(15) << spend_time
			 << std::setw(10) << calls[k]
			 << std::setw(10) << std::setprecision(2) << average_spend_time
			 << std::setw(10) << spend_time / cpu_second[0] * 100 << "%" << std::endl;

	}
	std::cout<<" ----------------------------------------------------------------------------------------"<<std::endl;
	ofs <<" ----------------------------------------------------------------------------------------"<<std::endl;
	delete[] print_flag;
	return;
}
*/

#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace ModuleBase
{
/**
 * @brief Tracing computation time
 * @authors Fangwei, Mohan, Peize Lin
 *
 */
class timer
{
  public:
    struct Timer_One
    {
        double cpu_start;
        double cpu_second = 0.0;
        size_t calls = 0;
        size_t order = n_now++;
        bool start_flag = true;
    };

    static std::map<std::string, std::map<std::string, Timer_One>> timer_pool;

    /**
     * @brief Use twice at a time: the first time, set start_flag to false;
     * the second time, calculate the time duration
     *
     * @param class_name_in The class name for timing
     * @param name_in The compuational process for timing
     */
    static void tick(const std::string &class_name_in, const std::string &name_in);

    /**
     * @brief Start total time calculation
     *
     */
    static void start(void);

    /**
     * @brief Finish total time calculation and
     * print computational processes with duration > 0.1 s
     *
     * @param ofs The output file for print out timings
     * @param print_flag Print timings or not
     */
    static void finish(std::ofstream &ofs, const bool print_flag = 1);

    /**
     * @brief Enable time computation
     *
     */
    static void enable(void)
    {
        disabled = false;
    }

    /**
     * @brief Disable time computation
     *
     */
    static void disable(void)
    {
        disabled = true;
    }

    /**
     * @brief Write all computational processes to json file
     *
     * @param file_name The output file name
     */
    static void write_to_json(std::string file_name);

    /**
     * @brief Print all computational processes with during > 0.1 s
     *
     * @param ofs The output file for print out timings
     */
    static void print_all(std::ofstream &ofs);

    /**
     * @brief Stop total time calculation, print total time until now,
     * and then start total time calculation again
     *
     * @return long double
     */
    static long double print_until_now(void);

  private:
    /**
     * @brief Member variable: if disabled , timer can't work
     *
     */
    static bool disabled;

    /**
     * @brief Member variable: the index of clocks
     *
     */
    static size_t n_now;

    /**
     * @brief Member function: calculate time
     *
     * @return double
     */
    static double cpu_time(void);
};

} // namespace ModuleBase
#endif

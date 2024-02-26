#ifndef MEMORY_CUDA_H
#define MEMORY_CUDA_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <map>

namespace ModuleBase
{

/**
 * @brief Record memory consumption during computation.
 * @author Zhang Xiaoyang
 * @note 8 bit  = 1 Byte; 1024 Byte = 1 KB;
 * 1024 KB   = 1 MB; 1024 MB   = 1 GB
 *
 */
class Memory_CUDA
{
  public:
    Memory_CUDA();
    ~Memory_CUDA();

    /**
     * @brief Record memory consumed in gpu device during computation
     *
     * @param class_name The name of a class
     * @param name The name of a quantity
     * @param size The memory usage of this record
     * @param accumulate Useless, always set false
     * @return double
     */
    static double record(const std::string &class_name,
                         const std::string &name,
                         const long &size,
                         const bool accumulate = false);

    static double &get_total(void)
    {
        return total;
    }

    static void finish(std::ofstream &ofs);

    /**
     * @brief Print memory consumed in gpu (> 1 MB) in a file
     *
     * @param ofs The output file stream for print out memory records
     */
    static void print_all(std::ofstream &ofs);

    static void print(const int find_in);

  private:
    static double total;
    static double total_gpu;
    //static std::string *name;
    static std::map<int,double> name_mem_map;
    static std::map<int,double> name_mem_gpu_map;
    //static std::string *class_name;
    //static double *consume;
    static int n_memory;
    static int n_now;
    static int n_now_gpu;
    static bool init_flag;
    static bool gpu_init_flag;
};

} // namespace ModuleBase

#endif
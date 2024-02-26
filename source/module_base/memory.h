#ifndef MEMORY_H
#define MEMORY_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <map>
namespace ModuleBase
{

/**
 * @brief Record memory consumption during computation.
 * @author Mohan
 * @note 8 bit  = 1 Byte; 1024 Byte = 1 KB;
 * 1024 KB   = 1 MB; 1024 MB   = 1 GB
 *
 */
class Memory
{
  public:
    Memory();
    ~Memory();

    /**
     * @brief Record memory consumed during computation
     *
     * @param class_name The name of a class
     * @param name The name of a quantity
     * @param n The number of the quantity
     * @param type The type of data
     * @param accumulate Useless, always set false
     * @return double
     */
    static double record(const std::string &class_name,
                         const std::string &name,
                         const long &n,
                         const std::string &type,
                         const bool accumulate = false);

    /**
     * @brief Record memory consumed during computation
     *
     * @param name The name of a quantity
     * @param n The number of the quantity
     * @param accumulate Useless, always set false
     */
    static void record(
      const std::string &name_in,
      const size_t &n_in,
      const bool accumulate = false
    );

    /**
     * @brief Record memory consumed in gpu device during computation. This is transfered from original memory_cuda.h
     *
     * @param class_name The name of a class
     * @param name The name of a quantity
     * @param size The memory usage of this record
     * @param accumulate Useless, always set false
     */
    static double record_gpu(const std::string &class_name,
                         const std::string &name,
                         const long &size,
                         const bool accumulate = false);

    static double &get_total(void)
    {
        return total;
    }

    static void finish(std::ofstream &ofs);

    /**
     * @brief Print memory consumed (> 1 MB) in a file
     *
     * @param ofs The output file stream for print out memory records
     */
    static void print_all(std::ofstream &ofs);

    static void print(const int find_in);

    /**
     * @brief Calculate memory requirements for various
     * types of data
     *
     * @param n The number of a type of data
     * @param type The type of data
     * @return double
     */
    static double calculate_mem(const long &n, const std::string &type);

  private:
    static double total;
    static double total_gpu;
    //static std::string *name;
    static std::map<std::string,double> name_mem_map;
    static std::map<std::string,double> name_mem_gpu_map;
    static std::map<std::string,bool> name_print_flag_map;
    static std::map<std::string,bool> name_print_flag_gpu_map;
    //static std::string *class_name;
    //static double *consume;
    static bool init_flag;
    static bool init_flag_gpu;

    static int complex_matrix_memory; //(16 Byte)
    static int double_memory; //(8 Byte)
    static int int_memory; //(4 Byte)
    static int bool_memory;
    static int short_memory; //(2 Byte)
    static int float_memory; //(4 Byte)
};

} // namespace ModuleBase

#endif

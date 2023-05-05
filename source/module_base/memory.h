#ifndef MEMORY_H
#define MEMORY_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
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
     * @param name The name of a quantity
     * @param n The number of the quantity
     * @param accumulate Useless, always set false
     */
    static void record(
      const std::string &name_in,
      const size_t &n_in
    );

    /**
     * @brief Print memory consumed (> 1 MB) in a file
     *
     * @param ofs The output file stream for print out memory records
     */
    static void print_all(std::ofstream &ofs);

  private:
    static double total;
    static std::vector<std::string> name;
    static std::vector<double> consume;

    static void print(const int find_in);
};

} // namespace ModuleBase

#endif

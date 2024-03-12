//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
// LAST UPDATE : 2009-03-23 mohan add modify make_dir_out
//==========================================================
#ifndef GLOBAL_FILE_H
#define GLOBAL_FILE_H

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

//==========================================================
// namespace : Global_File_Operation
//==========================================================
namespace ModuleBase
{
namespace Global_File
{
	// called in input.cpp, after reading parameters.
	void make_dir_out(const std::string &suffix,
		const std::string &calculation,
        const bool &out_dir,
		const int rank,
        const bool &restart,
		const bool out_alllog = false); 
	
	void make_dir_atom(const std::string &label);
	void open_log ( std::ofstream &ofs, const std::string &fn, const std::string &calculation, const bool &restart);
	void close_log( std::ofstream &ofs, const std::string &fn);
	void close_all_log(const int rank, const bool out_alllog = false);

    /**
     * @brief delete tmperary files
     *
     */
    void delete_tmp_files();
}
}
#endif


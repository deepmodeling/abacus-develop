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
		const bool &out_wfc_dir,
		const int rank,
        const bool &restart,
		const bool out_alllog,
		const std::string &global_out_dir,
		const std::string &global_stru_dir,
		const std::string &global_matrix_dir,
		const std::string &global_wfc_dir,
		const std::string &global_mlkedf_descriptor_dir,
		const std::string &global_deepks_label_elec_dir,
		const std::string &log_file,
		const bool of_ml_gene_data,
		const bool deepks_out_freq_elec);

	void make_dir_atom(const std::string &label, const std::string &global_out_dir);
	void open_log ( std::ofstream &ofs, const std::string &fn, const std::string &calculation, const bool &restart, const std::string &global_out_dir);
	void close_log( std::ofstream &ofs, const std::string &fn);
    void close_all_log(const int rank, const bool out_alllog = false, const std::string& calculation = "md");
}
}
#endif


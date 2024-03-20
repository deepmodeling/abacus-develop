#ifndef RESTART_H
#define RESTART_H

#include <string>
#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#endif
class Restart
{
public:
	struct Info_Save
	{
		bool save_charge = false;
		bool save_H = false;
	};
	Info_Save info_save;
	
	struct Info_Load
	{
		bool load_charge = false;
		bool load_charge_finish = false;
		bool load_H = false;
		bool load_H_finish = false;
		bool restart_exx = false;
	};
	Info_Load info_load;
	
	std::string folder;
	
	void save_disk(const std::string mode, const int is, const int nrxx, double** rho) const;
	void load_disk(const std::string mode, const int is, const int nrxx, double** rho) const;
#ifdef __LCAO
    void save_disk(LCAO_Matrix &lm, const std::string mode, const int is, const int nrxx, double** rho) const;
    void load_disk(LCAO_Matrix &lm, const std::string mode, const int is, const int nrxx, double** rho) const;
#endif
private:
	void write_file1(const std::string &file_name, const void*const ptr, const size_t size) const;
	void read_file1(const std::string &file_name, void*const ptr, const size_t size) const;
	void write_file2(const std::string &file_name, const void*const ptr, const size_t size) const;
	void read_file2(const std::string &file_name, void*const ptr, const size_t size) const;
};

#endif
#ifndef RESTART_H
#define RESTART_H

#include <string>
#include "module_base/global_function.h"
#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#endif
class Restart
{
public:
	struct Info_Save
	{
		bool save_charge = false;
        bool save_H = false;    // save H means save Hexx now, will be changed in the future.
	};
	Info_Save info_save;
	
	struct Info_Load
	{
		bool load_charge = false;
		bool load_charge_finish = false;
		bool load_H = false;
		bool load_H_finish = false;
        bool restart_exx = false;   // to avoid the repeated load in MD/Relax
	};
	Info_Load info_load;
	
	std::string folder;

    template<typename T>
    bool save_disk(const std::string label, const int index, const int size, T* data, const bool error_quit = true) const
    {
        return write_file2(folder + label + "_" + ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK) + "_"
            + ModuleBase::GlobalFunc::TO_STRING(index),
            data,
            size * sizeof(T),
            error_quit);
    }
    template<typename T>
    bool load_disk(const std::string label, const int index, const int size, T* data, const bool error_quit = true) const
    {
        return read_file2(folder + label + "_" + ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK) + "_"
            + ModuleBase::GlobalFunc::TO_STRING(index),
            data,
            size * sizeof(T),
            error_quit);
    }
private:
	void write_file1(const std::string &file_name, const void*const ptr, const size_t size) const;
	void read_file1(const std::string &file_name, void*const ptr, const size_t size) const;
    bool write_file2(const std::string& file_name, const void* const ptr, const size_t size, const bool error_quit = true) const;
    bool read_file2(const std::string& file_name, void* const ptr, const size_t size, const bool error_quit = true) const;
};

#endif
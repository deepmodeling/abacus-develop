//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#ifndef READ_TXT_INIPUT_PROCESS_H
#define READ_TXT_INIPUT_PROCESS_H

#include "read_txt_input_list.h"

namespace Read_Txt_Input
{
	class Input_Process
	{
	public:
		Input_Process(Input_List &input_in) :input(input_in) {}
		void read_and_convert(const std::string &file_name);

	private:
		Input_List &input;

		void read(const std::string &file_name);
		void check_transform();
		void default_2();
		void out(const std::string &file_name) const;
		void bcast();
		void convert();

		void check_transform_global(std::map<std::string, Input_Item> &list);
		void default_2_global(std::map<std::string, Input_Item> &list);
		void convert_global(const std::map<std::string, Input_Item> &list);
	};
}

#endif
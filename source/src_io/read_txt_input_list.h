//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#ifndef READ_TXT_INIPUT_LIST_H
#define READ_TXT_INIPUT_LIST_H

#include "read_txt_input_item.h"

#include <vector>
#include <map>
#include <string>

namespace Read_Txt_Input
{
	class Input_List
	{
	public:
		void set_items();

	private:	
		void add_item(const Input_Item &input_item);
		//static void check_value_size(const Input_Item& input_item, const int size);
		//static void check_value_size(const Input_Item& input_item, const int size_low, const int size_up);

		std::map<std::string, Input_Item> list;
		std::vector<std::string> output_labels;

		void set_items_general();
		void set_items_pw();
		void set_items_spectrum();

		friend class Input_Process;
	};
}

#endif
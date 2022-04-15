//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#ifndef READ_TXT_INPUT_ITEM_H
#define READ_TXT_INPUT_ITEM_H

#include "read_txt_input_value.h"

#include <vector>
#include <map>
#include <string>
#include <functional>

#ifdef USE_CEREAL_SERIALIZATION
#include <cereal/access.hpp>
#endif

namespace Read_Txt_Input
{
	class Input_Item
	{
	public:
		// call these functions
		Input_Item(const std::string &label_in) :label(label_in) {}
		template<typename T>
		void default_1(const T &value);
		template <typename T_head, typename... T_tail>
		void default_1(const T_head &value_head, const T_tail... values_tail);

		// set these variables and functions
		std::string annotation;
		std::function<void(Input_Item&)> check_transform
			= [](Input_Item&self){};
		std::function<void(Input_Item&, const std::map<std::string, Input_Item>&)> default_2
			= [](Input_Item&self, const std::map<std::string, Input_Item>&list){};
		std::function<void(const Input_Item&)> convert
			= [](const Input_Item&item){};

		std::vector<Input_Value> values;

		void check_values_size(const int lower_limit, const int upper_limit)
		{ values_size_lower_limit=lower_limit;	values_size_upper_limit=upper_limit; }

	private:
		std::string label;
		std::vector<std::string> values_type;
		
		int values_size_read = -1;
		int values_size_lower_limit = 0;
		int values_size_upper_limit = 0;

		friend class Input_List;
		friend class Input_Process;

		template<typename T>
		void set_value(const T &value);

#ifdef USE_CEREAL_SERIALIZATION
	public:
		Input_Item() = default;	
	private:
		friend class cereal::access;
		template <class Archive> void serialize( Archive & ar );
#endif
	};	
}

#include "read_txt_input_item-template.h"

#endif
//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_item.h"

namespace Read_Txt_Input
{
	template<>
	void Input_Item::set_value(const bool &b)
	{
		this->values.push_back({});
		this->values.back().setb(b);
		this->values_type.push_back("b");
	}

	template<>
	void Input_Item::set_value(const int &i)
	{
		this->values.push_back({});
		this->values.back().seti(i);
		this->values_type.push_back("i");
	}

	template<>
	void Input_Item::set_value(const double &d)
	{
		this->values.push_back({});
		this->values.back().setd(d);
		this->values_type.push_back("d");
	}

	template<>
	void Input_Item::set_value(const std::string &s)
	{
		this->values.push_back({});
		this->values.back().sets(s);
		this->values_type.push_back("s");
	}

	template<>
	void Input_Item::set_value(const char*const &s)
	{
		this->values.push_back({});
		this->values.back().sets(std::string(s));
		this->values_type.push_back("s");
	}	
}
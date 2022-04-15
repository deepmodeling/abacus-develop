//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_value.h"

#include "src_io/read_txt_tools.h"
#include <stdexcept>

namespace Read_Txt_Input
{
	void Input_Value::setb(const bool &b_in)
	{
		this->s = std::to_string(b_in);
		this->d = static_cast<double>(b_in);
		this->i = static_cast<int>(b_in);
		this->b = b_in;
	}
	
	void Input_Value::seti(const int &i_in)
	{
		this->s = std::to_string(i_in);
		this->d = static_cast<double>(i_in);
		this->i = i_in;
		this->b = static_cast<bool>(i_in);
	}
	
	void Input_Value::setd(const double &d_in)
	{
		this->s = std::to_string(d_in);
		this->d = d_in;
		this->i = static_cast<int>(d_in);
		this->b = static_cast<bool>(d_in);
	}
	
	void Input_Value::sets(const std::string &s_in)
	{
		this->s = s_in;
		if(Read_Txt_Tools::in_set(s_in, Read_Txt_Tools::Preset::True))
		{
			this->b = true;
			this->i = 1;
			this->d = 1.0;
		}
		else if(Read_Txt_Tools::in_set(s_in, Read_Txt_Tools::Preset::False))
		{
			this->b = false;
			this->i = 0;
			this->d = 0.0;
		}
		else
		{
			try
			{
				this->d = std::stod(s_in);
				this->i = std::stoi(s_in);
				this->b = static_cast<bool>(std::stoi(s_in));
			}
			catch(const std::invalid_argument &e){}
		}
	}
}
//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#ifndef READ_TXT_INPUT_VALUE_H
#define READ_TXT_INPUT_VALUE_H

#include <string>

#ifdef USE_CEREAL_SERIALIZATION
#include <cereal/access.hpp>
#endif

namespace Read_Txt_Input
{
	class Input_Value
	{
	public:
		const bool & getb()const{ return b; }
		const int & geti()const{ return i; }
		const double & getd()const{ return d; }
		const std::string & gets()const{ return s; }

		void setb(const bool &b_in);
		void seti(const int &i_in);
		void setd(const double &d_in);
		void sets(const std::string &s_in);

	private:
		bool b = false;
		int i = 0;
		double d = 0.0;
		std::string s;

#ifdef USE_CEREAL_SERIALIZATION
	private:
		friend class cereal::access;
		template <class Archive> void serialize( Archive & ar );
#endif
	};
}

#endif
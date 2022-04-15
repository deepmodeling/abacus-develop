//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-08
//=======================

#include "read_txt_tools.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace Read_Txt_Tools
{
	// {
	//  	INPUT_PARAMETERS:       [1.0],
	//  	ecutwfc:        		[200, Ry]
	// }
	std::map<std::string, std::vector<std::string>> read_file_to_map(
		const std::string &file_name,
		const std::set<std::string> &comment_seps,
		const bool flag_check_duplicate)
	{
		std::string str;
		std::ifstream ifs(file_name);
		std::map<std::string, std::vector<std::string>> m;
		while(ifs.good())
		{
			std::getline(ifs,str);
			const std::vector<std::string> vec = Read_Txt_Tools::split_whitespace(Read_Txt_Tools::ignore_comment(str, comment_seps));
			if(vec.size()>0)
			{
				if(flag_check_duplicate)
					if(m.find(vec[0])!=m.end())
						throw std::invalid_argument(vec[0]+" is duplicated in INPUT");
				m[vec[0]] = std::vector<std::string>(vec.begin()+1, vec.end());
			}
		}
		return m;
	}



	// [
	//  	[INPUT_PARAMETERS, 1.0],
	//  	[ecutwfc, 200, Ry]
	// ]
	std::vector<std::vector<std::string>> read_file_to_vector(
		const std::string &file_name,
		const std::set<std::string> &comment_seps)
	{
		std::string str;
		std::ifstream ifs(file_name);
		std::vector<std::vector<std::string>> vec_all;
		while(ifs.good())
		{
			std::getline(ifs,str);
			const std::vector<std::string> vec = Read_Txt_Tools::split_whitespace(Read_Txt_Tools::ignore_comment(str, comment_seps));
			if(vec.size()>0)
				vec_all.push_back(vec);
		}
		return vec_all;
	}



	// [ecutwfc, 200, Ry]
	std::vector<std::string> split_whitespace(const std::string &str)
	{
		std::stringstream ss;
		ss<<str;
		std::vector<std::string> vec;
		std::string s;
		while(ss>>s)
			vec.push_back(s);
		return vec;
	}



	// str:
	//		ecutwfc 200 Ry # XXX
	// return:
	//		ecutwfc 200 Ry
	std::string ignore_comment(const std::string &str, const std::set<std::string> &comment_seps)
	{
		std::string str_new = str;
		for(const std::string &sep : comment_seps)
			str_new = str_new.substr(0, str_new.find(sep));
		return str_new;
	}



	// return_prefix + content_new = content_old
	std::vector<std::vector<std::string>> cut_paragraph(
		std::vector<std::vector<std::string>> &content,
		const std::set<std::string> &labels)
	{
		std::vector<std::vector<std::string>>::iterator ptr=content.begin()+1;
		for(; ptr<content.end(); ++ptr)
		{
			if(std::find(labels.begin(), labels.end(), (*ptr)[0]) != labels.end())
				break;
		}
		const std::vector<std::vector<std::string>> prefix = std::vector<std::vector<std::string>>(content.begin(), ptr);
		content = std::vector<std::vector<std::string>>(ptr, content.end());
		return prefix;
	}
}
//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_process.h"

#include "src_io/read_txt_tools.h"
#include "module_base/global_variable.h"

#ifdef USE_CEREAL_SERIALIZATION
#include "src_lcao/serialization_cereal.h"
#endif

namespace Read_Txt_Input
{
	void Input_Process::read_and_convert(const std::string &file_name)
	{
		if(GlobalV::MY_RANK==0)
		{
			this->read(file_name);
			this->check_transform();
			this->default_2();
			this->out(GlobalV::global_out_dir + file_name);
		}
	#ifdef __MPI
		this->bcast();
	#endif
		this->convert();
	}

	void Input_Process::read(const std::string &file_name)
	{
		const std::map<std::string, std::vector<std::string>> inputs_read
			= Read_Txt_Tools::read_file_to_map(file_name, {"#","\\"}, true);
		for(const auto & input_read : inputs_read)
		{
			const auto item_ptr = this->input.list.find(input_read.first);
			if(item_ptr==this->input.list.end())
				throw std::out_of_range("input_read.first");
			Read_Txt_Input::Input_Item &item = item_ptr->second;

			item.values_size_read = input_read.second.size();
			if(item.values.size()<input_read.second.size())
			{
				item.values.resize(input_read.second.size());
				item.values_type.resize(input_read.second.size(),"s");
			}

			for(size_t i=0; i<input_read.second.size(); ++i)
				item.values[i].sets(input_read.second[i]);
		}
	}

	void Input_Process::check_transform()
	{
		for(auto &tmp : this->input.list)
		{
			Read_Txt_Input::Input_Item &item = tmp.second;

			if( item.values_size_read>=0 &&
			   (item.values_size_read<item.values_size_lower_limit ||
			    item.values_size_read>item.values_size_upper_limit   ))
			{
				throw std::out_of_range(
					"size of INPUT "+item.label+" needs in ["+std::to_string(item.values_size_lower_limit)+","+std::to_string(item.values_size_upper_limit)+"]"
					+" but reads "+std::to_string(item.values_size_read));
			}

			item.check_transform(item);

			for(size_t i=0; i<item.values.size(); ++i)
			{
				if(item.values_type[i]=="b")
				{
					if(!Read_Txt_Tools::in_set( item.values[i].gets(), Read_Txt_Tools::Preset::Bool))
						throw std::invalid_argument("INPUT "+item.label+" must be bool");
				}
			}		
		}
		this->check_transform_global(this->input.list);
	}

	void Input_Process::default_2()
	{
		for(auto &item : this->input.list)
			item.second.default_2(item.second, this->input.list);
		this->default_2_global(this->input.list);
	}

	void Input_Process::out(const std::string &file_name) const
	{
		std::ofstream ofs(file_name);
		for(const std::string &label : this->input.output_labels)
		{
			const auto item_ptr = this->input.list.find(label);
			if(item_ptr==this->input.list.end())
			{
				ofs<<std::endl<<"# "<<label<<std::endl;
			}
			else
			{
				ofs<<label<<"\t";
				const Read_Txt_Input::Input_Item &item = item_ptr->second;
				for(size_t i=0; i<item.values.size(); ++i)
				{
					if(item.values_type[i]=="s")
						ofs<<item.values[i].gets()<<" ";
					else if(item.values_type[i]=="d")
						ofs<<std::to_string(item.values[i].getd())<<" ";
					else if(item.values_type[i]=="i")
						ofs<<std::to_string(item.values[i].geti())<<" ";
					else if(item.values_type[i]=="b")
					{
						if(item.values[i].getb())
							ofs<<"true"<<" ";
						else
							ofs<<"false"<<" ";
					}
					else
						throw std::invalid_argument("Input_Process::out() value_type["+std::to_string(i)+"]="+item.values_type[i]);
				}
				ofs<<"\t# "<<item.annotation<<std::endl;
			}
		}
	}

#ifdef __MPI
	void Input_Process::bcast()
	{
#ifdef USE_CEREAL_SERIALIZATION
		ModuleBase::bcast_data_cereal(this->input.list, MPI_COMM_WORLD, 0);
#else
#error Input_Process::bcast() needs cereal
#endif
	}
#endif

	void Input_Process::convert()
	{
		for(auto &item : this->input.list)
			item.second.convert(item.second);
		this->convert_global(this->input.list);
	}
}
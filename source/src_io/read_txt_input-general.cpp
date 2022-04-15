//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_list.h"

#include "src_io/read_txt_tools.h"
#include <stdexcept>

// for convert
#include "module_base/global_variable.h"
#include "src_pw/global.h"

namespace Read_Txt_Input
{
	void Input_List::set_items_general()
	{
		this->output_labels.push_back("Parameters (1.General)");

		{
			Input_Item item("nspin");
			item.default_1(1);
			item.annotation = "1: single spin; 2: up and down spin; 4: noncollinear spin";
			item.check_transform = [](Input_Item &self)
			{
				if(!Read_Txt_Tools::in_set(self.values[0].geti(), {1,2,4}))
					throw std::invalid_argument("nspin must be 1,2,4");
			};
			item.convert = [](const Input_Item &self)
			{
				GlobalV::NSPIN = self.values[0].geti();
			};
			this->add_item(item);
		}

		{
			Input_Item item("ntype");
			item.default_1(0);
			item.annotation = "atom species number";
			item.check_transform = [](Input_Item &self)
			{
				if(self.values[0].geti()<=0)
					throw std::invalid_argument("ntype must > 0");
			};
			item.convert = [](const Input_Item &self)
			{
				GlobalC::ucell.ntype = self.values[0].geti();
			};
			this->add_item(item);
		}

		{
			Input_Item item("symmetry");
			item.default_1(false);
			item.annotation = "turn symmetry on or off";
			item.convert = [](const Input_Item &self)
			{
				ModuleSymmetry::Symmetry::symm_flag = self.values[0].getb();
			};
			this->add_item(item);
		}
	}
}
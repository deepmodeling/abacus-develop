//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_list.h"

#include <limits>

namespace Read_Txt_Input
{
	void Input_List::set_items_spectrum()
	{
		this->output_labels.push_back("Parameters (15.Spectrum)");

		{
			Input_Item item("ocp_set");
			item.annotation = "set occupation number";
			item.check_values_size(0, std::numeric_limits<int>::max());
			item.convert = [](const Input_Item &self)
			{
				// ...
			};
			this->add_item(item);
		}
	}
}
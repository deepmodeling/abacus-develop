//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_list.h"

#include "module_base/constants.h"
#include <stdexcept>

// for convert
#include "module_base/global_variable.h"

namespace Read_Txt_Input
{
	void Input_List::set_items_pw()
	{
		this->output_labels.push_back("Parameters (2.Pw)");

		{
			Input_Item item("ecutwfc");
			item.default_1(100.0,"Ry");
			item.check_values_size(1,2);
			item.annotation = "energy cutoff for wave functions";
			item.check_transform = [](Input_Item &self)
			{
				if(self.values[0].getd()<=0)
					throw std::invalid_argument("ecutwfc must > 0");

				if(self.values[1].gets()=="Ry"){}
				else if(self.values[1].gets()=="eV")
					self.values[0].setd( self.values[0].getd() / ModuleBase::Ry_to_eV );
				else
					throw std::invalid_argument(self.values[1].gets());
				self.values[1].sets("Ry");
			};
			item.convert = [](const Input_Item &self)
			{
				// ?? GlobalC::pw.set()
			};
			this->add_item(item);
		}
	}
}
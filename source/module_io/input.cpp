#include "module_io/input.h"
#include "input_conv.h"
Input INPUT;

std::vector<int> Input::get_out_band_kb() const
{
    std::vector<int> out_band_kb;
    Input_Conv::parse_expression(bands_to_print_, out_band_kb);
    return out_band_kb;
}

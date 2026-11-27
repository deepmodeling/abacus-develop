//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-12-24
//==========================================================
#ifndef INPUT_CONVERT_H
#define INPUT_CONVERT_H

namespace Input_Conv
{

/**
 * @brief template bridge codes for converting string to other types
 *
 */
void tmp_convert();

/**
 * @brief Pass the data members from the INPUT instance(defined in
 * source_io/input.cpp) to GlobalV and GlobalC.
 */
void Convert();

} // namespace Input_Conv

#endif // Input_Convert

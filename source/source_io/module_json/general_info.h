#ifndef GENERAL_INFO_H
#define GENERAL_INFO_H

class Parameter;

/**
 * @brief Generate the general_info section of the JSON tree.
 */
namespace Json
{
#ifdef __JSON
void gen_general_info(const Parameter& param);
#endif
} // namespace Json

#endif

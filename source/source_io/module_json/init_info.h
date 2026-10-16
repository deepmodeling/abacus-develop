#ifndef INIT_INFO_H
#define INIT_INFO_H

class UnitCell;
struct Input_para;

/**
 * @brief Generate the init section of the JSON tree.
 */
namespace Json
{
#ifdef __JSON

void gen_init(UnitCell* ucell, const Input_para& inp);
void add_nkstot(int nkstot);
void gen_stru(UnitCell* ucell, const Input_para& inp);

#endif // __JSON
} // namespace Json

#endif

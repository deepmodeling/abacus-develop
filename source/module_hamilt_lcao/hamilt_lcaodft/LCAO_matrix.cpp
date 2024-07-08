#include "LCAO_matrix.h"

#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

LCAO_Matrix::LCAO_Matrix() {}

LCAO_Matrix::~LCAO_Matrix() {}

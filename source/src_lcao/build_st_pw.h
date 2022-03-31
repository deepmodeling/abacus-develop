#ifndef BUILD_ST_PW_H
#define BUILD_ST_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_lcao/LCAO_matrix.h"

class Build_ST_pw
{
	public:
	Build_ST_pw(LCAO_Matrix *lm);
	~Build_ST_pw();

	void set_ST(const int &ik, const char& dtype);
	void set_local(const int &ik);

private:
    LCAO_Matrix* LM;
};

#endif

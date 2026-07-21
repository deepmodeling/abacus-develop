#ifndef PARAM_UPDATE_H
#define PARAM_UPDATE_H

#include "source_cell/cal_atoms_info.h"

namespace elecstate {

class ParamUpdater {
public:
    static void update_from_atoms_info(const AtomsInfoResult& atoms_info);
};

}

#endif
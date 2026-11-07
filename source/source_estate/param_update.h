#ifndef PARAM_UPDATE_H
#define PARAM_UPDATE_H

#include "source_cell/cal_atoms_info.h"

class Parameter;

namespace elecstate {

class ParamUpdater {
public:
    static void update_from_atoms_info(const AtomsInfoResult& atoms_info);
    static void update_electron_count(double nelect, Parameter& parameters);
};

}

#endif

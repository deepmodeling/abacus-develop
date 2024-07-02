#include "parameter.h"

Parameter PARAM;

const Input_para& Parameter::get() const { return this->input; }
const MD_para& Parameter::get_mdp() const { return this->input.mdp; }
const System_para& Parameter::globalV() const { return this->gv; }

void Parameter::set_rank_nproc(const int& myrank, const int& nproc) {
    gv.myrank = myrank;
    gv.nproc = nproc;
}

void Parameter::set_start_time(const std::time_t& start_time) {
    gv.start_time = start_time;
}
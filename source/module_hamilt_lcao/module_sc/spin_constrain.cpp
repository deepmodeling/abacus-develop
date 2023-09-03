#include "spin_constrain.h"

SpinConstrain& SpinConstrain::getInstance() {
    static SpinConstrain instance; // Guaranteed to be created and destroyed only once
    return instance;
}

// set itia
void SpinConstrain::set_itia(const std::map<int, int>& itia) {
    this->itia = itia;
}

// get itia
std::map<int, int>& SpinConstrain::get_itia() {
    return this->itia;
}
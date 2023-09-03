#include "spin_constrain.h"

SpinConstrain& SpinConstrain::getInstance() {
    static SpinConstrain instance; // Guaranteed to be created and destroyed only once
    return instance;
}

void SpinConstrain::set_nat(int nat) {
    this->nat = nat;
}

int SpinConstrain::get_nat() {
    return this->nat;
}
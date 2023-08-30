#include "formatter_physfmt.h"
#include <cstring>

formatter::PhysicalFmt::PhysicalFmt(std::string context, Fmt* p_formatter) {
    context_ = context;
    if (p_formatter != nullptr) {
        this->p_formatter_ = p_formatter;
        this->decorator_mode_ = true;
    }
    else {
        this->p_formatter_ = new Fmt();
    }
    this->adjust_formatter();
}

formatter::PhysicalFmt::~PhysicalFmt() {

    if (this->p_formatter_ != nullptr && !this->decorator_mode_) {
        delete this->p_formatter_;
    }
}

void formatter::PhysicalFmt::adjust_formatter(bool left) {

    auto context = this->context_.c_str();
    if (strcmp(context, "none") == 0) {
        return;
    }
    else if (
        (strcmp(context, "int_w2") == 0)
      ||(strcmp(context, "kmesh") == 0)
      ||(strcmp(context, "constraint") == 0)
    ) {
        this->p_formatter_->set_width(2); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "int_w4") == 0)
      ||(strcmp(context, "scf_step") == 0)
    ) {
        this->p_formatter_->set_width(4); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(left);
    }
    else if (
        (strcmp(context, "int_w8") == 0)
    ) {
        this->p_formatter_->set_width(8); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w6_f1") == 0)
      ||(strcmp(context, "time") == 0)
        ) {
        this->p_formatter_->set_width(6); this->p_formatter_->set_precision(1);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w6_f2") == 0)
      ||(strcmp(context, "mass") == 0)
    ) {
        this->p_formatter_->set_width(6); this->p_formatter_->set_precision(2);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w8_f4_scientific") == 0)
      ||(strcmp(context, "threshold") == 0)
    ) {
        this->p_formatter_->set_width(8); this->p_formatter_->set_precision(4);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(false);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w8_f4_error") == 0)
      ||(strcmp(context, "charge") == 0)
    ) {
        this->p_formatter_->set_width(6); this->p_formatter_->set_precision(4);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left); this->p_formatter_->set_error(true);
    }
    else if (
        (strcmp(context, "double_w10_f2") == 0)
    ) {
        this->p_formatter_->set_width(10); this->p_formatter_->set_precision(2);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w16_f10") == 0)
      ||(strcmp(context, "coordinate") == 0)
      ||(strcmp(context, "position") == 0)
      ||(strcmp(context, "displacement") == 0)
      ||(strcmp(context, "lattice_constant") == 0)
      ||(strcmp(context, "lattice_vector") == 0)
      ||(strcmp(context, "lattice") == 0)
      ||(strcmp(context, "velocity") == 0)
    ) {
        this->p_formatter_->set_width(16); this->p_formatter_->set_precision(10);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w20_f10") == 0)
      ||(strcmp(context, "force") == 0)
      ||(strcmp(context, "energy") == 0)
    ) {
        this->p_formatter_->set_width(20); this->p_formatter_->set_precision(10);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "str_w4") == 0)
      ||(strcmp(context, "numbered_item") == 0)
    ) {
        this->p_formatter_->set_width(4); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "str_w30") == 0)
      ||(strcmp(context, "long_title") == 0)
    ) {
        this->p_formatter_->set_width(30); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(left);
    }
    else {
    }
}

void formatter::PhysicalFmt::set_context(std::string context) {
    context_ = context;
    this->adjust_formatter();
}

void formatter::PhysicalFmt::set_p_formatter(Fmt* p_formatter) {
    this->p_formatter_ = p_formatter;
    this->decorator_mode_ = true;
}

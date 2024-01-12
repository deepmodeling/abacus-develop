#include "formatter_fmt.h"


formatter::Fmt::Fmt() {
    width_ = 4;
    precision_ = 2;
    fillChar_ = ' ';
    fixed_ = true;
    right_ = true;
    error_ = false;
}

formatter::Fmt::Fmt(int width, int precision, char fillChar, bool fixed, bool right, bool error) {
    width_ = width;
    precision_ = precision;
    fillChar_ = fillChar;
    fixed_ = fixed;
    right_ = right;
    error_ = error;
}

formatter::Fmt::~Fmt() {}



void formatter::Fmt::reset() {
    width_ = 4;
    precision_ = 2;
    fillChar_ = ' ';
    fixed_ = true;
    right_ = true;
    error_ = false;
}

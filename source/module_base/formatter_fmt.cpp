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

// it's not good practice to overload such a function in a relative basic class, which will spoil the whole idea of inheritance
// so here I comment them two out
/*
template <typename T>
std::string formatter::Fmt::format(const std::vector<T>& value) {
    std::stringstream ss;
    for (auto v : value) {
        ss << this->format(v);
    }
    return ss.str();
}

template <typename T>
std::string formatter::Fmt::format(const T* value, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << this->format(value[i]);
    }
    return ss.str();
}
*/
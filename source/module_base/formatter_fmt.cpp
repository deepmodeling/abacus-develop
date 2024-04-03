#include "formatter_fmt.h"
#include <regex>

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

void formatter::Fmt::to_format(const std::string& strfmt,
                               int& width,
                               int& precision,
                               bool& fixed,
                               bool& right)
{
    std::regex re("^(%)(-)?(\\d+\\.)?(\\d+)?([dsfe])$");
    std::smatch m;
    if (std::regex_match(strfmt, m, re)) 
    {
        if (!m[1].matched) return; // do nothing 
        right = m[2].matched ? false : true;
        width = ((m[5] == "s")||(m[5] == "d")) ? std::stoi(m[4]) : std::stoi(m[3]);
        precision = ((m[5] == "s")||(m[5] == "d")) ? 0 : std::stoi(m[4]);
        fixed = (m[5] == "e") ? false : true;
    }
    else return;
}

void formatter::Fmt::reset() {
    width_ = 4;
    precision_ = 2;
    fillChar_ = ' ';
    fixed_ = true;
    right_ = true;
    error_ = false;
}

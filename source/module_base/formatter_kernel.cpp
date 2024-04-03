#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <cassert>

void to_format(const std::string& strfmt,
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

int main()
{
    // test to_format
    std::string strfmt = "%20.10f";
    int width, precision;
    bool fixed, right;
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 10);
    assert(fixed == true);
    assert(right == true);
    strfmt = "%-20.10f";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 10);
    assert(fixed == true);
    assert(right == false);
    strfmt = "%20.10e";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 10);
    assert(fixed == false);
    assert(right == true);
    strfmt = "%-20.10e";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 10);
    assert(fixed == false);
    assert(right == false);
    strfmt = "%20d";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 0);
    assert(fixed == true);
    assert(right == true);
    strfmt = "%-20d";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 0);
    assert(fixed == true);
    assert(right == false);
    strfmt = "%20s";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 0);
    assert(fixed == true);
    assert(right == true);
    strfmt = "%-20s";
    to_format(strfmt, width, precision, fixed, right);
    assert(width == 20);
    assert(precision == 0);
    assert(fixed == true);
    assert(right == false);
}
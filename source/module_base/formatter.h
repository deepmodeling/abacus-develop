#ifndef FORMATTER_H
#define FORMATTER_H

#include "formatter_contextfmt.h"
/*
* # ABACUS formatter library
* providing from basic the single variable formatting to complex vectors, tables formatting output
* Author: Yike HUANG, Daye ZHENG
* Institution: AISI
* ps. Zuxin JIN also helps solve dependency problems
* ## Classes
* - Fmt: single variable formatting (0 dimension)
* - PhysicalFmt: bundle with Fmt object, work as a decorator, adjusting the Fmt object according to the physical context, energy, force, coordinate, etc.
* - Table: output data structure, can be used to establish a table with titles, frames, etc.
* - ContextFmt: context formatting, providing the most sophisticated formatting
* ## Usage
* There is a Chinese Feishu document available for this library, please refer to it for more details:
* https://ucoyxk075n.feishu.cn/docx/Yym9dnm3aoTMfHxin8rcX9Rvnmb?theme=LIGHT&contentTheme=DARK
* ## Structure
* - File structure:
* formatter.h: header file
* formatter.cpp: source file
* test/formatter_fmt_test.cpp: unit test for Fmt class
* test/formatter_physfmt_test.cpp: unit test for PhysicalFmt class
* test/formatter_table_test.cpp: unit test for Table class
* test/formatter_contextfmt_test.cpp: unit test for ContextFmt class
* - Class structure:
*               Fmt: single variable formatting (0 dimension)
*                |
*                |
*                |
*          PhysicalFmt: decorator of Fmt
*                |
*                |
*                |
*           ContextFmt <----(inheritance)---- Table
*/
namespace formatter
{
}
/// @brief global variable, declared in module_base/formatter.h and defined in module_base/formatter_contextfmt.cpp for convenience, the same as INPUT
extern formatter::ContextFmt context;
#endif

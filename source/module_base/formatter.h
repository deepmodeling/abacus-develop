/**
 * @file formatter.h
 * @author HUANG Yike (huangyk@aisi.ac.cn)
 * @brief a new formatter library for formatting data
 * @version 2.0
 * @date 2024-04-23
 */
#ifndef FORMATTER_H
#define FORMATTER_H

#include <string>
#include <iostream>
#include <type_traits>
#include <utility>
#include <cstdio>
#include <numeric>
#include "module_base/ndarray.h"
/**
 * @brief 
 * 
 * In C++20, the std::format library is introduced. However, it is not supported under restriction of ABACUS development that not later than C++11. Plus in ABACUS the formatting-output demands is not quite general but more specific, therefore, a simple alternative is proposed here.
 * To use:
 * 1. Use the static function format() to format data like `FmtCore::format("%d", 1);`
 * 2. Use the class FmtCore to format data like `FmtCore fmt("%d"); fmt.format(1);`.
 * The first way is more flexible while the second way is more efficient. The format string can be reset by reset() function. If empty, the format string is empty, otherwise it will be updated.
 */
class FmtCore
{
public:
    FmtCore(const std::string& fmt): fmt_(fmt) {};
    ~FmtCore() {};
    /**
     * @brief static function to format data
     * 
     * @tparam Ts datatype of the data
     * @param fmt format string
     * @param args data to format
     * @return std::string 
     */
    template<typename... Ts>
    static inline std::string format(const char* fmt, const Ts&... args)
    {
        const int size = snprintf(nullptr, 0, fmt, FmtCore::filter(args)...) + 1;
        std::string dst(size, ' ');
        const int size_filled = snprintf(&dst[0], size, fmt, FmtCore::filter(args)...);
        dst.resize(size_filled);
        return dst;
    }
    /**
     * @brief std::string overload of the varadic template function
     * 
     * @param fmt 
     * @param arg 
     * @return std::string 
     */
    template<typename... Ts>
    std::string format(const Ts&... args) { return FmtCore::format(fmt_.c_str(), args...); }
    /**
     * @brief reset the format string (std::string overloads)
     * 
     * @param fmt 
     */
    void reset(const std::string& fmt = "") { fmt_ = fmt; }
    /**
     * @brief get the format string
     *  
     * @return std::string 
     */
    const std::string& fmt() { return fmt_; }
    /**
     * Python-style string functions will be implemented here as toolbox
     */
    /**
     * @brief split a string with a delimiter, return uncollapse vector
     * @param in string to split
     * @param delim delimiter
     * @return std::vector<std::string> 
    */
    static std::vector<std::string> split(const std::string& in, const std::string& delim)
    {
        std::vector<std::string> dst;
        std::string::size_type beg = 0, end;
        while((end = in.find(delim, beg)) != std::string::npos)
        {
            dst.push_back(in.substr(beg, end - beg));
            beg = end + delim.size();
        }
        dst.push_back(in.substr(beg));
        return dst;
    }
    /**
     * @brief split a string with a delimiter, return only non-empty elements
     * 
     * @param in string to split
     * @return std::vector<std::string> 
     */
    static std::vector<std::string> split(const std::string& in)
    {
        std::vector<std::string> dst;
        std::string::size_type beg = 0, end = 0;
        while((beg = in.find_first_not_of(" ", end)) != std::string::npos)
        {
            end = in.find_first_of(" ", beg);
            dst.push_back(in.substr(beg, end - beg));
        }
        return dst;
    }
    static bool startswith(const std::string& in, const std::string& prefix)
    {
        return (in.size() >= prefix.size()) && (in.substr(0, prefix.size()) == prefix);
    }
    static bool endswith(const std::string& in, const std::string& suffix)
    {
        return (in.size() >= suffix.size()) && (in.substr(in.size() - suffix.size()) == suffix);
    }
    static std::string strip(const std::string& in, const std::string& chars = " ")
    {
        std::string::size_type beg = in.find_first_not_of(chars);
        return (beg == std::string::npos)? "": in.substr(beg, in.find_last_not_of(chars) - beg + 1);
    }
    static std::string center(const std::string& in, const size_t& width, const char& fillchar = ' ')
    {
        if (in.size() >= width) {
            return in;
        }
        const size_t nwhitespaces = width - in.size();
        const size_t nleft = nwhitespaces / 2;
        const size_t nright = nwhitespaces - nleft;
        return std::string(nleft, fillchar) + in + std::string(nright, fillchar);
    }
    static std::string replace(const std::string& in, const std::string& old, const std::string& new_)
    {
        std::string dst = in;
        std::string::size_type pos = dst.find(old);
        while(pos != std::string::npos)
        {
            dst.replace(pos, old.size(), new_);
            pos = dst.find(old, pos + new_.size());
        }
        return dst;
    }
    static std::string join(const std::string& delim, const std::vector<std::string>& src)
    {
        return (src.empty())? "": std::accumulate(src.begin() + 1, src.end(), src[0], 
            [&delim](const std::string& acc, const std::string& s) { return acc + delim + s; });
    }

    static std::string upper(const std::string& in)
    {
        std::string dst = in;
        std::transform(dst.begin(), dst.end(), dst.begin(), ::toupper);
        return dst;
    }

    static std::string lower(const std::string& in)
    {
        std::string dst = in;
        std::transform(dst.begin(), dst.end(), dst.begin(), ::tolower);
        return dst;
    }

private:
    std::string fmt_;
    template<typename T>
    static typename std::enable_if<std::is_same<T, std::string>::value, const char*>::type filter(const T& s) { return s.c_str(); }
    template<typename T>
    static typename std::enable_if<!std::is_same<T, std::string>::value, const T&>::type filter(const T& s) { return s; }
};

class FmtTable
{
public:
    enum class Align{LEFT, RIGHT, CENTER};
private:
    typedef FmtCore core;
    struct Alignments{
        Alignments(const Align& val = Align::RIGHT, const Align& title = Align::CENTER): val_(val), title_(title) {};
        Align val_, title_; // value and title alignments
    } aligns_;
    struct Frames{
        Frames(char up = '-', char mid = '-', char dw = '-', char l = ' ', char r = ' '): up_(up), mid_(mid), dw_(dw), l_(l), r_(r) {};
        char up_, mid_ , dw_, l_, r_; // up, middle, down, left, right frames. up: the frame above title, middle: the one between title and data
    } frames_;                        // down: the frame below data, left: the frame at left, right: the frame at right
    struct Delimiters{
        Delimiters(char h = '-', char v = ' '): h_(h), v_(v) {};
        char h_, v_; // horizontal and vertical delimiters
    } delimiters_;
public:
    /**
     * @brief Construct a new Fmt Table object
     * 
     * @param titles titles, its size should be the same as the number of columns
     * @param nrows number of rows
     * @param aligns Alignments instance, can be constructed with initializer_list<char> like {'r', 'c'}, for right and center alignment for values and titles
     * @param frames Frames instance, can be constructed with initializer_list<char> like {'-', '-', '-', ' ', ' '}, for up, middle, down, left and right frames
     * @param delimiters Delimiters instance, can be constructed with initializer_list<char> like {'-', ' '}, for horizontal and vertical delimiters
     */
    FmtTable(const std::vector<std::string>& titles, 
             const size_t& nrows, 
             const std::vector<std::string>& fmts,
             const Alignments& aligns = {},
             const Frames& frames = {},
             const Delimiters& delimiters = {}): titles_(titles), data_(nrows, titles.size()), fmts_(fmts), aligns_(aligns), frames_(frames), delimiters_(delimiters)
    { assert(titles.size() == fmts.size()); };
    ~FmtTable() {};
    /**
     * @brief import data from std::vector
     * 
     * @tparam T datatype of the data
     * @param src source data
     * @return FmtTable& itself
     */
    template<typename T>
    FmtTable& operator<<(const std::vector<T>& src)
    {
        // create a copy of source data, then format
        std::vector<std::string> data(src.size());
        for(size_t i = 0UL; i < src.size(); i++) { data[i] = core::format(fmts_[j_].c_str(), src[i]); }
        set_value(0, j_, 'v', data);
        j_ = (j_ + 1) % titles_.size();
        return *this;
    }
    /**
     * @brief Set the value object
     * 
     * @tparam T datatype of the data
     * @param i row index
     * @param j col index
     * @param value value to set
     */
    template<typename T>
    void set_value(const size_t& i, const size_t& j, const T& value) { data_(i, j) = core::format(fmts_[j].c_str(), value); }
    /**
     * @brief adjust the width of each column
     * 
     * @param col col to relax, organized as std::vector<std::string>
     * @param title title of the column
     * @param vlyot value layout, can be Align::LEFT, Align::RIGHT, Align::CENTER
     * @param tlyot title layout, can be Align::LEFT, Align::RIGHT, Align::CENTER
     * @return std::vector<std::string> newly relaxed column
     */
    std::vector<std::string> relax_col_width(const std::vector<std::string>& col,
                                             const std::string& title = "",
                                             const Align valign = Align::RIGHT,  // because enum type would be the smallest integral type, so it is safe to pass by value
                                             const Align talign = Align::CENTER)
    {
        size_t max_width = title.size();
        for(const std::string& s : col) { max_width = std::max(max_width, s.size()); }
        std::vector<std::string> new_col(col.size() + 1); // the first is column title
        for(size_t i = 0; i < col.size() + 1; i++)
        {
            new_col[i] = (i == 0)? FmtCore::strip(title): FmtCore::strip(col[i - 1]);
            const size_t nwhitespaces = max_width - new_col[i].size();
            switch((i == 0)? talign: valign)
            {
                case Align::RIGHT: new_col[i] = std::string(nwhitespaces, ' ') + new_col[i]; break;
                case Align::LEFT: new_col[i] += std::string(nwhitespaces, ' '); break;
                case Align::CENTER: new_col[i] = FmtCore::center(new_col[i], max_width); break;
            }
        }
        return new_col;
    }

    /**
     * @brief concatenate titles into a string
     * 
     * @param titles titles to concatenate
     * @return std::string 
     */
    std::string concat_title(const std::vector<std::string>& titles) const
    {
        std::string dst;
        // first sum width of all titles
        size_t width = std::accumulate(titles.begin(), titles.end(), 0, [](const size_t& acc, const std::string& s) { return acc + s.size(); });
        // add width of delimiters
        width += titles.size() - 1;
        // add width of left and right frames
        width += 2;
        dst += std::string(width, frames_.up_) + "\n" + std::string(1, frames_.l_);
        for(size_t i = 0; i < titles.size(); i++)
        {
            dst += titles[i];
            if (i != titles.size() - 1) {
                dst += delimiters_.v_;
            }
        }
        dst += std::string(1, frames_.r_) + "\n" + std::string(width, frames_.mid_) + "\n";
        return dst;
    }
    /**
     * @brief concatenate a row into a string
     * 
     * @param row row to concatenate
     * @param pos position, can be 't' for top, 'b' for bottom, 'n' for normal
     * @return std::string 
     */
    std::string concat_row(const std::vector<std::string>& row, const char& pos) const
    {
        std::string dst = "";
        // first sum width of all elements of the row
        size_t width = std::accumulate(row.begin(), row.end(), 0, [](const size_t& acc, const std::string& s) { return acc + s.size(); });
        // for the delimiters
        width += row.size() - 1;
        // for the left and right frame
        width += 2;
        if (pos == 't') {
            dst += std::string(width, frames_.up_) + "\n";
        }
        dst += std::string(1, frames_.l_);
        for(size_t i = 0; i < row.size(); i++)
        {
            dst += row[i];
            if (i != row.size() - 1) {
                dst += delimiters_.v_;
            }
        }
        dst += std::string(1, frames_.r_) + "\n";
        if (pos == 'b') {
            dst += std::string(width, frames_.dw_) + "\n";
        }
        return dst;
    }
    /**
     * @brief to get the table as a string
     * 
     * @return std::string 
     */
    std::string str()
    {
        std::string dst = "";
        const size_t nrows = data_.shape()[0];
        const size_t ncols = data_.shape()[1];
        // if not all titles are empty, then with_title boolean will be true
        bool with_title = false;
        for (auto& title: titles_) {
            if (!title.empty()) {
                with_title = true;
                break;
            }
        }
        // first to relax each column
        for(size_t j = 0UL; j < ncols; j++)
        {
            std::vector<std::string> col(nrows);
            for (size_t i = 0UL; i < nrows; i++) {
                col[i] = data_(i, j);
            }
            col = relax_col_width(col, titles_[j], aligns_.val_, aligns_.title_);
            titles_[j] = col[0UL];
            std::vector<std::string> col_new(col.begin() + 1, col.end());
            set_value(0UL, j, 'v', col_new);
        }
        // then print titles
        if (with_title) {
            dst += concat_title(titles_);
        }
        // then print contents
        for(size_t i = 0UL; i < nrows; i++)
        {
            std::vector<std::string> row(ncols);
            for (size_t j = 0; j < ncols; j++) {
                row[j] = data_(i, j);
            }
            dst += concat_row(row, ((i == 0UL)&&!with_title)? 't': (i == nrows - 1)? 'b': 'n');
        }
        return dst;
    }
    void str(const std::string& s) {};
    // reuse
    void iter_set(const size_t val) { j_ = val; }
private:
    /**
     * @brief Set the value object from std::vector
     * 
     * @tparam T datatype of the data
     * @param i row index
     * @param j column index
     * @param dir direction, if 'v' then vertical, if 'h' then horizontal
     * @param src source data
     */
    template<typename T>
    void set_value(const size_t& i, const size_t& j, const char& dir, const std::vector<T>& src)
    {
        if (dir == 'v') {
            for (size_t k = 0UL; k < src.size(); k++) {
                data_(i + k, j) = src[k];
            }
        } else if (dir == 'h') {
            for (size_t k = 0UL; k < src.size(); k++) {
                data_(j, i + k) = src[k];
            }
        }
    }
    // iterator support indices
    size_t j_ = 0;

    std::vector<std::string> titles_;
    NDArray<std::string> data_; // data
    std::vector<std::string> fmts_; // format strings for each column
};

#endif
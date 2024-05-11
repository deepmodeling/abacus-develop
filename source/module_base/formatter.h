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
#include "./ndarray.h"
#include <iostream>
#include <type_traits>
#include <utility>
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
        std::string dst(size, '\0');
        snprintf(&dst[0], size, fmt, FmtCore::filter(args)...);
        dst.pop_back();
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

private:
    std::string fmt_;
    template<typename T>
    static typename std::enable_if<std::is_same<T, std::string>::value, const char*>::type filter(const T& s) { return s.c_str(); }
    template<typename T>
    static typename std::enable_if<!std::is_same<T, std::string>::value, const T&>::type filter(const T& s) { return s; }
};

class FmtTable
{
private:
    typedef FmtCore core;
    struct Frames{
        Frames(char upfrm = '-', char mdfrm = '-', char dwfrm = '-', char lfrm = ' ', char rfrm = ' '): upfrm_(upfrm), mdfrm_(mdfrm), dwfrm_(dwfrm), lfrm_(lfrm), rfrm_(rfrm) {};
        char upfrm_, mdfrm_ , dwfrm_, lfrm_, rfrm_;
    } frames_; 
    struct Delimiters{
        Delimiters(char hdlmt = '-', char vdlmt = ' '): hdlmt_(hdlmt), vdlmt_(vdlmt) {};
        char hdlmt_, vdlmt_;
    } delimiters_;
    struct Alignments{
        Alignments(char valign = 'r', char talign = 'c'): valign_(valign), talign_(talign) {};
        char valign_, talign_;
    } aligns_;

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
             const Delimiters& delimiters = {}): titles_(titles), fmts_(fmts), data_(nrows, titles.size()), aligns_(aligns), frames_(frames), delimiters_(delimiters)
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
     * @param vlyot value layout, can be 'l', 'r', 'c' for left, right and center
     * @param tlyot title layout, can be 'l', 'r', 'c' for left, right and center
     * @return std::vector<std::string> newly relaxed column
     */
    std::vector<std::string> relax_col_width(const std::vector<std::string>& col,
                                             const std::string& title = "",
                                             const char& vlyot = 'r',
                                             const char& tlyot = 'c')
    {
        size_t max_width = 0;
        for(const std::string& s : col) { max_width = std::max(max_width, s.size()); }
        if(!title.empty()) { max_width = std::max(max_width, title.size()); }
        std::vector<std::string> new_col(col.size() + static_cast<int>(!title.empty()));
        if(!title.empty())
        {
            const size_t nwhitespaces = max_width - title.size();
            const char tl_ = (tlyot == 'u') ? 'c' : tlyot;
            std::string dst = title;
            if(tl_ == 'r') dst = std::string(nwhitespaces, ' ') + title;
            else if(tl_ == 'l') dst = title + std::string(nwhitespaces, ' ');
            else if(tl_ == 'c')
            {
                const size_t nleft = nwhitespaces / 2;
                const size_t nright = nwhitespaces - nleft;
                dst = std::string(nleft, ' ') + title + std::string(nright, ' ');
            }
            new_col[0] = dst;
        }
        for(size_t i = 0; i < col.size(); i++)
        {
            const std::string src = col[i];
            // delete all whitespaces at left and right of src
            size_t ltr = src.find_first_not_of(' ');
            size_t rtr = src.find_last_not_of(' ');
            if(ltr == std::string::npos) ltr = 0;
            if(rtr == std::string::npos) rtr = src.size() - 1;
            std::string dst = src.substr(ltr, rtr - ltr + 1);
            const size_t nwhitespaces = max_width - src.size();
            if(vlyot == 'r') dst = std::string(nwhitespaces, ' ') + src;
            else if(vlyot == 'l') dst = src + std::string(nwhitespaces, ' ');
            else if(vlyot == 'c')
            {
                const size_t nleft = nwhitespaces / 2;
                const size_t nright = nwhitespaces - nleft;
                dst = std::string(nleft, ' ') + src + std::string(nright, ' ');
            }
            new_col[i + static_cast<int>(!title.empty())] = dst;
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
        std::string dst = "";
        // first sum width of all titles
        size_t width = std::accumulate(titles.begin(), titles.end(), 0, [](const size_t& acc, const std::string& s) { return acc + s.size(); });
        // for the delimiters
        width += titles.size() - 1;
        // for the left and right frame
        width += 2;
        dst += std::string(width, frames_.upfrm_) + "\n" + std::string(1, frames_.lfrm_);
        for(size_t i = 0; i < titles.size(); i++)
        {
            dst += titles[i];
            if(i != titles.size() - 1) dst += delimiters_.vdlmt_;
        }
        dst += std::string(1, frames_.rfrm_) + "\n" + std::string(width, frames_.mdfrm_) + "\n";
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
        if(pos == 't') dst += std::string(width, frames_.upfrm_) + "\n";
        dst += std::string(1, frames_.lfrm_);
        for(size_t i = 0; i < row.size(); i++)
        {
            dst += row[i];
            if(i != row.size() - 1) dst += delimiters_.vdlmt_;
        }
        dst += std::string(1, frames_.rfrm_) + "\n";
        if(pos == 'b') dst += std::string(width, frames_.dwfrm_) + "\n";
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
        // first to relax each column
        for(size_t j = 0UL; j < ncols; j++)
        {
            std::vector<std::string> col(nrows);
            for(size_t i = 0UL; i < nrows; i++) col[i] = data_(i, j);
            col = relax_col_width(col, titles_[j], aligns_.valign_, aligns_.talign_);
            std::string title = (titles_[j].empty())? "": col[0UL];
            titles_[j] = title;
            std::vector<std::string> col_new(col.begin() + static_cast<int>(!title.empty()), col.end());
            set_value(0UL, j, 'v', col_new);
        }
        // if not all titles are empty, then with_title boolean will be true
        bool with_title = false;
        for(auto& title : titles_) if(!title.empty()) { with_title = true; break; }
        if(with_title) dst += concat_title(titles_);
        // then print contents
        for(size_t i = 0UL; i < nrows; i++)
        {
            std::vector<std::string> row(ncols);
            for(size_t j = 0; j < ncols; j++) row[j] = data_(i, j);
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
        if(dir == 'v') for(size_t k = 0UL; k < src.size(); k++) data_(i + k, j) = src[k];
        else if(dir == 'h') for(size_t k = 0UL; k < src.size(); k++) data_(j, i + k) = src[k];
    }
    // iterator support indices
    size_t j_ = 0;

    std::vector<std::string> fmts_; // format strings for each column
    std::vector<std::string> titles_;
    NDArray<std::string> data_; // data
};

#endif
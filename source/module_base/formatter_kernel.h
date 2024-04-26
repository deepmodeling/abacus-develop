/**
 * @file formatter_kernel.cpp
 * @author HUANG Yike (huangyk@aisi.ac.cn)
 * @brief a new formatter library for formatting data
 * @version 2.0
 * @date 2024-04-23
 */
#ifndef FORMATTER_KERNEL_H
#define FORMATTER_KERNEL_H

#include <string>
#include "./ndarray.h"
#include <iostream>

/**
 * @brief 
 * 
 * In C++20, the std::format library is introduced. However, it is not supported under restriction of ABACUS development that not later than C++11. Plus in ABACUS the formatting-output demands is not quite general but more specific, therefore, a simple alternative is proposed here.
 * To use:
 * 1. Use the static function format() to format data like `ABACUSFormatter::format("%d", 1);`
 * 2. Use the class ABACUSFormatter to format data like `ABACUSFormatter fmt("%d"); fmt.format(1);`.
 * The first way is more flexible while the second way is more efficient. The format string can be reset by reset() function. If empty, the format string is empty, otherwise it will be updated.
 */
class ABACUSFormatter
{
public:
    ABACUSFormatter(const std::string& fmt): fmt_(fmt) {};
    ~ABACUSFormatter() {};
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
        // for there will not explicit fault can be detected by compiler,
        // but there seems some incompatilibility between std::string and const char*
        size_t buf_size = snprintf(nullptr, 0, fmt, args...);
        char* buf = new char[buf_size + 1];
        snprintf(buf, buf_size + 1, fmt, args...);
        std::string str(buf);
        delete[] buf;
        return str;
    }
    /**
     * @brief std::string overload of the varadic template function
     * 
     * @param fmt 
     * @param arg 
     * @return std::string 
     */
    static inline std::string format(const char* fmt, const std::string& arg)
    {
        return format(fmt, arg.c_str());
    }
    template<typename... Ts>
    std::string format(const Ts&... args)
    {
        int buf_size = snprintf(nullptr, 0, fmt().c_str(), args...);
        char* buf = new char[buf_size + 1];
        snprintf(buf, buf_size + 1, fmt().c_str(), args...);
        std::string str(buf);
        delete[] buf;
        return str;
    }
    /**
     * @brief reset the format string (const char* overloads)
     * 
     * @param fmt 
     */
    void reset(const char* fmt) { fmt_ = fmt; }
    /**
     * @brief reset the format string (std::string overloads)
     * 
     * @param fmt 
     */
    void reset(const std::string& fmt) { fmt_ = fmt; }
    /**
     * @brief clear the format string
     * 
     */
    void reset() { fmt_.clear(); }
    /**
     * @brief get the format string
     * 
     * @return std::string 
     */
    std::string fmt() const { return fmt_; }

private:
    std::string fmt_;
};

/**
 * @brief ABACUSTable class, for formatting data into a table
 * Due to implementation reasons, just declared here as key in pass-key idiom implementation
 * of other two subsidiary classes.
 */
class ABACUSTable;

class ABACUSTableStyle
{
public:
    // pass-key idiom
    ABACUSTableStyle(const ABACUSTable& table) {};
    ~ABACUSTableStyle() {};
    /**
     * @brief set the style of the table
     * 
     * @param attrb attribute name
     * @param value attribute value
     */
    void set(const std::string& attrb, const char& value)
    {
        if(attrb == "upfrm") upfrm_ = value;
        else if(attrb == "mdfrm") mdfrm_ = value;
        else if(attrb == "dwfrm") dwfrm_ = value;
        else if(attrb == "lfrm") lfrm_ = value;
        else if(attrb == "rfrm") rfrm_ = value;
        else if(attrb == "hdlmt") hdlmt_ = value;
        else if(attrb == "vdlmt") vdlmt_ = value;
        else if(attrb == "valign") valign_ = value;
        else if(attrb == "talign") talign_ = value;
    }
    /**
     * @brief relax the column width to the maximal width of the elements in the column
     * 
     * @param col source data, std::vector<std::string>
     * @param title optional, title of the column
     * @param vlyot optional, value layout, 'l' for left, 'r' for right, 'c' for center
     * @param tlyot optional, title layout, 'l' for left, 'r' for right, 'c' for center
     * @return std::vector<std::string> 
     */
    static std::vector<std::string> relax_col_width(const std::vector<std::string>& col,
                                                    const std::string& title = "",
                                                    const char& vlyot = 'r',
                                                    const char& tlyot = 'c')
    {
        size_t max_width = 0;
        for(const std::string& s : col) max_width = std::max(max_width, s.size());
        if(!title.empty()) max_width = std::max(max_width, title.size());
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
            std::string dst = src;
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
    std::string title(const std::vector<std::string>& titles,
                      const bool& dlmt = false) const
    {
        std::string dst = "";
        // first sum width of all titles
        size_t width = std::accumulate(titles.begin(), titles.end(), 0, [](const size_t& acc, const std::string& s) { return acc + s.size(); });
        // for the delimiters
        width += titles.size() - 1;
        // for the left and right frame
        width += 2;
        dst += upfrm(width);
        dst += "\n";
        dst += lfrm();
        for(size_t i = 0; i < titles.size(); i++)
        {
            dst += titles[i];
            if(i != titles.size() - 1) dst += dlmt? vdlmt(): ' ';
        }
        dst += rfrm();
        dst += "\n";
        dst += mdfrm(width);
        dst += "\n";
        return dst;
    }
    
    /**
     * @brief insert frame and delimiters on given row
     * 
     * @param src source data
     * @param pos position of present row, 't' for top, 'b' for bottom, omit otherwise
     * @return std::string 
     */
    std::string row(const std::vector<std::string>& src,
                    const char& pos) const
    {
        std::string dst = "";
        // first sum width of all elements of the row
        size_t width = std::accumulate(src.begin(), src.end(), 0, [](const size_t& acc, const std::string& s) { return acc + s.size(); });
        // for the delimiters
        width += src.size() - 1;
        // for the left and right frame
        width += 2;
        if(pos == 't') dst += upfrm(width) + "\n";
        dst += lfrm();
        for(size_t i = 0; i < src.size(); i++)
        {
            dst += src[i];
            if(i != src.size() - 1) dst += vdlmt();
        }
        dst += rfrm();
        dst += "\n";
        if(pos == 'b') dst += dwfrm(width) + "\n";
        return dst;
    }
    /**
     * @brief insert frame and delimiters on given column
     * 
     * @param src source data
     * @param pos position of present column, 'l' for left, 'r' for right, omit otherwise
     * @return std::string 
     */
    std::string col(const std::vector<std::string>& src,
                    const std::string& title,
                    const char& pos = 'n') const
    {
        std::string dst = "";
        throw std::runtime_error("Not implemented yet.");
        return dst;
    }
    char valign() const { return valign_; } // value alignment, 'l' for left, 'r' for right, 'c' for center
    char talign() const { return talign_; } // title alignment, 'l' for left, 'r' for right, 'c' for center
    char upfrm() const { return upfrm_; } // upframe, the frame above title
    std::string upfrm(const size_t& n) const { return std::string(n, upfrm_); }
    char mdfrm() const { return mdfrm_; } // midframe, the frame between title and content
    std::string mdfrm(const size_t& n) const { return std::string(n, mdfrm_); }
    char dwfrm() const { return dwfrm_; } // downframe, the frame at the end of table
    std::string dwfrm(const size_t& n) const { return std::string(n, dwfrm_); }
    char lfrm() const { return lfrm_; } // leftframe, the frame at the left of the table
    char rfrm() const { return rfrm_; } // rightframe, the frame at the right of the table
    char hdlmt() const { return hdlmt_; } // horizontal delimiter, between lines of content
    std::string hdlmt(const size_t& n) const { return std::string(n, hdlmt_); }
    char vdlmt() const { return vdlmt_; } // vertical delimiter, between columns of content
private:
    ABACUSTableStyle() {};
    // alignment
    char valign_ = 'r';
    char talign_ = 'c';
    // table style control
    char upfrm_ = '-';
    char mdfrm_ = '-';
    char dwfrm_ = '-';
    char lfrm_ = ' ';
    char rfrm_ = ' ';
    char hdlmt_ = '-';
    char vdlmt_ = ' ';
};

/**
 * @brief ABACUSTableContainer class, for storing the data of the table
 * - Only basic read/write operations are allowed, ModuleBase::NDArray is used to store the data
 * - Conceptually RAII, should define table size at the beginning.
 * - Cannot be created without ABACUSTable as the "key" in pass-key idiom.
 */
class ABACUSTableContainer
{
public:
    // pass-key idiom constructor, RAII
    ABACUSTableContainer(const ABACUSTable& table, const size_t nrows, const size_t ncols): nrows_(nrows), ncols_(ncols), titles_(ncols), table_(nrows, ncols) {};
    ~ABACUSTableContainer() {};
    /**
     * @brief set one title of the table
     * 
     * @param title title contents
     * @param i optional, the index of the title, default is 0
     */
    void assign_title(const std::string& title, const size_t j = 0) { titles_[j] = title; }
    /**
     * @brief set multiple titles of the table
     * 
     * @param titles title contents
     * @param i starting index of the titles, default is 0
     */
    void assign_title(const std::vector<std::string>& titles, const size_t j = 0)
    {
        // starting from j, to overwrite the titles
        #ifdef __DEBUG
        assert(j + titles.size() <= ncols_);
        #endif
        std::transform(titles.begin(), titles.end(), titles_.begin() + j, [](const std::string& title) { return title; });
    }
    /**
     * @brief assign a value to the table
     * 
     * @param value the value
     * @param i row index
     * @param j column index
     */
    void assign_value(const std::string& value, const size_t i = 0, const size_t j = 0) { table_(i, j) = value; }
    /**
     * @brief assign multiple values to the table
     * 
     * @param values values
     * @param direction 'v' for vertical, 'h' for horizontal
     * @param i starting row index
     * @param j column index
     */
    void assign_value(const std::vector<std::string>& values, const char& direction, const size_t i = 0, const size_t j = 0)
    {
        // if direction is v, means vertical, will update values to one column
        // if direction is h, means horizontal, will update values to one row
        // (i, j) defines the starting position
        if(direction == 'v')
        {
            // should assert i + values.size() <= nrows_
            #ifdef __DEBUG
            assert(i + values.size() <= nrows_);
            #endif
            for(size_t k = 0; k < values.size(); k++) { table_(i + k, j) = values[k]; }
        }
        else if(direction == 'h')
        {
            // should assert j + values.size() <= ncols_
            #ifdef __DEBUG
            assert(j + values.size() <= ncols_);
            #endif
            for(size_t k = 0; k < values.size(); k++) { table_(i, j + k) = values[k]; }
        }
    }

    std::string title(const size_t i) const { return titles_[i]; }
    std::vector<std::string> titles() const { return titles_; }
    std::string value(const size_t i, const size_t j) const { return table_(i, j); }

    /**
     * @brief get copy of row indiced by i of the table
     * 
     * @param i row index
     * @return std::vector<std::string> 
     */
    std::vector<std::string> row(const size_t i) const
    {
        std::vector<std::string> row;
        for(size_t j = 0; j < ncols_; j++) { row.push_back(table_(i, j)); }
        return row;
    }
    /**
     * @brief get copy of column indiced by j of the table
     * 
     * @param j column index
     * @return std::vector<std::string> 
     */
    std::vector<std::string> col(const size_t j) const
    {
        std::vector<std::string> col;
        for(size_t i = 0; i < nrows_; i++) { col.push_back(table_(i, j)); }
        return col;
    }
    size_t nrows() const { return nrows_; }
    size_t ncols() const { return ncols_; }
private:
    /**
     * @brief Construct a new abacus table object
     * This constructor is set as private to prevent the user from creating an object without
     * specifying the number of rows and columns.
     */
    ABACUSTableContainer() {};

    // attributes
    size_t nrows_ = 0;
    size_t ncols_ = 0;

    // container
    std::vector<std::string> titles_;
    NDArray<std::string> table_;
};

class ABACUSTable
{
private:
    // name is too long, use alias
    typedef ABACUSFormatter f;
    typedef ABACUSTableStyle stylizer;
public:
    ABACUSTable(const size_t& nrows, const size_t& ncols): table_(*this, nrows, ncols), style_(*this), fmts_(ncols) {};
    ABACUSTable(): table_(*this, 0, 0), style_(*this)
    {
        #ifdef __DEBUG
        printf("Default un-parameterized constructor is called.\n");
        assert(false);
        #endif
    };
    ~ABACUSTable() {};
    // step 1: set format
    // use a varadic parameter list to dynamically set format of each column, change value of fmts_
    template<typename... Ts>
    void fix_fmt(const Ts&... fmts) { fmts_ = {fmts...}; }
    void fix_fmt(const std::vector<std::string>& fmts) { fmts_ = fmts; }
    // step 2: import data and titles
    template<typename T>
    ABACUSTable& operator<<(const std::vector<T>& data)
    {
        for(size_t i = 0; i < data.size(); i++) // i is forever the row index
        {
            std::string s = f::format(fmts_[jval_].c_str(), data[i]); // need to cast std::string to const char*
            table_.assign_value(s, i, jval_);
        }
        jval_ = (jval_ + 1) % table_.ncols();
        return *this;
    }
    ABACUSTable& operator<<(const std::string& title)
    {
        table_.assign_title(title, jtitle_);
        jtitle_ = (jtitle_ + 1) % table_.ncols();
        return *this;
    }
    // step 3: set table frame, style, delimiters, ...
    ABACUSTableStyle& style() { return style_; }
    // step 4: output the table
    // unify column width. First calculate the maximal length of elements in each column
    // then for titles and table, fill the elements with spaces to make them have the same width
    std::string str()
    {
        std::string dst = "";
        // first to relax each column
        for(size_t j = 0; j < table_.ncols(); j++)
        {
            std::vector<std::string> col = table_.col(j);
            col = stylizer::relax_col_width(col, table_.title(j), style_.valign(), style_.talign());
            std::string title = (table_.title(j).empty())? "": col[0UL];
            std::vector<std::string> col_new(col.begin() + static_cast<int>(!title.empty()), col.end());
            table_.assign_value(col_new, 'v', 0UL, j);
            table_.assign_title(title, j);
        }
        // if not all titles are empty, then with_title boolean will be true
        bool with_title = false;
        for(auto& title : table_.titles()) if(!title.empty()) { with_title = true; break; }
        if(with_title) dst += style_.title(table_.titles());
        // then print contents
        for(size_t i = 0; i < table_.nrows(); i++)
        {
            dst += style_.row(table_.row(i), ((i == 0)&&!with_title)? 't': (i == table_.nrows() - 1)? 'b': 'n');
        }
        return dst;
    }
    void str(const std::string& s) {};
    // reuse
    void iter_reset() { jtitle_ = 0; jval_ = 0; }
    void iter_set(const char& domain, const size_t val) { if(domain == 't') jtitle_ = val; else jval_ = val; }
private:
    // iterator support indices
    size_t jtitle_ = 0;
    size_t jval_ = 0;
    // members
    ABACUSTableContainer table_; // container
    ABACUSTableStyle style_; // table formatting
    std::vector<std::string> fmts_; // format strings for each column
};

#endif
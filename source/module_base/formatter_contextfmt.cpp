#include "formatter_contextfmt.h"
#include <iostream>
#include <cstring>

formatter::ContextFmt::ContextFmt() {
    this->p_phys_fmt_ = new formatter::PhysicalFmt("none", &(this->fmt_));
    // because the table is called inside ContextFmt, so the mode of the table is set to 'data'
    this->disable_title();
}

formatter::ContextFmt::~ContextFmt() {
    delete this->p_phys_fmt_;
}

void formatter::ContextFmt::set_context(std::string context) {
    int iterative = this->iterative_;
    if (strcmp(context.c_str(), this->context_.c_str()) == 0) {
        this->iterative_ = iterative;
    }
    this->reset();
    this->context_ = context;
    auto it = this->predefined_phys_fmt.find(context);
    if (it != this->predefined_phys_fmt.end()) {
        this->known_context_ = true;
        for (auto fmt : it->second.second) {
            this->phys_fmt_.push_back(fmt);
            this->ncol_++;
        }
        if (it->second.first == 1) {
            if (this->iterative_ >= 1) {
                // iterative_ >= 1 means it is iterative mode
                this->enable_title();
            }
            else {
                // else for known context, the title is forbidden to input by crisscrossly calling << operator, say context << title << data
                // instead, use set_titles(std::vector<std::string> titles)
                this->disable_title();
            }
        }
        else if (it->second.first == 0) {
            // there is a title defined in this context, so it is needed to distinguish the title and the data
            this->enable_title();
        }
        else if (it->second.first == -1) {
            // I have not implemented this yet
            this->only_title();
        }
    }
    else {
        // unknown context, undefined behavior
        // warning: you would be better to set the context you like in formatter_contextfmt.h predefined_phys_fmt before using it
        this->known_context_ = false;
    }
}

void formatter::ContextFmt::set_context(int ncol, std::string* phys_fmt, int* nrows) {
    // this must be an unknown context
    this->reset();
    this->context_ = "";
    this->ncol_ = ncol;
    for (int icol = 0; icol < this->ncol_; ++icol) {
        this->nrows_.push_back(nrows[icol]);
        this->phys_fmt_.push_back(phys_fmt[icol]);
    }
}

void formatter::ContextFmt::set_context(std::string context, int ncol, int* &nrows) {
    this->set_context(context);
    this->ncol_ = ncol;
    for (int i = 0; i < ncol; ++i) {
        this->nrows_.push_back(nrows[i]);
    }
}

void formatter::ContextFmt::set_context(std::vector<std::string> phys_fmt) {
    // this must be an unknown context
    this->reset();
    this->context_ = "";
    this->ncol_ = phys_fmt.size();
    for (int icol = 0; icol < this->ncol_; ++icol) {
        this->phys_fmt_.push_back(phys_fmt[icol]);
    }
}

template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <double>(double const&);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <int>(int const&);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <char>(char const&);

template <>
formatter::ContextFmt& formatter::ContextFmt::operator<<(const std::string& value) {

    if (this->title_switch_%2 == 0) {
        // it is a title
        this->cache_title_ = value;  
    }
    else {
        // it is a data
        if (this->known_context_) {
            // known context, will be formatted
            this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
        }
        else {
            // unknown context, will not be formatted
            this->p_phys_fmt_->set_context("none");
        }
        Table::add_col(this->cache_title_, (std::vector<std::string>){this->fmt_.format(value)});
        this->icol_++;
        this->fmt_.reset();
    }
    if (Table::get_mode() == 1) {
        this->title_switch_ += 2;
    }
    else {
        this->title_switch_++;
    }
    return *this;
}

formatter::ContextFmt& formatter::ContextFmt::operator<< (char const* value)
{
    std::string value_(value);
    return *this << value_;
}

template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <double>(std::vector<double> const& value);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <std::string>(std::vector<std::string> const& value);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <int>(std::vector<int> const& value);

template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <double>(double*& value);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <int>(int*& value);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <std::string>(std::string*& value);
template<> formatter::ContextFmt& formatter::ContextFmt::operator<< <char>(char*& value);

void formatter::ContextFmt::reset() {
    Table::reset();
    //this->context_ = "";
    this->phys_fmt_.clear();
    this->default_phys_fmt_ = "energy";
    this->cache_title_ = "";
    this->with_title_ = true;

    this->ncol_ = 0;
    this->icol_ = 0;
    this->nrows_.clear();
    this->iterative_ = 0;

    this->fmt_.reset();
    this->title_switch_ = 0;

    this->known_context_ = false;

    this->disable_title(); // do what constructor does
}

void formatter::ContextFmt::context_refresh() {
    this->reset();
    this->set_context(this->context_);
}

std::string formatter::ContextFmt::str(bool next_line) {
    if (this->iterative_ == 1)
    {
        Table::set_mode(0);
        Table::disable_down_frame();
        this->iterative_ += 1;
    }
    else if (this->iterative_ > 1)
    {
        Table::set_mode(1);
        this->iterative_ += 1;
    }
    std::string str = this->print_table();
    if (!next_line) {
        str.pop_back();
    }
    this->context_refresh();
    return str;
}

void formatter::ContextFmt::print_status() const {
    std::cout << "-----------------status-----------------" << std::endl;
    std::cout << "context: " << this->context_ << std::endl;
    std::cout << "known_context: " << this->known_context_ << std::endl;
    std::cout << "ncol: " << this->ncol_ << std::endl;
    std::cout << "nrows: ";
    for (auto n : this->nrows_) {
        std::cout << n << " ";
    }
    std::cout << std::endl;
    std::cout << "phys_fmt: ";
    for (auto fmt : this->phys_fmt_) {
        std::cout << fmt << " ";
    }
    std::cout << std::endl;
    std::cout << "cache_title: " << this->cache_title_ << std::endl;
    std::cout << "icol: " << this->icol_ << std::endl;
    std::cout << "title_switch: " << this->title_switch_ << std::endl;
    std::cout << "iterative: " << this->iterative_ << std::endl;
    std::cout << "fmt: " << std::endl;
    std::cout << "width: " << this->fmt_.get_width() << std::endl;
    std::cout << "precision: " << this->fmt_.get_precision() << std::endl;
    std::cout << "fillChar: " << this->fmt_.get_fillChar() << std::endl;
    std::cout << "fixed: " << this->fmt_.get_fixed() << std::endl;
    std::cout << "right: " << this->fmt_.get_right() << std::endl;
    std::cout << "error: " << this->fmt_.get_error() << std::endl;
    std::cout << "p_phys_fmt: " << this->p_phys_fmt_ << std::endl;
    std::cout << "mode: " << Table::get_mode() << std::endl;
    std::cout << "col_delimiter: " << Table::get_col_delimiter() << std::endl;
    std::cout << "frame_switches: ";
    for (auto s : Table::get_frame_switches()) {
        std::cout << s << " ";
    }
    std::cout << std::endl;
    std::cout << "frame_delimiters: ";
    for (auto d : Table::get_frame_delimiters()) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
    std::cout << "flexible_width: " << Table::get_flexible_width() << std::endl;
    std::cout << "frame_mid_switch: " << Table::get_frame_mid_switch() << std::endl;
    std::cout << "frame_mid_delimiter: " << Table::get_frame_mid_delimiter() << std::endl;
    std::cout << "col_widths: ";
    for (auto w : Table::get_col_widths()) {
        std::cout << w << " ";
    }
    std::cout << std::endl;
    std::cout << "col_max_width: " << Table::get_col_max_width() << std::endl;
    std::cout << "total_width: " << Table::get_total_width() << std::endl;
    std::cout << "----------------------------------------" << std::endl;
}

formatter::ContextFmt context;
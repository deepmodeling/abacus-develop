#include "formatter_table.h"
#include <iomanip>
#include <sstream>

int formatter::Table::add_col(std::string new_title, std::vector<std::string> new_col) {
    this->ncol_++;
    this->titles_.push_back(new_title);
    this->data_.push_back(new_col);
    return this->ncol_-1;
}

void formatter::Table::adjust_col_width() {

    this->col_widths_.clear();
    int max_width = 0;
    for (int icol = 0; icol < this->ncol_; ++icol) {
        int s1 = this->titles_[icol].size();
        int s2 = 0;
        for (auto row : this->data_[icol]) {
            if (row.size() > s2) {
                s2 = row.size();
            }
        }
        if (s1 > s2) {
            if (flexible_width_) {
                this->col_widths_.push_back(s1);
            }
            else {
                if (s1 > max_width) {
                    max_width = s1;
                }
            }
        }
        else {
            if (flexible_width_) {
                this->col_widths_.push_back(s2);
            }
            else {
                if (s2 > max_width) {
                    max_width = s2;
                }
            }
        }
    }
    this->col_max_width_ = max_width;
    //initialize the vector col_widths_ by all the same number max_width
    for (int icol = 0; icol < this->ncol_; ++icol) {
        if (!flexible_width_) {
            this->col_widths_.push_back(max_width);
            this->total_width_ += max_width;
        }
        else {
            this->total_width_ += this->col_widths_[icol];
            if (this->col_max_width_ < this->col_widths_[icol]) {
                this->col_max_width_ = this->col_widths_[icol];
            }
        }
    }
    this->total_width_ += this->ncol_+1; // add the width of the delimiters
}

void formatter::Table::centerize_title() {
    for (int icol = 0; icol < this->ncol_; ++icol) {
        int s = this->titles_[icol].size();
        int n = this->col_widths_[icol];
        int n1 = (n-s)/2;
        int n2 = n-s-n1;
        std::string title = "";
        for (int i = 0; i < n1; ++i) {
            title += " ";
        }
        title += this->titles_[icol];
        for (int i = 0; i < n2; ++i) {
            title += " ";
        }
        this->titles_[icol] = title;
    }
}

void formatter::Table::reset() {
    this->mode_ = 0;

    this->overall_title = "";
    
    this->ncol_ = 0;
    this->col_delimiter_ = ' ';
    this->frame_switches_ = {1, 1, 0, 0};
    this->frame_delimiters_ = {'-', '-', '|', '|'};
    this->frame_mid_switch_ = 1;
    this->frame_mid_delimiter_ = '-';
    this->flexible_width_ = true;
    this->col_max_width_ = 0;
    this->total_width_ = 0;

    this->title_position_ = 0;
    this->clean();
}

void formatter::Table::clean() {
    this->titles_.clear();
    this->col_widths_.clear();
    this->data_.clear();
}

std::string formatter::Table::print_table() {

    this->adjust_col_width();
    if (this->title_position_ == 0) {
        this->centerize_title();
    }
    std::stringstream ss;

    // the numbers of columns are allowed to be different, get the maximum number of rows
    int nrow_max = this->data_[0].size();
    std::vector<int> nrows;
    int ncol = this->ncol_;
    for (int icol = 0; icol < ncol; ++icol) {
        int s = this->data_[icol].size();
        nrows.push_back(s);
        if (nrow_max < s) {
            nrow_max = s;
        }
    }
    if (this->mode_ != 1) {
        if (this->overall_title.size() > 0) {
            ss << this->overall_title << std::endl;
        }
        // print the top frame
        if (this->frame_switches_[0]) {
            for (int iw = 0; iw < this->total_width_; ++iw) {
                ss << this->frame_delimiters_[0];
            }
            ss << std::endl;
        }
        // print the title
        if (this->frame_switches_[2]) {
            ss << this->frame_delimiters_[2];
        }
        for (int icol = 0; icol < ncol; ++icol) {
            ss << std::setw(this->col_widths_[icol]) 
            << std::setfill(' ');
            if(this->title_position_ <= 0) {
                ss << std::left;
            }
            else {
                ss << std::right;
            }
            ss << this->titles_[icol];
            if (icol != ncol-1) {
                ss << this->col_delimiter_;
            }
        }
        if (this->frame_switches_[3]) {
            ss << this->frame_delimiters_[3];
        }
        ss << std::endl;
        // print the middle frame
        if (this->frame_mid_switch_) {
            for (int iw = 0; iw < this->total_width_; ++iw) {
                ss << this->frame_mid_delimiter_;
            }
            ss << std::endl;
        }
    }
    // print the data
    if (this->mode_ >= 0) {
        for (int irow = 0; irow < nrow_max; ++irow) {
            if (this->frame_switches_[2]) {
                ss << this->frame_delimiters_[2];
            }
            for (int icol = 0; icol < ncol; ++icol) {
                if (irow < nrows[icol]) {
                    ss << std::setw(this->col_widths_[icol]) 
                    << std::setfill(' ') 
                    << std::left 
                    << this->data_[icol][irow];
                }
                else {
                    ss << std::setw(this->col_widths_[icol]) 
                    << std::setfill(' ') 
                    << std::left 
                    << "";
                }
                if (icol != ncol-1) {
                    ss << this->col_delimiter_;
                }
            }
            if (this->frame_switches_[3]) {
                ss << this->frame_delimiters_[3];
            }
            ss << std::endl;
        }
    }
    // finally the bottom frame
    if (this->mode_ == 0) {
        if (this->frame_switches_[1]) {
            for (int iw = 0; iw < this->total_width_; ++iw) {
                ss << this->frame_delimiters_[1];
            }
            ss << std::endl;
        }
    }
    this->reset();
    return ss.str();
}

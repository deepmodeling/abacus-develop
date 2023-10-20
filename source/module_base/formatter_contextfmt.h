#ifndef FORMATTER_CONTEXTFMT_H
#define FORMATTER_CONTEXTFMT_H

#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <utility>
#include "formatter_fmt.h"
#include "formatter_physfmt.h"
#include "formatter_table.h"

namespace formatter
{
    class ContextFmt : public Table
    {
        public:
            /// @brief default constructor, for default values see declarations
            ContextFmt();
            /// @brief destructor
            ~ContextFmt();

            // physical formatter default context: energy
            /// @brief set context by context name, if it is already a predefined type, then directly initialize phys_fmt_ vector, set number of columns and title.
            /// @param context context name, see this->predefined_phys_fmt
            /// @note this function is mostly used for inputing data in 0 dimension one-by-one
            /// @attention for pointer pass to << operator, you MUST NOT USE THIS but to provide another int* nrows!
            void set_context(std::string context); // most simple case, rely on default values pre-defined
            /// @brief (pointer overloaded version) set context by context name, , if it is already a predefined type, then directly initialize phys_fmt_ vector, set number of columns and title.
            /// @param context context name, see this->predefined_phys_fmt
            /// @param ncol number of columns, needed to tranverse the last param nrows
            /// @param nrows (pointer case specifically compulsory) rows of every column
            void set_context(std::string context, int ncol, int* &nrows);
            /// @brief (pointer overloaded version) set context directly by properties of each column stored in std::string*, also set number of rows for each column
            /// @param ncol number of columns
            /// @param phys_fmt physical format of each column
            /// @param nrows number of rows for each column
            void set_context(int ncol, std::string* phys_fmt, int* nrows); // general case, for array input
            /// @brief (STL vector overloaded version) set context directly by properties of each column stored in std::vector<std::string>, the number of rows for each column is automatically contained in input data vector.
            /// @param phys_fmt physical format of each column
            void set_context(std::vector<std::string> phys_fmt); // general case, for vector input
            
            /// @brief setter the default physical context, for cases the context is not pre defined, like vector3d, vector3d_i, etc.
            /// @param default_phys_fmt default physical context, default is "energy" now
            void set_default_phys_fmt(std::string default_phys_fmt) { default_phys_fmt_ = default_phys_fmt; }
            /// @brief setter the cache title, to be used for adding new column to table
            /// @param cache_title cache title
            void set_cache_title(std::string cache_title) { cache_title_ = cache_title; }
            /// @brief setter the number of columns, only used for pointer input data
            /// @param ncol number of columns
            void set_ncol(int ncol) { ncol_ = ncol; }
            /// @brief setter the current column index
            /// @param icol current column index
            void set_icol(int icol) { icol_ = icol; }

            /// @brief setter the title_switch
            /// @param title_switch title switch, 0: not print title, 1: print title
            void set_title_switch(int title_switch) { title_switch_ = title_switch; }

            /// @brief overloaded operator<<, for inputing data in 0 dimension one-by-one
            /// @tparam T int, double, std::string
            /// @param value input value
            /// @return *this
            /// @attention if enable_title(), you can input title and then the data one-by-one, but if disable_title(), you can only input data one-by-one, or unexpected result will be produced
            template <typename T>
            formatter::ContextFmt& operator<<(const T& value) {
                if (this->title_switch_%2 == 0) {
                    // it is a title
                    this->cache_title_ = std::to_string(value);
                }
                else {
                    // it is a data
                    if (this->known_context_) {
                        this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
                    }
                    else {
                        this->p_phys_fmt_->set_context(this->default_phys_fmt_);
                    }
                    Table::add_col(this->cache_title_, (std::vector<std::string>){this->fmt_.format(value)});
                    this->cache_title_ = "";
                    this->fmt_.reset();
                    this->icol_++;
                }
                if (Table::get_mode() == 1) {
                    this->title_switch_ += 2;
                }
                else {
                    this->title_switch_++;
                }
                return *this;
            }
            /// @brief overloaded operator<<, for inputing data in 1 dimension std::vector<T>
            /// @tparam T int, double and std::string
            /// @param value input value
            /// @return *this
            /// @attention if enable_title(), you can input title and then the data one-by-one, but if disable_title(), you can only input data one-by-one, or unexpected result will be produced
            template <typename T>
            formatter::ContextFmt& operator<<(const std::vector<T>& value) {
                // well if a whole vector is input, it is definitely a data... or presently I can't think of a case that a vector is a title
                std::vector<std::string> value_;
                if (this->known_context_) {
                    this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
                }
                else if (this->icol_ < this->phys_fmt_.size()) {
                    this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
                }
                else {
                    this->p_phys_fmt_->set_context(this->default_phys_fmt_);
                }

                for (auto v : value) {
                    value_.push_back(this->fmt_.format(v));
                }
                Table::add_col(this->cache_title_, value_);
                this->cache_title_ = "";
                this->fmt_.reset();
                this->icol_++;
                this->title_switch_++;
                return *this;
            }
            /// @brief overloaded operator<<, for inputing data in 1 dimension T*
            /// @tparam T int, double and std::string
            /// @param value input value
            /// @return *this
            /// @attention if enable_title(), you can input title and then the data one-by-one, but if disable_title(), you can only input data one-by-one, or unexpected result will be produced
            template <typename T>
            formatter::ContextFmt& operator<<(T*& value) {
                std::vector<std::string> value_;
                if (this->known_context_) {
                    this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
                }
                else if (this->icol_ < this->phys_fmt_.size()) {
                    this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
                }
                else {
                    this->p_phys_fmt_->set_context(this->default_phys_fmt_);
                }
                for (int i = 0; i < this->nrows_[this->icol_]; ++i) {
                    value_.push_back(this->fmt_.format(value[i]));
                }
                Table::add_col(this->cache_title_, value_);
                this->cache_title_ = "";
                this->fmt_.reset();
                this->icol_++;
                this->title_switch_++;
                return *this;
            }

            /// @brief overloaded operator<<, for inputing data in char*
            /// @param value input value
            /// @return *this
            formatter::ContextFmt& operator<<(char const* value);

            /// @brief reset all
            void reset();
            void context_refresh();
            /// @brief print the table
            std::string str(bool next_line = false);
            /// @brief set to mode in which title will be manually input and certainly will be output
            void enable_title() { Table::set_mode(0); this->set_title_switch(0); }
            /// @brief set to mode in which title will not be manually input and certainly will not be output
            void disable_title() { Table::set_mode(1); this->set_title_switch(1); }
            /// @brief actually it is something strange and the function is not recommended to use
            void only_title() { Table::set_mode(-1); this->set_title_switch(0); }
            /// @brief enable the iterative mode, in which the data will be output in iterative way, for details see online document whose link is provided in the header of this file
            void enable_iterative() { iterative_ = 1; }
            /// @brief disable the iterative mode
            void disable_iterative() { iterative_ = 0; }

            /// @brief getter the context
            /// @return context
            std::string get_context() const { return context_; }
            /// @brief getter the physical format of each column
            /// @return std::vector<std::string> of physical format
            std::vector<std::string> get_phys_fmt() const { return phys_fmt_; }
            /// @brief getter the pointer to the physical formatter
            /// @return pointer to the physical formatter
            formatter::PhysicalFmt* get_p_phys_fmt() { return p_phys_fmt_; }
            /// @brief getter the default physical context
            /// @return default physical context
            std::string get_default_phys_fmt() const { return default_phys_fmt_; }
            /// @brief getter the cache title
            /// @return cache title
            std::string get_cache_title() const { return cache_title_; }
            /// @brief getter the number of columns, only used for pointer input data
            /// @return number of columns
            int get_ncol() const { return ncol_; }
            /// @brief getter the current column index
            /// @return current column index
            int get_icol() const { return icol_; }
            /// @brief getter the number of rows for each column
            /// @return std::vector<int> of number of rows
            std::vector<int> get_nrows() const { return nrows_; }
            /// @brief getter the formatter
            /// @return formatter
            formatter::Fmt& get_fmt() { return fmt_; }
            /// @brief getter the stringstream (did I use it? I don't remember)
            /// @return stringstream
            std::stringstream& get_ss() { return ss_; }
            /// @brief getter the title switch
            /// @return title switch
            int get_title_switch() const { return title_switch_; }
            /// @brief getter the iterative mode
            /// @return iterative mode
            bool get_iterative() const { return iterative_; }
            // for debug, no need to test
            void print_status() const;
        private:
            /// @brief context
            std::string context_;
            /// @brief physical format of each column
            std::vector<std::string> phys_fmt_;
            /// @brief default physical context
            std::string default_phys_fmt_ = "energy";
            /// @brief cache title, to be used for adding new column to table
            std::string cache_title_ = "";
            /// @brief whether the title is manually input
            bool with_title_ = true;

            /// @brief number of columns, only used for pointer input data
            int ncol_ = 0;
            /// @brief current column index
            int icol_ = 0;
            /// @brief number of rows for each column
            std::vector<int> nrows_; // for operator<< the T* input case
            /// @brief boolean controlling whether the iterative mode is enabled
            int iterative_ = 0;

            /// @brief formatter
            formatter::Fmt fmt_;
            /// @brief pointer to the physical formatter
            formatter::PhysicalFmt* p_phys_fmt_;
            /// @brief stringstream
            std::stringstream ss_;
            /// @brief title switch, 0: not print title, 1: print title
            int title_switch_ = 0;
            // values of title_switch_:
            // title: 0 2 4 6 8 ... , for mode = 1, stepsize = 2, otherwise 1
            // data:  1 3 5 7 9 ...

            /// @brief boolean controlling whether the context is known
            bool known_context_ = false;
            /// @brief database of predefined physical format
            std::unordered_map<std::string, std::pair<int, std::vector<std::string>>> predefined_phys_fmt = {
                {"vector3d", std::make_pair(1,
                std::vector<std::string>{"coordinate", "coordinate", "coordinate"})}, // one can define general context (recommended)
                {"vector3d_i", std::make_pair(1,
                std::vector<std::string>{"constraint", "constraint", "constraint"})}, // vector3d will be position, vectors and for this, it is constraint, kmesh, ...
                {"scf", std::make_pair(0,
                std::vector<std::string>{"str_w4", "int_w4", "energy", "energy", "energy", "time"})}, // but for scf it is really a special case
                {"time_statistics", std::make_pair(0,
                std::vector<std::string>{"mid_title", "mid_title", "time", "int_w8", "time", "percentage"})}, // so is time statistics
                {"atomic_species", std::make_pair(1,
                std::vector<std::string>{"str_w4", "mass", "str_w30"})}, // so is ATOMIC_SPECIES
                {"lattice_vectors", std::make_pair(1,
                std::vector<std::string>{"coordinate", "coordinate", "coordinate"})} // but it is not true for LATTICE_VECTORS
            };
    };
}

#endif // FORMATTER_CONTEXTFMT_H
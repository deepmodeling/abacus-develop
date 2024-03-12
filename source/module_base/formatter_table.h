#ifndef FORMATTER_TABLE_H
#define FORMATTER_TABLE_H

#include <string>
#include <vector>

namespace formatter
{
    class Table 
    {
        public:
            /// @brief default constructor
            Table(){};
            /// @brief destructor
            ~Table(){};
            /// @brief setter the mode
            /// @param new_mode 0: table, 1: data, -1: header
            void set_mode(int new_mode) { mode_ = new_mode; }
            /// @brief setter the overall title
            /// @param new_title overall title
            void set_overall_title(std::string new_title) { overall_title = new_title; }
            /// @brief setter the titles of each column
            /// @param new_titles std::vector<std::string> of titles
            void set_titles(std::vector<std::string> new_titles) { titles_ = new_titles; }
            /// @brief setter the title for one specific column
            /// @param icol column index
            /// @param new_title new title
            void set_title(int icol, std::string new_title) { titles_[icol] = new_title; }
            /// @brief setter the column delimiter
            /// @param new_delimiter new delimiter
            void set_col_delimiter(char new_delimiter) { col_delimiter_ = new_delimiter; }
            /// @brief setter the frame switches
            /// @param new_switches std::vector<int> of switches, 1: print frame, 0: not print frame
            void set_frame_switches(std::vector<int> new_switches) { frame_switches_ = new_switches; }
            /// @brief setter the frame switch for one specific frame
            /// @param frame frame index, 0: top, 1: down, 2: left, 3: right
            /// @param new_switch new switch, 1: print frame, 0: not print frame
            void set_frame_switch(int frame, int new_switch) { frame_switches_[frame] = new_switch; }
            /// @brief setter the frame delimiters
            /// @param new_delimiters std::vector<char> of delimiters
            void set_frame_delimiters(std::vector<char> new_delimiters) { frame_delimiters_ = new_delimiters; }
            /// @brief setter the frame delimiter for one specific frame
            /// @param iframe frame index, 0: top, 1: down, 2: left, 3: right
            /// @param new_delimiter new delimiter
            void set_frame_delimiter(int iframe, char new_delimiter) { frame_delimiters_[iframe] = new_delimiter; }
            /// @brief setter the frame mid switch
            /// @param new_switch new switch, 1: print frame, 0: not print frame
            void set_frame_mid_switch(int new_switch) { frame_mid_switch_ = new_switch; }
            /// @brief setter the frame mid delimiter
            /// @param new_delimiter new delimiter
            void set_frame_mid_delimiter(char new_delimiter) { frame_mid_delimiter_ = new_delimiter; }

            /// @brief enable the flexible width, if set true, for each column different width will be used
            void enable_flexible_width() { flexible_width_ = true; }
            /// @brief disable the flexible width
            void disable_flexible_width() { flexible_width_ = false; }
            /// @brief print frame for top
            void enable_up_frame() { frame_switches_[0] = 1; }
            /// @brief not print frame for top
            void disable_up_frame() { frame_switches_[0] = 0; }
            /// @brief print frame for down
            void enable_down_frame() { frame_switches_[1] = 1; }
            /// @brief not print frame for down
            void disable_down_frame() { frame_switches_[1] = 0; }
            /// @brief print frame for left
            void enable_left_frame() { frame_switches_[2] = 1; }
            /// @brief not print frame for left
            void disable_left_frame() { frame_switches_[2] = 0; }
            /// @brief print frame for right
            void enable_right_frame() { frame_switches_[3] = 1; }
            /// @brief not print frame for right
            void disable_right_frame() { frame_switches_[3] = 0; }
            /// @brief print frame between title and data
            void enable_mid_frame() { frame_mid_switch_ = 1; }
            /// @brief not print frame between title and data
            void disable_mid_frame() { frame_mid_switch_ = 0; }

            /// @brief centerize the title
            void center_title() { title_position_ = 0; }
            /// @brief make title left
            void left_title() { title_position_ = -1; }
            /// @brief make title right
            void right_title() { title_position_ = 1; }
            /// @brief add new column of data to the table
            /// @param new_title title of the new column
            /// @param new_col data, stored in std::vector<std::string>
            /// @return return the index of the new column
            int add_col(std::string new_title, std::vector<std::string> new_col);
            int add_cols(std::vector<std::string> new_titles, std::vector<std::vector<std::string>> new_cols); // later
            void del_col(int col);  // not implemented
            void del_col(std::string title); // not implemented
            void del_cols(std::vector<std::string> titles); // not implemented

            /// @brief adjust the width of each column according to the data and title
            void adjust_col_width();
            /// @brief DO NOT CALL IT, it is called by adjust_col_width()
            void centerize_title();
            /// @brief clean the data
            void clean();
            /// @brief reset the table
            void reset();
            /// @brief print the table
            std::string print_table();

            /// @brief getter the mode
            /// @return mode, 0: table, 1: data, -1: header
            int get_mode() const { return mode_; }
            /// @brief getter the overall title
            /// @return overall title
            std::string get_overall_title() const { return overall_title; }
            /// @brief getter the number of columns
            /// @return number of columns
            int get_ncol() const { return ncol_; } 
            /// @brief getter the titles of each column
            /// @return std::vector<std::string> of titles
            std::vector<std::string> get_titles() const { return titles_; }
            /// @brief getter the delimiter between columns
            /// @return delimiter
            char get_col_delimiter() const { return col_delimiter_; }
            /// @brief getter the frame switches
            /// @return std::vector<int> of switches, 1: print frame, 0: not print frame
            std::vector<int> get_frame_switches() const { return frame_switches_; }
            /// @brief getter the frame delimiters
            /// @return std::vector<char> of delimiters
            std::vector<char> get_frame_delimiters() const { return frame_delimiters_; }
            /// @brief getter the frame mid switch
            /// @return frame mid switch, 1: print frame, 0: not print frame
            int get_frame_mid_switch() const { return frame_mid_switch_; }
            /// @brief getter the frame mid delimiter
            /// @return frame mid delimiter
            char get_frame_mid_delimiter() const { return frame_mid_delimiter_; }
            /// @brief getter the boolean controlling flexible width switch
            /// @return flexible width switch
            bool get_flexible_width() const {return flexible_width_;}
            /// @brief getter the width of each column
            /// @return std::vector<int> of width
            std::vector<int> get_col_widths() const { return col_widths_; }
            /// @brief getter the maximum width of each column
            /// @return maximum width
            int get_col_max_width() const { return col_max_width_; }
            /// @brief getter the total width of the table, will be the sum of column width and delimiter width plus 3 (delimiters)
            /// @return total width
            int get_total_width() const { return total_width_; }
            /// @brief getter the data of the table
            /// @return std::vector<std::vector<std::string>> of data
            std::vector<std::vector<std::string>> get_data() const { return data_; }
        private:
            int mode_ = 0; // 0: table, 1: data, -1: header

            /// @brief overall title
            std::string overall_title;

            /// @brief number of columns
            int ncol_ = 0;
            /// @brief titles of each column
            std::vector<std::string> titles_;
            /// @brief delimiter between columns
            char col_delimiter_ = ' ';
            /// @brief frame switches, 1: print frame, 0: not print frame, positional info.: 0: top, 1: down, 2: left, 3: right
            std::vector<int> frame_switches_ = {1, 1, 0, 0};
            /// @brief frame delimiters, positional info.: 0: top, 1: down, 2: left, 3: right
            std::vector<char> frame_delimiters_ = {'-', '-', '|', '|'};
            /// @brief frame mid switch, the mid frame means the frame between title and data, 1: print frame, 0: not print frame
            int frame_mid_switch_ = 1;
            /// @brief frame mid delimiter
            char frame_mid_delimiter_ = '-';
            /// @brief if set true, for each column different width will be used
            bool flexible_width_ = true;
            /// @brief width of each column
            std::vector<int> col_widths_;
            /// @brief maximum width of each column
            int col_max_width_ = 0;
            /// @brief total width of the table
            int total_width_ = 0;
            /// @brief data of the table
            int title_position_ = -1;
            // warning: the col data is stored here row-by-row
            std::vector<std::vector<std::string>> data_;
    };
}
#endif
#ifndef FORMATTER_FMT_H
#define FORMATTER_FMT_H

#include <string>
#include <sstream>
#include <iomanip>
#include <complex>

namespace formatter
{
    class Fmt {
            
        public:
            /// @brief default constructor, set default values: width = 4, precision = 2, fillChar = ' ', fixed = true, right = true, error = false
            Fmt();
            /// @brief constructor with input values
            /// @param width width of the output string
            /// @param precision precision of the output string
            /// @param fillChar fill character of the output string
            /// @param fixed whether the output string keeps decimal places
            /// @param right whether the output string is right aligned
            /// @param error whether the output string has "+" sign if positive
            Fmt(int width, int precision, char fillChar, bool fixed, bool right, bool error);
            /// @brief destructor
            ~Fmt();
            /// @brief setter the width
            /// @param width 
            void set_width(int width) { width_ = width; }
            /// @brief setter the precision
            /// @param precision
            void set_precision(int precision) { precision_ = precision; }
            /// @brief setter the fillChar
            /// @param fillChar
            void set_fillChar(char fillChar) { fillChar_ = fillChar; }
            /// @brief setter the fixed
            /// @param fixed
            void set_fixed(bool fixed) { fixed_ = fixed; }
            /// @brief setter the right
            /// @param right
            void set_right(bool right) { right_ = right; }
            /// @brief setter the error
            /// @param error
            void set_error(bool error) { error_ = error; }
            /// @brief reset the format to default values
            void reset();
            /// @brief format the input value to string
            /// @tparam T double, int or std::string
            /// @param value input value
            /// @return std::string, the formatted string
            template <typename T> std::string format(const T& value) {
                std::stringstream ss_;
                if (std::is_signed<T>::value) {
                    if (value >= 0) {
                        if (error_) ss_ << "+";
                    } else {}
                }
                ss_ << std::setprecision(precision_);
                std::stringstream ss;
                if (fixed_) {
                    ss_ << std::fixed;
                } else {
                    ss_ << std::scientific;
                }
                ss_ << (double)value;

                if (!right_) ss << std::left;
                ss << std::setw(width_) << std::setfill(fillChar_);
                ss << ss_.str();
                return ss.str();
            }
            std::string format(const std::complex<double>& value){
                return "(" + format(value.real()) + "," + format(value.imag()) + ")";
            }
            std::string format(const std::complex<float>& value){
                return "(" + format(value.real()) + "," + format(value.imag()) + ")";
            }
            std::string format(const std::string& value) {
                std::stringstream ss;
                if (!right_) ss << std::left;
                ss << std::setw(width_) << std::setfill(fillChar_);
                ss << value;
                return ss.str();
            }
            /// @brief getter the width
            /// @return width
            int get_width() const { return width_; }
            /// @brief getter the precision
            /// @return precision
            int get_precision() const { return precision_; }
            /// @brief getter the fillChar
            /// @return fillChar
            char get_fillChar() const { return fillChar_; }
            /// @brief getter the fixed
            /// @return fixed
            bool get_fixed() const { return fixed_; }
            /// @brief getter the right
            /// @return right
            bool get_right() const { return right_; }
            /// @brief getter the error
            /// @return error
            bool get_error() const { return error_; }

        private:
            /// @brief width of the output string
            int width_ = 4;
            /// @brief precision of the output string
            int precision_ = 2;
            /// @brief fill character of the output string
            char fillChar_ = ' ';
            /// @brief whether the output string keeps decimal places
            bool fixed_ = true;
            /// @brief whether the output string is right aligned
            bool right_ = true;
            /// @brief whether the output string has "+" sign if positive
            bool error_ = false;
    };
}
#endif // FORMATTER_FMT_H
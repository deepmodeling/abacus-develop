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

            /// @brief format the input value to string
            /// @tparam T double, int or std::string
            /// @param value input value
            /// @return std::string, the formatted string
            template <typename T> std::string format(const T& value) {
                std::stringstream ss_;
                if((std::is_same<T, std::string>::value)||(std::is_same<T, const char*>::value)
                    ||(std::is_same<T, char*>::value)||(std::is_same<T, const char>::value)||(std::is_same<T, char>::value))
                {
                    // if is string, no need to set precision
                }
                else
                {
                    // if is double, float, int, can set error flag
                    if (error_) {
                        if (value >= 0) ss_ << "+";
                    }
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

            template <typename T> std::string format(const std::string& strfmt, const T& value) {
                int width, precision;
                bool fixed, right;
                to_format(strfmt, width, precision, fixed, right);
                set_width(width);
                set_precision(precision);
                set_fixed(fixed);
                set_right(right);
                return format<T>(value);
            }
            /// @brief for compatibility with printf fashion format string
            /// @param strfmt [in] "%20.10f"-like format string
            /// @param width [out] width of the output string
            /// @param precision [out] precision of the output string
            /// @param fixed [out] whether the output string keeps decimal places
            /// @param right [out] whether the output string is right aligned
            void to_format(const std::string& strfmt,
                        int& width,
                        int& precision,
                        bool& fixed,
                        bool& right);

            // getter and setter
            int width() const { return width_; }
            int precision() const { return precision_; }
            char fillChar() const { return fillChar_; }
            bool fixed() const { return fixed_; }
            bool right() const { return right_; }
            bool error() const { return error_; }

            void set_width(int width) { width_ = width; }
            void set_precision(int precision) { precision_ = precision; }
            void set_fillChar(char fillChar) { fillChar_ = fillChar; }
            void set_fixed(bool fixed) { fixed_ = fixed; }
            void set_right(bool right) { right_ = right; }
            void set_error(bool error) { error_ = error; }
            /// @brief reset the format to default values
            void reset();

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
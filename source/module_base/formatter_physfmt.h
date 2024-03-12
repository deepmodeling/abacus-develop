#ifndef FORMATTER_PHYSFMT_H
#define FORMATTER_PHYSFMT_H

#include <string>

#include "formatter_fmt.h"

namespace formatter
{
    class PhysicalFmt {
        public:
            /// @brief default constructor, set default values: context = "none", p_formatter = nullptr
            /// @param context available choices see function adjust_formatter()
            /// @param p_formatter pointer to the formatter, can be passed with pointer to Fmt object or nullptr
            PhysicalFmt(std::string context = "none", Fmt* p_formatter = nullptr);
            /// @brief destructor
            ~PhysicalFmt();

            /// @brief adjust the formatter according to the context
            /// @param left whether the output string is left aligned
            void adjust_formatter(bool left = false);
            /// @brief adjust formatter but provide interface to external for more flexibly control precision
            /// @param decisive_length for double, it is precision, for int, it is width
            /// @param width_ratio for double specifically, the length ratio len(integer)/len(decimal) part. If a negative value is given, then there will not be alignment.
            /// @param fillchar most of time it is okay to keep it as default, set as like otherwise
            /// @param left whether the output string is left aligned
            void adjust_formatter_flexible(const int& decisive_length, 
                                           const double& width_ratio = 0.0, 
                                           const bool& scientific = false,
                                           const bool& left = false,
                                           const char& fillchar = ' '
                                           );
            /// @brief setter the context
            /// @param context available choices see function adjust_formatter()
            void set_context(std::string context);
            /// @brief setter the pointer to the formatter
            /// @param p_formatter pointer to the formatter, can be passed with pointer to Fmt object or nullptr
            void set_p_formatter(Fmt* p_formatter);
            /// @brief setter the decorator mode, only for debug and unittest purpose
            /// @param decorator_mode if external Fmt object is set to bundle with PhysicalFmt object, then decorator_mode should be set true, otherwise false
            void set_decorator_mode(bool decorator_mode) { decorator_mode_ = decorator_mode; }
            /// @brief getter the context
            /// @return context
            std::string get_context() const { return context_; }
            /// @brief getter the pointer to the formatter
            /// @return pointer to the formatter
            Fmt* get_p_formatter() const { return p_formatter_; }
            /// @brief getter the decorator mode, only for debug and unittest purpose
            /// @return decorator_mode
            bool get_decorator_mode() const { return decorator_mode_; }
        private:
            /// @brief context, available choices see function adjust_formatter()
            std::string context_ = "none";
            /// @brief pointer to the formatter, can be passed with pointer to Fmt object or nullptr
            Fmt* p_formatter_;

            /// @brief decorator mode, indicating whether the external Fmt object is set to bundle with PhysicalFmt object
            bool decorator_mode_ = false;
    };
}
#endif 
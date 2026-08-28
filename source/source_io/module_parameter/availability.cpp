#include "source_io/module_parameter/availability.h"

#include <cctype>
#include <stdexcept>
#include <utility>

namespace ModuleIO
{
namespace
{

bool is_identifier_start(const char c)
{
    return std::isalpha(static_cast<unsigned char>(c)) != 0;
}

bool is_identifier_char(const char c)
{
    return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_';
}

std::string format_value(const std::string& value)
{
    for (const char c : value)
    {
        if (std::isspace(static_cast<unsigned char>(c)) || c == '(' || c == ')'
            || c == '[' || c == ']' || c == ',' || c == '=' || c == '<'
            || c == '>' || c == '!')
        {
            return "\"" + value + "\"";
        }
    }
    return value;
}

class AvailabilityParser
{
  public:
    explicit AvailabilityParser(const std::string& text) : text_(text) {}

    AvailabilityExpr parse()
    {
        if (at_end())
        {
            return AvailabilityExpr();
        }
        skip_space();
        if (at_end())
        {
            fail("expected parameter identifier");
        }

        AvailabilityExpr result = parse_or();
        skip_space();
        if (!at_end())
        {
            fail("unexpected token");
        }
        return result;
    }

  private:
    const std::string& text_;
    std::size_t pos_ = 0;

    bool at_end() const
    {
        return pos_ == text_.size();
    }

    void skip_space()
    {
        while (!at_end() && std::isspace(static_cast<unsigned char>(text_[pos_])))
        {
            ++pos_;
        }
    }

    [[noreturn]] void fail(const std::string& message) const
    {
        throw std::invalid_argument("Invalid availability expression at column "
                                    + std::to_string(pos_ + 1) + ": " + message
                                    + " in '" + text_ + "'");
    }

    bool consume_char(const char expected)
    {
        skip_space();
        if (!at_end() && text_[pos_] == expected)
        {
            ++pos_;
            return true;
        }
        return false;
    }

    bool consume_keyword(const char* keyword)
    {
        skip_space();
        const std::size_t length = std::char_traits<char>::length(keyword);
        if (text_.compare(pos_, length, keyword) != 0)
        {
            return false;
        }
        const std::size_t end = pos_ + length;
        if ((pos_ > 0 && is_identifier_char(text_[pos_ - 1]))
            || (end < text_.size() && is_identifier_char(text_[end])))
        {
            return false;
        }
        pos_ = end;
        return true;
    }

    std::string parse_identifier()
    {
        skip_space();
        if (at_end() || !is_identifier_start(text_[pos_]))
        {
            fail("expected parameter identifier");
        }
        const std::size_t begin = pos_++;
        while (!at_end() && is_identifier_char(text_[pos_]))
        {
            ++pos_;
        }
        return text_.substr(begin, pos_ - begin);
    }

    std::string parse_scalar_value()
    {
        skip_space();
        if (!at_end() && text_[pos_] == '"')
        {
            ++pos_;
            const std::size_t content_begin = pos_;
            while (!at_end() && text_[pos_] != '"')
            {
                ++pos_;
            }
            if (at_end())
            {
                fail("expected closing quote");
            }
            if (pos_ == content_begin)
            {
                fail("expected non-empty quoted value");
            }
            const std::string value = text_.substr(content_begin, pos_ - content_begin);
            ++pos_;
            return value;
        }

        const std::size_t begin = pos_;
        while (!at_end() && !std::isspace(static_cast<unsigned char>(text_[pos_]))
               && text_[pos_] != '(' && text_[pos_] != ')' && text_[pos_] != '['
               && text_[pos_] != ']' && text_[pos_] != ',')
        {
            if (text_[pos_] == '=' || text_[pos_] == '<' || text_[pos_] == '>'
                || text_[pos_] == '!')
            {
                fail("invalid character in value");
            }
            ++pos_;
        }
        if (begin == pos_)
        {
            fail("expected value");
        }
        return text_.substr(begin, pos_ - begin);
    }

    AvailabilityExpr parse_condition()
    {
        AvailabilityExpr result;
        skip_space();
        if (at_end())
        {
            fail("expected parameter identifier");
        }

        result.condition.param = parse_identifier();
        skip_space();

        static const char* comparison_ops[] = {"==", "!=", ">=", "<=", ">", "<"};
        for (const char* op : comparison_ops)
        {
            const std::size_t length = std::char_traits<char>::length(op);
            if (text_.compare(pos_, length, op) == 0)
            {
                pos_ += length;
                result.condition.op = op;
                result.condition.values.push_back(parse_scalar_value());
                return result;
            }
        }

        if (consume_keyword("contains"))
        {
            result.condition.op = "contains";
            result.condition.values.push_back(parse_scalar_value());
            return result;
        }

        if (consume_keyword("in"))
        {
            if (!consume_char('['))
            {
                fail("expected '[' after 'in'");
            }
            result.condition.op = "in";
            result.condition.values.push_back(parse_scalar_value());
            while (consume_char(','))
            {
                result.condition.values.push_back(parse_scalar_value());
            }
            if (!consume_char(']'))
            {
                fail("expected ']' after list");
            }
            if (result.condition.values.size() == 1)
            {
                fail("use '==' for a single value");
            }
            return result;
        }

        fail("expected comparison operator");
    }

    AvailabilityExpr parse_primary()
    {
        skip_space();
        if (consume_char('('))
        {
            AvailabilityExpr result = parse_or();
            if (!consume_char(')'))
            {
                fail("expected ')' after grouped expression");
            }
            return result;
        }
        return parse_condition();
    }

    AvailabilityExpr parse_and()
    {
        AvailabilityExpr result = parse_primary();
        std::vector<AvailabilityExpr> children;
        children.push_back(result);
        while (consume_keyword("and"))
        {
            children.push_back(parse_primary());
        }
        if (children.size() == 1)
        {
            return result;
        }
        AvailabilityExpr node;
        node.op = "and";
        node.children = std::move(children);
        return node;
    }

    AvailabilityExpr parse_or()
    {
        AvailabilityExpr result = parse_and();
        std::vector<AvailabilityExpr> children;
        children.push_back(result);
        while (consume_keyword("or"))
        {
            children.push_back(parse_and());
        }
        if (children.size() == 1)
        {
            return result;
        }
        AvailabilityExpr node;
        node.op = "or";
        node.children = std::move(children);
        return node;
    }
};

} // namespace

std::string AvailabilityCondition::to_string() const
{
    if (values.empty())
    {
        return {};
    }
    if (op == "in")
    {
        std::string result = param + " in [";
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            if (i)
            {
                result += ", ";
            }
            result += format_value(values[i]);
        }
        return result + "]";
    }
    if (op == "contains")
    {
        return param + " contains " + format_value(values[0]);
    }
    return param + op + format_value(values[0]);
}

std::string AvailabilityExpr::to_string() const
{
    if (is_leaf())
    {
        return condition.to_string();
    }
    const std::string separator = op == "or" ? " or " : " and ";
    std::string result;
    for (std::size_t i = 0; i < children.size(); ++i)
    {
        if (i)
        {
            result += separator;
        }
        if (children[i].is_leaf())
        {
            result += children[i].to_string();
        }
        else
        {
            result += "(" + children[i].to_string() + ")";
        }
    }
    return result;
}

AvailabilityExpr parse_availability(const std::string& raw)
{
    AvailabilityExpr expression = AvailabilityParser(raw).parse();
    const std::string canonical = expression.to_string();
    if (raw != canonical)
    {
        throw std::invalid_argument("Non-canonical availability expression '" + raw
                                    + "'; expected '" + canonical + "'");
    }
    return expression;
}

} // namespace ModuleIO

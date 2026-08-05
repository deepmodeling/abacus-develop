#include "source_io/module_parameter/availability.h"

#include <cctype>

namespace ModuleIO
{
namespace
{

std::string trim_copy(const std::string& s)
{
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b])))
    {
        ++b;
    }
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
    {
        --e;
    }
    return s.substr(b, e - b);
}

std::string lower_copy(const std::string& s)
{
    std::string r = s;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        r[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(r[i])));
    }
    return r;
}

/// Split `text` on the keyword `kw` (e.g. "and"/"or", space separated, or the
/// plain "," separator) occurring at bracket/paren depth 0. Respects () and []
/// so value lists are not split. Returns trimmed, non-empty parts.
std::vector<std::string> split_top_keyword(const std::string& text,
                                           const std::string& kw)
{
    std::vector<std::string> out;
    int depth = 0;
    std::size_t start = 0;
    const std::size_t n = text.size();
    std::size_t i = 0;
    while (i < n)
    {
        char c = text[i];
        if (c == '(' || c == '[')
        {
            ++depth;
            ++i;
            continue;
        }
        if (c == ')' || c == ']')
        {
            if (depth > 0)
            {
                --depth;
            }
            ++i;
            continue;
        }
        if (depth == 0)
        {
            // "," is a plain separator; word keywords need word boundaries
            const bool bound_ok = (kw == ",") || (i == 0) ||
                std::isspace(static_cast<unsigned char>(text[i - 1]));
            const bool next_ok = (kw == ",") ||
                (i + kw.size() >= n ||
                 std::isspace(static_cast<unsigned char>(text[i + kw.size()])) ||
                 text[i + kw.size()] == '(' || text[i + kw.size()] == '[');
            if (bound_ok && next_ok &&
                i + kw.size() <= n &&
                text.compare(i, kw.size(), kw) == 0)
            {
                std::string part = trim_copy(text.substr(start, i - start));
                if (!part.empty() && part != kw)
                {
                    out.push_back(part);
                }
                start = i + kw.size();
                i = start;
                continue;
            }
        }
        ++i;
    }
    std::string part = trim_copy(text.substr(start));
    if (!part.empty())
    {
        out.push_back(part);
    }
    return out;
}

/// Split a comma-separated list inside `[...]`.
std::vector<std::string> split_commas(const std::string& inner)
{
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= inner.size(); ++i)
    {
        if (i == inner.size() || inner[i] == ',')
        {
            std::string v = trim_copy(inner.substr(start, i - start));
            if (!v.empty())
            {
                out.push_back(v);
            }
            start = i + 1;
        }
    }
    return out;
}

bool parse_atom(const std::string& text, AvailabilityCondition& cond)
{
    std::string t = trim_copy(text);
    if (t.empty())
    {
        return false;
    }
    // canonical comparison operators (longest/most specific first)
    static const char* OPS[] = {"==", "!=", ">=", "<=", ">", "<"};
    for (const char* op: OPS)
    {
        std::string opstr(op);
        const std::size_t pos = t.find(opstr);
        if (pos != std::string::npos)
        {
            std::string param = trim_copy(t.substr(0, pos));
            std::string rhs = trim_copy(t.substr(pos + opstr.size()));
            if (param.empty() || rhs.empty())
            {
                return false;
            }
            cond.param = param;
            cond.op = opstr;
            if (opstr == "==" && rhs.find('/') != std::string::npos)
            {
                // equality may carry a slash-separated value list
                std::size_t beg = 0;
                while (beg <= rhs.size())
                {
                    std::size_t sl = rhs.find('/', beg);
                    std::string v = trim_copy(rhs.substr(beg,
                        sl == std::string::npos ? std::string::npos : sl - beg));
                    if (!v.empty())
                    {
                        cond.values.push_back(v);
                    }
                    if (sl == std::string::npos)
                    {
                        break;
                    }
                    beg = sl + 1;
                }
                return !cond.values.empty();
            }
            cond.values.push_back(rhs);
            return true;
        }
    }
    // canonical containment: param contains value (e.g. a Vector holds an element)
    {
        const std::string marker = " contains ";
        const std::size_t cp = t.find(marker);
        if (cp != std::string::npos)
        {
            std::string param = trim_copy(t.substr(0, cp));
            std::string rhs = trim_copy(t.substr(cp + marker.size()));
            if (param.empty() || rhs.empty())
            {
                return false;
            }
            cond.param = param;
            cond.op = "contains";
            cond.values.push_back(rhs);
            return true;
        }
    }
    // canonical in-list: param in [v1, v2]
    std::size_t i = 0;
    while (i < t.size())
    {
        if (std::isspace(static_cast<unsigned char>(t[i])))
        {
            std::size_t j = i + 1;
            while (j < t.size() && std::isspace(static_cast<unsigned char>(t[j])))
            {
                ++j;
            }
            if (j + 2 <= t.size() && t.compare(j, 2, "in") == 0 &&
                (j + 2 == t.size() || std::isspace(static_cast<unsigned char>(t[j + 2]))))
            {
                std::string param = trim_copy(t.substr(0, i));
                std::string rest = trim_copy(t.substr(j + 2));
                if (param.empty() || rest.size() < 2 ||
                    rest[0] != '[' || rest[rest.size() - 1] != ']')
                {
                    return false;
                }
                cond.param = param;
                cond.op = "in";
                cond.values = split_commas(rest.substr(1, rest.size() - 2));
                return !cond.values.empty();
            }
            i = j;
        }
        else
        {
            ++i;
        }
    }
    return false;
}

AvailabilityExpr parse_or(const std::string& text);
AvailabilityExpr parse_and(const std::string& text);
AvailabilityExpr parse_unary(const std::string& text);

AvailabilityExpr parse_or(const std::string& text)
{
    std::vector<std::string> parts = split_top_keyword(text, "or");
    if (parts.size() == 1)
    {
        return parse_and(text);
    }
    AvailabilityExpr node;
    node.op = "or";
    for (std::size_t k = 0; k < parts.size(); ++k)
    {
        node.children.push_back(parse_and(parts[k]));
    }
    return node;
}

AvailabilityExpr parse_and(const std::string& text)
{
    std::vector<std::string> parts;
    {
        std::vector<std::string> a = split_top_keyword(text, "and");
        for (std::size_t k = 0; k < a.size(); ++k)
        {
            std::vector<std::string> c = split_top_keyword(a[k], ",");
            for (std::size_t m = 0; m < c.size(); ++m)
            {
                if (!c[m].empty() && c[m] != "," && c[m] != "and")
                {
                    parts.push_back(c[m]);
                }
            }
        }
    }
    if (parts.size() == 1)
    {
        return parse_unary(parts[0]);
    }
    AvailabilityExpr node;
    node.op = "and";
    for (std::size_t k = 0; k < parts.size(); ++k)
    {
        node.children.push_back(parse_unary(parts[k]));
    }
    return node;
}

AvailabilityExpr parse_unary(const std::string& text)
{
    std::string t = trim_copy(text);
    if (t.size() >= 2 && t[0] == '(' && t[t.size() - 1] == ')')
    {
        return parse_or(t.substr(1, t.size() - 2));
    }
    AvailabilityExpr leaf;
    if (parse_atom(t, leaf.condition))
    {
        return leaf;
    }
    return leaf;
}

} // namespace

std::string AvailabilityCondition::to_string() const
{
    if (values.empty())
    {
        return {};
    }
    if (op == "in")
    {
        std::string s = param + " in [";
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            if (i)
            {
                s += ", ";
            }
            s += values[i];
        }
        return s + "]";
    }
    if (op == "contains")
    {
        return param + " contains " + values[0];
    }
    // single-value comparison (==, !=, >, >=, <, <=)
    return param + op + values[0];
}

std::string AvailabilityExpr::to_string() const
{
    if (is_leaf())
    {
        return condition.to_string();
    }
    std::string sep = (op == "or") ? " or " : " and ";
    std::string s;
    for (std::size_t i = 0; i < children.size(); ++i)
    {
        if (i)
        {
            s += sep;
        }
        if (children[i].is_leaf())
        {
            s += children[i].to_string();
        }
        else
        {
            s += "(" + children[i].to_string() + ")";
        }
    }
    return s;
}

AvailabilityExpr parse_availability(const std::string& raw)
{
    std::string t = trim_copy(raw);
    if (t.empty())
    {
        return AvailabilityExpr();
    }
    return parse_or(t);
}

} // namespace ModuleIO

#include <cstring>
#include <iostream>

class SimpleString
{
  private:
    char m_data[80];
    size_t m_length;

  public:
    // Default constructor
    SimpleString() : m_length(0)
    {
    }

    // constructor
    SimpleString(const char* str)
    {
        m_length = std::strlen(str);
        std::strcpy(m_data, str);
    }

    // Copy constructor
    SimpleString(const SimpleString& other)
    {
        m_length = other.m_length;
        std::strcpy(m_data, other.m_data);
    }

    // Get string length
    size_t length() const
    {
        return m_length;
    }

    // Get string content
    const char* c_str() const
    {
        return m_data;
    }

    // Overload the assignment operator
    SimpleString& operator=(const SimpleString& other)
    {
        if (this != &other)
        {
            m_length = other.m_length;
            std::strcpy(m_data, other.m_data);
        }
        return *this;
    }
};

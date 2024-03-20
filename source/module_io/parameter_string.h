#include <cstring>
#include <iostream>

/**
 * @brief   This class is a simple string class that contains an array of characters and a variable that represents the
 * length of the string. It has a default constructor, a copy constructor, and two public functions that return the
 * length and content of the string, respectively. The assignment operator is also overridden to be able to copy one
 * SimpleString object to another.
 *
 * @param m_data store the value of string
 * @param m_length the length of string
 */
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

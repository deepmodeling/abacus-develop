#include <iostream>
#include <cstring>

class SimpleString {
private:
    char m_data[80];
    size_t m_length;

public:
    // 默认构造函数
    SimpleString() : m_length(0) {}

    // 构造函数
    SimpleString(const char* str) {
        m_length = std::strlen(str);
        std::strcpy(m_data, str);
    }

    // 拷贝构造函数
    SimpleString(const SimpleString& other) {
        m_length = other.m_length;
        std::strcpy(m_data, other.m_data);
    }


    // 获取字符串长度
    size_t length() const {
        return m_length;
    }

    // 获取字符串内容
    const char* c_str() const {
        return m_data;
    }

    // 重载赋值运算符
    SimpleString& operator=(const SimpleString& other) {
        if (this != &other) {
            m_length = other.m_length;
            std::strcpy(m_data, other.m_data);
        }
        return *this;
    }
};


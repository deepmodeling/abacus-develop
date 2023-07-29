#include <iostream>
#include <vector>

/**
 * @brief   This code defines a template class, SimpleVector, which is a simple dynamic array that can store any type of
 *          data. It contains private members data_ and size_, which represent the data and size of the dynamic array,
 *          respectively. It provides a default constructor and a constructor with an initializer list. It also provides
 *          the push_back() function to add elements to dynamic arrays, and overrides the [] operator to fetch elements
 *          by index. The size() function returns the size of the dynamic array.
 *
 * @param data_ store the value of vetor
 * @param size_ the size of vetor
 */
template <typename T>
class SimpleVector
{
  private:
    T data_[100];
    size_t size_;

  public:
    SimpleVector() : size_(0)
    {
    }

    SimpleVector(std::initializer_list<T> values) : size_(0)
    {
        for (auto it = values.begin(); it != values.end(); ++it)
        {
            push_back(*it);
        }
    }
    void push_back(const T& value)
    {
        data_[size_] = value;
        ++size_;
    }

    T& operator[](size_t index)
    {
        return data_[index];
    }

    const T& operator[](size_t index) const
    {
        return data_[index];
    }

    size_t size() const
    {
        return size_;
    }
};

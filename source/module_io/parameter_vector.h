#include <iostream>
#include <vector>

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

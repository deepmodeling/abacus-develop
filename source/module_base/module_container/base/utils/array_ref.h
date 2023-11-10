#ifndef BASE_UTILS_ARRAY_REF_H_
#define BASE_UTILS_ARRAY_REF_H_

#include <cstddef>
#include <cstdint>
#include <base/macros/macros.h>

namespace base {
namespace utils {

template <typename T>
class array_ref final {

  private:
    T* data_;
    size_t length_;

  public:
    /* implicit */ constexpr array_ref() : data_(nullptr), length_(0) {}
    /* implicit */ constexpr array_ref(T* data, size_t length) : data_(data), length_(length) {}
    /* implicit */ constexpr array_ref(T* begin, T* end) : data_(begin), length_(end - begin) {}
    explicit constexpr array_ref(const T& item) : data_(&item), length_(1) {}

    // Construct from a std::vector.
    template <typename A>
    /* implicit */ array_ref(const std::vector<T, A>& vec)
      : data_(vec.data()), length_(vec.size()) {
        static_assert(
            !std::is_same<T, bool>::value,
            "array_ref<bool> cannot be constructed from a std::vector<bool> bitfield.");
    }

    // Construct from a std::array.
    template <size_t N>
    /* implicit */ constexpr array_ref(const std::array<T, N>& arr) : data_(arr.data()), length_(N) {}

    // Construct from a std::initializer_list.
    /* implicit */ constexpr array_ref(const std::initializer_list<T>& list)
      : data_(list.begin()), length_(list.size()) {}

    constexpr const T* begin() const { return data_; }
    constexpr const T* end()   const { return data_ + length_; }

    constexpr bool empty() const {
        return length_ == 0;
    }

    constexpr const T* data() const {
        return data_;
    }

    constexpr size_t size() const {
        return length_;
    }

    constexpr const T& front() const {
        REQUIRES_OK(!empty(), "Cannot call front() on an empty array_ref.");
        return data_[0];
    }

    constexpr const T& back() const {
        REQUIRES_OK(!empty(), "Cannot call back() on an empty array_ref.");
        return data_[length_ - 1];
    }

    constexpr bool equals(const array_ref<T>& rhs) const {
        return length_ == rhs.size() && std::equal(begin(), end(), rhs.begin());
    }

    constexpr const T& operator[](size_t index) const {
        return data_[index];
    }

    template <typename U>
    typename std::enable_if<std::is_same<U, T>::value, array_ref<T>>::type&
    operator=(U&& Temporary) = delete;

    template <typename U>
    typename std::enable_if<std::is_same<U, T>::value, array_ref<T>>::type&
    operator=(std::initializer_list<U>) = delete;

    std::vector<T> vec() const {
        return std::vector<T>(data_, data_ + length_);
    }
};

template <typename T>
std::ostream& operator<<(std::ostream& out, array_ref<T> arr) {
    int ii = 0;
    out << "[";
    for (const auto& item : arr) {
        if (ii++ > 0)
            out << ", ";
        out << item;
    }
    out << "]";
    return out;
}

template <typename T>
bool operator==(array_ref<T> a1, array_ref<T> a2) {
    return a1.equals(a2);
}

template <typename T>
bool operator!=(array_ref<T> a1, array_ref<T> a2) {
    return !a1.equals(a2);
}


} // namespace utils
} // namespace base

using int_array_ref = ::base::utils::array_ref<int64_t>;

#endif // BASE_UTILS_ARRAY_REF_H_
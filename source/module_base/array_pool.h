#ifndef ARRAY_POOL_H
#define ARRAY_POOL_H


namespace ModuleBase
{
    /**
     * @brief Array_Pool is a class designed for dynamically allocating a two-dimensional array
     *  with all its elements contiguously arranged in memory. Compared to a two-dimensional vector,
     *  it offers better data locality because all elements are stored in a continuous block of memory.
     *  
     * @tparam T 
     */
    template <typename T>
    class Array_Pool
    {
    public:
        Array_Pool();
        Array_Pool(const int nr, const int nc);
        Array_Pool(Array_Pool<T>&& other);
        Array_Pool& operator=(Array_Pool<T>&& other);
        ~Array_Pool();
        Array_Pool(const Array_Pool<T>& other) = delete;
        Array_Pool& operator=(const Array_Pool& other) = delete;

        T** get_ptr_2D() const { return ptr_2D; }
        T* get_ptr_1D() const { return ptr_1D; }
        int get_nr() const { return nr; }
        int get_nc() const { return nc; }
        T* operator[](const int ir) const { return ptr_2D[ir]; }
    private:
        T** ptr_2D;
        T* ptr_1D;
        int nr;
        int nc;
    };

    template <typename T>
    Array_Pool<T>::Array_Pool() : ptr_2D(nullptr), ptr_1D(nullptr), nr(0), nc(0)
    {
    }

    template <typename T>
    Array_Pool<T>::Array_Pool(const int nr, const int nc) // Attention: uninitialized
    {
        this->nr = nr;
        this->nc = nc;
        ptr_1D = new T[nr * nc];
        ptr_2D = new T*[nr];
        for (int ir = 0; ir < nr; ++ir)
            ptr_2D[ir] = &ptr_1D[ir * nc];
    }

    template <typename T>
    Array_Pool<T>::~Array_Pool()
    {
        delete[] ptr_2D;
        delete[] ptr_1D;
    }

    template <typename T>
    Array_Pool<T>::Array_Pool(Array_Pool<T>&& other)
    {
        ptr_2D = other.ptr_2D;
        ptr_1D = other.ptr_1D;
        nr = other.nr;
        nc = other.nc;
        other.ptr_2D = nullptr;
        other.ptr_1D = nullptr;
        other.nr = 0;
        other.nc = 0;
    }

    template <typename T>
    Array_Pool<T>& Array_Pool<T>::operator=(Array_Pool<T>&& other)
    {
        if (this != &other)
        {
            delete[] ptr_2D;
            delete[] ptr_1D;
            ptr_2D = other.ptr_2D;
            ptr_1D = other.ptr_1D;
            nr = other.nr;
            nc = other.nc;
            other.ptr_2D = nullptr;
            other.ptr_1D = nullptr;
            other.nr = 0;
            other.nc = 0;
        }
        return *this;
    }

}
#endif
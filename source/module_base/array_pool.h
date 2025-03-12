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
        Array_Pool() = default;
        Array_Pool(const int nr_in, const int nc_in);
        Array_Pool(Array_Pool<T>&& other);
        Array_Pool& operator=(Array_Pool<T>&& other);
        ~Array_Pool();
        Array_Pool(const Array_Pool<T>& other) = delete;
        Array_Pool& operator=(const Array_Pool& other) = delete;

        T** get_ptr_2D() const { return this->ptr_2D; }
        T* get_ptr_1D() const { return this->ptr_1D; }
        int get_nr() const { return this->nr; }
        int get_nc() const { return this->nc; }
        T* operator[](const int ir) const { return this->ptr_2D[ir]; }
    private:
        T** ptr_2D = nullptr;
        T* ptr_1D = nullptr;
        int nr = 0;
        int nc = 0;
    };

    template <typename T>
    Array_Pool<T>::Array_Pool(const int nr_in, const int nc_in) // Attention: uninitialized
        : nr(nr_in),
          nc(nc_in)
    {
        this->ptr_1D = new T[nr * nc]();
        this->ptr_2D = new T*[nr];
        for (int ir = 0; ir < nr; ++ir)
            this->ptr_2D[ir] = &this->ptr_1D[ir * nc];
    }

    template <typename T>
    Array_Pool<T>::~Array_Pool()
    {
        delete[] this->ptr_2D;
        delete[] this->ptr_1D;
    }

    template <typename T>
    Array_Pool<T>::Array_Pool(Array_Pool<T>&& other)
        : ptr_2D(other.ptr_2D),
          ptr_1D(other.ptr_1D),
          nr(other.nr),
          nc(other.nc)
    {
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
            delete[] this->ptr_2D;
            delete[] this->ptr_1D;
            this->ptr_2D = other.ptr_2D;
            this->ptr_1D = other.ptr_1D;
            this->nr = other.nr;
            this->nc = other.nc;
            other.ptr_2D = nullptr;
            other.ptr_1D = nullptr;
            other.nr = 0;
            other.nc = 0;
        }
        return *this;
    }

}
#endif
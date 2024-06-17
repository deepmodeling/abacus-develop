#ifndef PARALLEL_K2D_H
#define PARALLEL_K2D_H


/***
 * This is a class to realize k-points parallelism in LCAO code.
 * It is now designed only to work with 2D Eigenvalue solver parallelism.
 */

class Parallel_K2D
{
public:
    /**
     * Define Parallel_K2D as a signleton class
     */
    static Parallel_K2D& get_instance()
    {
        static Parallel_K2D instance;
        return instance;
    }
    /// delete copy constructor
    Parallel_K2D(const Parallel_K2D&) = delete;
    /// delete copy assignment
    Parallel_K2D& operator=(const Parallel_K2D&) = delete;

    /**
     * Public member functions
     */
    /// set the number of k-points
    void set_kpar(int kpar) { kpar_ = kpar; }
    /// get the number of k-points
    int get_kpar() const { return kpar_; }

private:
    /**
     * Public member functions
     */
    /// private constructor
    Parallel_K2D();
    /// private destructor
    ~Parallel_K2D();

private:
    /**
     * Private member variables
     */
    int kpar_ = 1;

};

#endif
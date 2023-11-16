template <typename T>
class Associate_Laguerre
{
    public:
        Associate_Laguerre();
        ~Associate_Laguerre();
        
        T value(const int& n, const int& alpha, T x);
        void generate(const int& n, const int& alpha, T* x, T* Laguerre, const int& size);

    private:
        int alpha = 0;
        int n = 0; // order of the polynomial
};
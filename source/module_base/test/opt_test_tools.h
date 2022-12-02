class LinearEqu
{
public:
    LinearEqu();
    ~LinearEqu();
    void get_Ap(double **A, double *p, double *Ap, int nx = 3, int ny = 3);
    double func(double *x);
    void dfuncdx(double *x, double *gradient);
    double dfuncdstp(double *x, double *p);
    double *b;
    double **A;
    int nx = 3;
};

class MinFunc
{
public:
    MinFunc(){};
    ~MinFunc(){};
    double func(double *x);
    void dfuncdx(double *x, double *gradient);
    double dfuncdstp(double *x, double *p);
    double *x;
};

class TestTools
{
public:
    TestTools()
    {
        this->nx = le.nx;
    }
    double func(double *x, int func_label)
    {
        double result = 0.;
        if (func_label==0) result = le.func(x);
        else if (func_label==1) result = mf.func(x);
        return result;
    }
    void dfuncdx(double *x, double *gradient, int func_label)
    {
        if (func_label==0) le.dfuncdx(x, gradient);
        else if (func_label==1) mf.dfuncdx(x, gradient);
    }
    double dfuncdstp(double *x, double *p, int func_label)
    {
        double result = 0.;
        if (func_label==0) result = le.dfuncdstp(x, p);
        else if (func_label==1) result = mf.dfuncdstp(x, p);
        return result;
    }

    int nx = 0;
    LinearEqu le;
    MinFunc mf;
};

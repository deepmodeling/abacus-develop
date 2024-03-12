#ifndef MATRIX_H
#define MATRIX_H
class matrix
{
	/* data */
public:
	double *c;    /* Holds the data */
	int nr;
	int nc;   /* Number of rows and columns */
	int hermetian;
	static int mCount;

	/* Constructors and destructor */
	matrix(const int nrows = 1,const int ncols = 1);
	matrix(const matrix &m1); /* copy constructor */
	~matrix();

	void create(const int nrow,const int ncol);
	void operator=(const matrix &m1); /* Nonstandard: returns void */

	double &operator()(const int ir,const int ic)
	{ return c[nc*ir + ic]; };

	const double &operator()(const int ir,const int ic)const
	{ return c[nc*ir + ic]; };

//	inline double &operator()(int i,int j) const
//	    { return c[nc*i+j]; }

	friend matrix operator+(const matrix &m1, const matrix &m2);
	friend matrix operator-(const matrix &m1, const matrix &m2);
	friend matrix operator*(const matrix &m1, const matrix &m2);
	friend matrix operator*(const double &s, const matrix &m);
	friend matrix operator*(const matrix &m, const double &s);
	void operator*=(const double &s);
	void operator+=(const matrix &m);
	void operator-=(const matrix &m);

	/* member function: */
	matrix Inverse(void);

	double det(void);

	/* Does the memory allocations of the constructor */
	void init(int nrows, int ncols);

	/* Free up memory */
	void freemem(void);

	/* zero out all the entries */
	void zero_out(void);
};

matrix transpose(matrix m);

#endif // MATRIX_H

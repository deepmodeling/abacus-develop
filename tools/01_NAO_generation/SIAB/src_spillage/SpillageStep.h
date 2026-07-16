#ifndef SPILLAGE_STEP
#define SPILLAGE_STEP

#include "common.h"
#include "Step_Data.h"
#include "Type_Information.h"

// use orbital index to find the information
// of each orbital.
struct Way2Data // dimension : nwfc2
{
	int type;
	int i; // atom index
	int L;
	int N;
	int m; // mohan add 2010-06-16

	int iw0; // to label <J|Bloch>, which not include 'N', but not used now,
	int iw00; // to <J|J>
	int ic;

	bool average; // mohan add 2009-07-12
};

struct Way2c4 // dimension : nchi
{
	int type;
	int L;
	int N;
	char spd; //mohan add 2010-04-17, the orbital type.
};

// (1) It's a memeber of MultiZeta.
// (2) define the operations for each step.
// (3) Especially it defines many complicated index.
// (4) It also defines the main routine to calculate the spillage value.
// (5) Also the update of Q and S matrix using complicated index.
// (6) It also contains four calss: Type_Information; Step_Data; Way2c4; Way2Data;
class SpillageStep
{
	public:
	SpillageStep();
	~SpillageStep();

	Type_Information *info; // dimension : type
	Step_Data *data;	
	Way2c4 *wayc;
	Way2Data *wayd;

	void set_info(const int &ntype_in);
	void set_nwfc(void);
	void set_iw00(void);
	void set_nchi(void);

	// allocate data for current level.
	void allocate_data( const int &istep);

	void init_QS_matrix( const int &is);
	double get_spillage( const int &is, const int &ic, const int &ie);
	void updateQS( const int &is );

	int ntype;
	int nchi; // total radial wave functions(T,L,N);
	int nwfc2;
	int ne;
	int ilevel;

	private:

	void newQ( const int &is, const int &ic, const int &ie, const int &ik,
			const double &c4_now, const double &c4_old);
	void newS( const int &is, const int &ic, const int &ie, const int &ik,
			const double &c4_now ,const double &c4_old);

	int test;
	int nwfc;
};
#endif

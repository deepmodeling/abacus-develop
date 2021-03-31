#ifndef MD_FUNC_H
#define MD_FUNC_H

class MD_func
{
    public:
    MD_func(){};
    ~MD_func(){};
    void initMD();
	bool RestartMD();
    void mstout(int step);
//	void RemoveMovementOfCenterOfMass();
	double GetAtomKE();
//	void MakeIntCoeff();
	void InitVelocity();
//	void MonitorMeanSquareDisplacement(int step);
//	void OutputIonVelocityAndStats(int iter);
	void OutputMDHeader();
	void OutputMDGeom(int iter);
	void ReadNewTemp(int step);
	string intTurnTostring(int iter,string path);
	void connection0();
	void connection1();
	void connection2();
	void callforce();
    void moveatoms(int step);
    void printpos(string file,int iter);
//        void md_release();
    void scalevel();
	double MAXVALF();
	double Conserved(double KE, double PE);
};
#endif
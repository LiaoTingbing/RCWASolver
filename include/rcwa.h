#pragma once
#include "common.h"
#include "userfunc.h"
#include "convulationmatrix.h"

class RCWA
{
	vector<cx_mat> Index;
	vec LayerPos, x, y, z;
	int ku, kv;
	double  lambda, theta, phi, n_lower, n_upper;

	double Rs, Rp, Ts, Tp;

public:

	RCWA();
 
	RCWA(DataRCWA & In);

	~RCWA();
 
	void Run();

	double getRs();
	double getRp();
	double getTs();
	double getTp();
 
	void set_lambda(double In) ;
	void set_theta(double In);
	void set_phi(double In);
	void set_n_lower(double In);
	void set_n_upper(double In);
	void set_Index(vector<cx_mat> In);


};


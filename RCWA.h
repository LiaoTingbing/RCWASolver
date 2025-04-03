#pragma once
#include "common.h"
#include "userFunc.h"
#include "convulationMatrix.h"

class RCWA
{
	cx_mat* Index;
	vec LayerPos, x, y, z;
	int ku, kv;
	double  lambda, theta, phi ;
	//cx_double ER_inc, UR_inc, ER_ref, UR_ref, ER_trn, UR_trn;
	double n_lower, n_upper;

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

	void set

};


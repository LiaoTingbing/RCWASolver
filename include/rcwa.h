#pragma once
#include "common.h"
#include "userfunc.h"
#include "convulationmatrix.h"

class RCWA
{
	vector<cx_mat> Index_;
	vec layer_pos_, x_pos_, y_pos_, z_pos_;
	int ku_, kv_;
	double  lambda_, theta_, phi_angle_, n_lower_value, n_upper_value;

	double Rs, Rp, Ts, Tp;

public:

	RCWA();
 
	RCWA(DataRCWA & In);

	~RCWA();
 
	void Run();

	//	计算结果
	double getRs();
	double getRp();
	double getTs();
	double getTp();
 
	//	
	void set_lambda(double In) ;
	void set_theta(double In);
	void set_phi(double In);
	void set_n_lower(double In);
	void set_n_upper(double In);
	void set_Index(vector<cx_mat> In);


};


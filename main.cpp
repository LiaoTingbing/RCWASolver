// RCWASolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "RCWA.h"
 
DataFile loadDATA();

int main()
{
  
	DataFile In = loadDATA();
	RCWA rcwa(In);
	rcwa.Run();
 
	return 0;
}

DataFile loadDATA()
{
	DataFile DATA;
	string filepath = "input/";
	DATA.LayerPos.load(filepath + "LayerPos.txt");
	loadTXT(DATA.ku, filepath + "ku.txt");
	loadTXT(DATA.kv, filepath + "kv.txt");
	loadTXT(DATA.lambda, filepath + "lambda.txt");
	loadTXT(DATA.theta, filepath + "theta.txt");
	loadTXT(DATA.phi, filepath + "phi.txt");
	loadTXT(DATA.n_lower, filepath + "n_lower.txt");
	loadTXT(DATA.n_upper, filepath + "n_upper.txt");

	DATA.x.load(filepath + "x.txt");
	DATA.y.load(filepath + "y.txt");
	DATA.z.load(filepath + "z.txt");
	DATA.layersNum = DATA.LayerPos.size();

	DATA.ER_inc = DATA.n_lower * DATA.n_lower;
	DATA.UR_inc = 1;
	DATA.ER_trn = DATA.n_upper * DATA.n_upper;
	DATA.UR_trn = 1;
	DATA.ER_ref = DATA.ER_inc;
	DATA.UR_ref = DATA.UR_inc;

	string sreal, simag;
	DATA.Index = new cx_mat[DATA.layersNum];
	mat IndexReal, IndexImag;
	for (int i = 0; i < DATA.layersNum;i++)
	{
		sreal = filepath + "Index_real_" + "z" + to_string(i + 1) + ".txt";
		simag = filepath + "Index_imag_" + "z" + to_string(i + 1) + ".txt";

		IndexReal.load(sreal);
		IndexImag.load(simag);
		DATA.Index[i] = IndexReal + iI * IndexImag;
	}

	



	return DATA;
}

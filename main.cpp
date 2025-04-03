// RCWASolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "RCWA.h"

DataFile loadDATA();
//void multiLambda();


int main()
{

	//multiLambda();
	DataFile dataIn = loadDATA();
	DataRCWA rcwaIn;

	rcwaIn.Index = dataIn.Index;
	rcwaIn.LayerPos = dataIn.LayerPos.col(0);
	rcwaIn.x = dataIn.x.col(0);
	rcwaIn.y = dataIn.y.col(0);
	rcwaIn.z = dataIn.z.col(0);
	rcwaIn.ku = round (dataIn.ku(0) );
	rcwaIn.kv =round( dataIn.kv(0) );
	rcwaIn.theta =  (dataIn.theta(0));
	rcwaIn.phi =  (dataIn.phi(0));
	rcwaIn.lambda=  (dataIn.lambda(0));
	rcwaIn.n_upper =  (dataIn.n_upper(0));
	rcwaIn.n_lower =   (dataIn.n_lower(0));
 
	RCWA rcwa(rcwaIn);
	rcwa.Run();

	return 0;
}

DataFile loadDATA()
{
	DataFile DATA;

	string filepath = "input/";
	DATA.LayerPos.load(filepath + "LayerPos.txt");
	DATA.ku.load(filepath + "ku.txt");
	DATA.kv.load(filepath + "kv.txt");
	DATA.theta.load(filepath + "theta.txt");
	DATA.phi.load(filepath + "phi.txt");
	DATA.lambda.load(filepath + "lambda.txt");
	DATA.n_lower.load(filepath + "n_lower.txt");
	DATA.n_upper.load(filepath + "n_upper.txt");
	DATA.x.load(filepath + "x.txt");
	DATA.y.load(filepath + "y.txt");
	DATA.z.load(filepath + "z.txt");
	//DATA.layersNum = DATA.LayerPos.size();

	string sreal, simag;
	DATA.Index = new cx_mat[DATA.LayerPos.size()];
	mat IndexReal, IndexImag;
	for (int i = 0; i < DATA.LayerPos.size();i++)
	{
		sreal = filepath + "Index_real_" + "z" + to_string(i + 1) + ".txt";
		simag = filepath + "Index_imag_" + "z" + to_string(i + 1) + ".txt";

		IndexReal.load(sreal);
		IndexImag.load(simag);
		DATA.Index[i] = IndexReal + iI * IndexImag;
	}
	return DATA;
}

 
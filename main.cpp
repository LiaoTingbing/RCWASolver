// RCWASolver.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "RCWA.h"

DataFile loadDATA();


int main()
{

	DataFile dataIn = loadDATA();

	DataRCWA rcwaIn;

	// 器件数据
	rcwaIn.Index = dataIn.IndexS[0];           // 第一个波长数据
	rcwaIn.LayerPos = dataIn.LayerPos.col(0);
	rcwaIn.x = dataIn.x.col(0);
	rcwaIn.y = dataIn.y.col(0);
	rcwaIn.z = dataIn.z.col(0);
	rcwaIn.n_upper = (dataIn.n_upper(0));
	rcwaIn.n_lower = (dataIn.n_lower(0));
	// 谐波数
	rcwaIn.ku = round(dataIn.ku(0));
	rcwaIn.kv = round(dataIn.kv(0));
	// 波源数据
	rcwaIn.theta = (dataIn.theta(0));
	rcwaIn.phi = (dataIn.phi(0));
	rcwaIn.lambda = (dataIn.lambda(0));


	RCWA rcwa(rcwaIn);

	vec Rs(dataIn.lambda.size());
	vec Rp(dataIn.lambda.size());
	vec Ts(dataIn.lambda.size());
	vec Tp(dataIn.lambda.size());

	for (size_t i = 0; i < dataIn.lambda.size(); i++)
	{
		rcwa.set_lambda(dataIn.lambda(i));
		rcwa.set_n_lower(dataIn.n_lower(i));
		rcwa.set_n_upper(dataIn.n_upper(i));
		rcwa.set_Index(dataIn.IndexS[i]);
		rcwa.Run();
		Rs(i) = rcwa.getRs();
		Rp(i) = rcwa.getRp();
		Ts(i) = rcwa.getTs();
		Tp(i) = rcwa.getTp();

	}

	//Rs.print();
	//Rp.print();
	//Ts.print();
	//Tp.print();

	Rs.save("output/Rs.txt", arma::raw_ascii);
	Rp.save("output/Rp.txt", arma::raw_ascii);
	Ts.save("output/Ts.txt", arma::raw_ascii);
	Tp.save("output/Tp.txt", arma::raw_ascii);
	dataIn.lambda.save("output/lambda.txt", arma::raw_ascii);

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
	mat IndexReal, IndexImag;
	cx_mat tmp;

	DATA.IndexS.resize(DATA.lambda.size());

	for (size_t j = 0; j < DATA.lambda.size(); j++)  //
	{
		for (size_t i = 0; i < DATA.LayerPos.size();i++)  // f
		{
			DATA.IndexS[j].resize(DATA.LayerPos.size());
			sreal = filepath + "Index_real_" + "z_" + to_string(i + 1) + "_" + to_string(j+1) + ".txt";
			simag = filepath + "Index_imag_" + "z_" + to_string(i + 1) + "_" + to_string(j+1) + ".txt";

			IndexReal.load(sreal);
			IndexImag.load(simag);
			tmp = IndexReal + iI * IndexImag;
			DATA.IndexS[j][i] = tmp;
		}
	}

	return DATA;
}


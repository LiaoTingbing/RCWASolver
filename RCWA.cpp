#pragma once
#include "RCWA.h"

RCWA::RCWA()
{
}

RCWA::~RCWA()
{
}

RCWA::RCWA(Dev& dev, Source& sc ,rcwaPara Para)
{
	this->dev = dev;
	this->sc = sc;
	this->Para = Para;
}

void RCWA::initiliaze()
{
	cout << "初始化\n";
	//clock_t t1, t2, t3, t4;
	//t1 = clock();

	DATA.Nx = dev.x.size();
	DATA.Ny = dev.y.size();
	DATA.Ncx = dev.x.size() - 1;
	DATA.Ncy = dev.y.size() - 1;
	DATA.Nc = DATA.Ncx * DATA.Ncy;

	DATA.Lx = dev.x(DATA.Nx - 1) - dev.x(0);
	DATA.Ly = dev.y(DATA.Ny - 1) - dev.y(0);
	DATA.thetaArc = sc.theta * pi / 180.0;
	DATA.phiArc = sc.phi * pi / 180.0;
	DATA.incident_direction = {
		sin(DATA.thetaArc) * cos(DATA.phiArc),
		sin(DATA.thetaArc) * sin(DATA.phiArc),
		cos(DATA.thetaArc)
	};
	DATA.k0 = 2 * pi / sc.lambda;

	cout << "\t计算入射波矢:";
	cout << "\ttheta=" << to_string(sc.theta);
	cout << "\tphi=" << to_string(sc.phi);
	cout << endl;

	DATA.ER_inc = 1;
	DATA.UR_inc = 1;
	DATA.k_inc = sqrt(DATA.ER_inc) * DATA.incident_direction;
	DATA.kx_inc = DATA.k_inc(0);
	DATA.ky_inc = DATA.k_inc(1);
	DATA.kz_inc = DATA.k_inc(2);

	cout << "\t计算周期矢量";
	cout << "\tLx = " << to_string(DATA.Lx);
	cout << "\tLy = " << to_string(DATA.Ly);
	cout << endl;

	DATA.Tx = 2.0 * pi / DATA.Lx / DATA.k0;
	DATA.Ty = 2.0 * pi / DATA.Ly / DATA.k0;

	cout << "\t计算谐波展开:";
	cout << "\tku=" << to_string(Para.ku);
	cout << "\tkv=" << to_string(Para.kv);
	cout << endl;



	DATA.m = linspace(-Para.ku, Para.ku, 2 * Para.ku + 1);
	DATA.n = linspace(-Para.kv, Para.kv, 2 * Para.kv + 1);

	DATA.M = DATA.m.size();
	DATA.N = DATA.n.size();
	DATA.Nh = DATA.M * DATA.N;


	cx_mat kx_mn(DATA.M, DATA.N);
	cx_mat ky_mn(DATA.M, DATA.N);
	
	DATA.kx_mn = cx_mat(DATA.M, DATA.N);
	DATA.ky_mn = cx_mat(DATA.M, DATA.N);


	for (size_t i = 0; i < DATA.M; i++)
	{
		for (size_t j = 0; j < DATA.N; j++)
		{
			DATA.kx_mn(i, j) = DATA.kx_inc - DATA.m(i) * DATA.Tx;
			DATA.ky_mn(i, j) = DATA.ky_inc - DATA.n(j) * DATA.Ty;
		}
	}

	//real(kx_mn).print();
	//cout << endl;
	//real(ky_mn).print();
	//kx_mn.reshape(Nh, 1).print();
	DATA.Kx = diagmat(DATA.kx_mn.reshape(DATA.Nh, 1));
	  DATA.Ky = diagmat(DATA.ky_mn.reshape(DATA.Nh, 1));


	//diagvec(Kx).print();

	//cout << "\t初始化散射矩阵\n";
	  I = cx_mat(DATA.Nh, DATA.Nh, fill::eye);
	  II = cx_mat(2 * DATA.Nh, 2 * DATA.Nh, fill::eye);
	  Z = cx_mat(DATA.Nh, DATA.Nh);
	  ZZ = cx_mat(DATA.Nh * 2, DATA.Nh * 2);

	G.S11 = DATA.ZZ;
	G.S12 = DATA.II;
	G.S21 = DATA.II;
	G.S22 = DATA.ZZ;
}

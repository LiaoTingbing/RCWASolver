#pragma once
#include "common.h"
#include "userFunc.h"
#include "convulationMatrix.h"

class Device
{
	cx_mat* Index;
	vec LayerPos, x, y, z;
	int ku, kv;
	double  lambda, theta, phi, thetaArc, phiArc;

	size_t layersNum;
	size_t Nx;
	size_t Ny;

	size_t Ncx;
	size_t Ncy;

	size_t Nc;

	double Lx;
	double Ly;

	Smatrix G;

	double Rs, Ts, Rp, Tp;

public:

	Device();
	~Device();

	void initiliaze();
	void refSlove();
	void layerSole();
	void trnSolve();
	void RTsolve();

	void Run() {
		initiliaze();
		refSlove();
		layerSole();
		trnSolve();
		RTsolve();
	}

	void RCWA();

	double get_Rs();
	double get_Rp();
	double get_Ts();
	double get_Tp();

};


#pragma once
#include "common.h"
class Device
{
	cx_mat* IndexZ;
	vec Izpos, x, y, z;
	int ku, kv;
	double  lambda, theta, phi,thetaArc,phiArc;
	size_t layersNum;
	
	//
	size_t Nx  ;
	size_t Ny  ;


	size_t Ncx ;
	size_t Ncy;

	size_t Nc  ;

	double Lx  ;
	double Ly ;

	//
	Smatrix G;


public:
 

	Device();


	~Device();


	void RCWA();
 



};


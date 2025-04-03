#pragma once
#include "common.h"

class RCWA
{
	Smatrix ST , SR , G;
	Source sc;
	Dev dev;
	rcwaDATA DATA;

	rcwaPara Para;
	
	cx_mat I ;
	cx_mat II ;
	cx_mat Z ;
	cx_mat ZZ;

public:
	RCWA();
	RCWA(Dev& dev, Source& sc, rcwaPara Para);
	~RCWA();

	void initiliaze();
	void gapSolve();
	void refSlove();
	void layerSole();
	void trnSolve();
	void RTsolve();

	void Run() {
		initiliaze();
		gapSolve();
		refSlove();
		layerSole();
		trnSolve();
		RTsolve();
	}
};



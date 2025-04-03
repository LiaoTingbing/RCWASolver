#pragma once
#include "RCWA.h"

void loadTXT(double& s, string filename);

void loadTXT(int & s, string filename);

cx_mat MatrixConnect(cx_mat  A, cx_mat B, cx_mat  C, cx_mat  D);

Smatrix SconnectRight(Smatrix& G, Smatrix& S);


 




#include "Device.h"


 


void loadTXT(double& s, string filename);

void loadTXT(int & s, string filename);

cx_mat matrix_to_vector_row(cx_mat & K);

WVLmatrix Gap(cx_mat & Kx, cx_mat& Ky, cx_mat& I, cx_mat& Z);

cx_mat Matrix_merging(cx_mat  A, cx_mat  B, cx_mat  C, cx_mat  D);

WVLmatrix solution_in_Homogeneous_Layers(cx_mat  & Kx, cx_mat & Ky, cx_double & ER, cx_double & UR);

Smatrix SconnectRight(Smatrix& G, Smatrix& S);
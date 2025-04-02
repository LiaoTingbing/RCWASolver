#pragma once
#include "Device.h"


 


void loadTXT(double& s, string filename);

void loadTXT(int & s, string filename);

WVLmatrix Gap(cx_mat & Kx, cx_mat& Ky, cx_mat& I, cx_mat& Z);

cx_mat MatrixConnect(cx_mat  A, cx_mat B, cx_mat  C, cx_mat  D);

WVLmatrix solution_in_Homogeneous_Layers(
	cx_mat  & Kx, cx_mat & Ky, cx_double   ER, cx_double   UR);

Smatrix SconnectRight(Smatrix& G, Smatrix& S);


Mesh ndgrid(vec x , vec y);

cx_mat Liu_conv_acc(
	cx_mat & Kx, cx_mat &Ky, cx_mat &Z, size_t Nh, 
	size_t Ncx, size_t Ncy, cx_mat &ER,double k0, mat X, mat Y);

WVLmatrix solution_in_inHomogeneous_Layers(
	cx_mat &Kx, cx_mat& Ky, cx_mat& ERC, cx_mat& URC);



cx_mat  invD( const cx_mat & In);  //∂‘Ω«æÿ’Û«ÛƒÊ
mat  invD(const mat& In);  //∂‘Ω«æÿ’Û«ÛƒÊ




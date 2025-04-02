#pragma once
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
 
const double pi = 3.1415926;

const cx_double iI(0, 1);

struct  Smatrix
{
	cx_mat S11;
	cx_mat S12;
	cx_mat S21;
	cx_mat S22;
};

struct  WVLmatrix
{
	cx_mat W; 
	cx_mat V;
	cx_mat LAM;
};

struct Mesh 
{
	mat X;
	mat Y;
	mat Z;
};
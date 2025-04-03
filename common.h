#pragma once
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

const double pi = acos(-1.0);

const cx_double iI(0, 1);

struct  Smatrix
{
	cx_mat S11;
	cx_mat S12;
	cx_mat S21;
	cx_mat S22;
};

struct DataRCWA
{
	vector<cx_mat> Index;
	vec LayerPos, x, y, z;
	int ku, kv;
	double  lambda, theta, phi, n_upper, n_lower;
};

struct DataFile
{
	vector<cx_mat> Index;   // x-y-pos-lambda
	mat LayerPos, x, y, z, ku, kv, lambda, theta, phi, n_upper, n_lower;
};
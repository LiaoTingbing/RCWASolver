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

struct Source
{
	double theta;
	double phi;
	double lambda;
	double frequency;
	int ku;
	int kv;
};

struct Dev
{
	cx_mat* Index;
	vec LayerPos, x, y, z;
};


struct rcwaDATA
{
 
	clock_t t1, t2, t3, t4;
 

	size_t Nx  ;
	size_t Ny  ;
	size_t Ncx  ;
	size_t Ncy  ;
	size_t Nc  ;

	double Lx  ;
	double Ly  ;
	double thetaArc ;
	double phiArc  ;
	vec incident_direction;
	double k0  ;

 

	cx_double ER_inc ;
	cx_double UR_inc  ;
	cx_vec k_inc  ;
	cx_double kx_inc ;
	cx_double ky_inc  ;
	cx_double kz_inc  ;

 

	double Tx  ;
	double Ty  ;
 



	vec m  ;
	vec n  ;

	 size_t M ;
	 size_t N  ;
	 size_t Nh  ;


	cx_mat kx_mn;
	cx_mat ky_mn;

  
	cx_mat Kx  ;
	cx_mat Ky ;

 
	cx_mat I ;
	cx_mat II ;
	cx_mat Z ;
	cx_mat ZZ ;
 
};

struct rcwaPara
{
	int ku;
	int kv;

};
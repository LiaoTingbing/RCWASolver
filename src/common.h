#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>


//using std::cout;
//using std::endl;
//using std::string;
//using std::ifstream;
//using arma::cx_mat;
//using arma::mat;
//using arma::colvec;
using namespace std;
using namespace arma;


const double pi = 3.1415926;

const cx_double iI(0, 1);

//struct DataIn
//{
//	int i;
//	//cx_mat Index;
//	//colvec Lnode;
//};


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
	cx_mat L;
};
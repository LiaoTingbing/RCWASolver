#include "userFunc.h"

void loadTXT(double& s, string filename)
{
	mat A;
	A.load(filename);
	s = A(0);
}

void loadTXT(int& s, string filename)
{
	mat A;
	A.load(filename);
	s = int(A(0));
}

cx_mat MatrixConnect(cx_mat   A, cx_mat   B, cx_mat    C, cx_mat    D)
{
	return join_cols(join_rows(A, B), join_rows(C, D));
}

Smatrix SconnectRight(Smatrix& G, Smatrix& S)
{
	cx_mat Identity(S.S11.n_rows, S.S11.n_cols, fill::eye);
	cx_mat D = G.S12 * inv(Identity - S.S11 * G.S22);
	cx_mat F = S.S21 * inv(Identity - G.S22 * S.S11);
	G.S11 = G.S11 + D * S.S11 * G.S21;
	G.S12 = D * S.S12;
	G.S21 = F * G.S21;
	G.S22 = S.S22 + F * G.S22 * S.S12;
	return G;
}

 
 

 


 

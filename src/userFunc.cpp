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
	s = int( A(0) );
}

cx_mat matrix_to_vector_row(cx_mat& K)
{
	size_t MN = K.n_rows * K.n_cols;
	size_t col = 0;
	size_t row = 0;
	cx_mat R(MN, MN, fill::zeros);
	for (size_t i = 0;i < MN;i++)
	{
		row = i / K.n_cols;
		col = i % K.n_cols;
		R(i, i) = K(row, col);
	}
	return R ;
}

cx_mat Matrix_merging(cx_mat  A, cx_mat  B, cx_mat  C, cx_mat  D)
{
	return join_cols(join_rows(A, B), join_rows(C, D));
}

WVLmatrix solution_in_Homogeneous_Layers(cx_mat& Kx, cx_mat& Ky, cx_double& ER, cx_double& UR)
{
	size_t Nh = Kx.n_rows;
	cx_mat I(Nh, Nh, fill::eye);
	cx_mat II(2 * Nh, 2 * Nh, fill::eye);
	cx_mat Z(Nh, Nh, fill::zeros);
	cx_mat ZZ(Nh * 2, Nh * 2, fill::zeros);

	cx_mat P= 1.0 / ER * Matrix_merging(Kx * Kx, UR * ER * I - Kx * Kx, Ky * Ky - UR * ER * I, -Ky * Kx);
	cx_mat Q = ER / UR * P;

	cx_mat W = II;
	cx_mat Kz =conj( sqrt(conj(UR) * conj(ER) * I - Kx * Kx - Ky * Ky) );
	cx_mat LAM = Matrix_merging(iI * Kz , Z ,Z ,iI*Kz );
	cx_mat V = Q * inv(LAM);
	return WVLmatrix{ W , V , LAM };
}

WVLmatrix Gap(cx_mat& Kx, cx_mat& Ky, cx_mat& I, cx_mat& Z)
{
	cx_mat Kz = sqrt(I - Kx * Kx - Ky * Ky);
	cx_mat Q = Matrix_merging( Kx * Ky ,  I - Kx * Kx ,  Ky * Ky - I ,  - Kx * Ky );
	cx_mat W0 = join_cols(join_rows(I, Z), join_rows(Z, I));
	cx_mat LAM = join_cols(join_rows(iI*Kz, Z), join_rows(Z, iI*Kz));
	cx_mat V0 = Q * inv(LAM);
	return WVLmatrix{ W0 , V0 };
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
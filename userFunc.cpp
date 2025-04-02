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

WVLmatrix solution_in_Homogeneous_Layers(
	cx_mat& Kx, cx_mat& Ky, cx_double ER, cx_double UR)
{
	size_t Nh = Kx.n_rows;
	cx_mat I(Nh, Nh, fill::eye);
	cx_mat II(2 * Nh, 2 * Nh, fill::eye);
	cx_mat Z(Nh, Nh, fill::zeros);
	cx_mat ZZ(Nh * 2, Nh * 2, fill::zeros);

	cx_mat P = 1.0 / ER * MatrixConnect(Kx * Ky, UR * ER * I - Kx * Kx, Ky * Ky - UR * ER * I, -Ky * Kx);
	cx_mat Q = ER / UR * P;

	cx_mat W = II;
	cx_mat Kz = conj(sqrt(conj(UR) * conj(ER) * I - Kx * Kx - Ky * Ky));
	cx_mat LAM = MatrixConnect(iI * Kz, Z, Z, iI * Kz);
	cx_mat V = Q * inv(LAM);
	WVLmatrix R;
	R.W = W;
	R.V = V;
	R.LAM = LAM;
	return R;
}

WVLmatrix Gap(cx_mat& Kx, cx_mat& Ky, cx_mat& I, cx_mat& Z)
{
	//Z.diag().print();
	cx_mat Kz = sqrt(I - Kx * Kx - Ky * Ky).t();
	cx_mat Q = MatrixConnect(Kx * Ky, I - Kx * Kx, Ky * Ky - I, -Kx * Ky);
	cx_mat W0 = MatrixConnect(I,Z,Z,I);
	cx_mat LAM = MatrixConnect(iI * Kz, Z ,  Z, iI * Kz);
	cx_mat V0 = Q * inv(LAM);
	WVLmatrix R;
	R.W = W0;
	R.V = V0;
	//V0.diag().print();
	return R;
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

Mesh ndgrid (vec x, vec y)
{
	size_t nx = x.size();
	size_t ny = y.size();

	mat X(nx, ny, fill::zeros);
	mat Y(nx, ny, fill::zeros);

	for (size_t i = 0;i < ny;i++)
	{
		X.col(i) = x;
	}
	for (size_t i = 0;i < nx;i++)
	{
		Y.row(i) = y.st();
	}
	return { X, Y };
}

cx_mat Liu_conv_acc(
	cx_mat& Kx, cx_mat& Ky, cx_mat& Z, size_t Nh,
	size_t Ncx, size_t Ncy, cx_mat& ER, double k0, mat X, mat Y)
{
	cx_mat GG = join_rows(Kx.diag(), Ky.diag());
	cx_mat B = Z;

	//GG.print();


	cx_rowvec GM, GN, GGG;

	cx_mat A;
	for (size_t i = 0;i < Nh;i++)
	{
		for (size_t j = i;j < Nh;j++)
		{
			GM = GG(i, 0, size(1, 2));
			GN = GG(j, 0, size(1, 2));
			GGG = GM - GN;
			A = 1.0 / Ncx / Ncy * ER % exp(-iI *
				(GGG(0) * k0 * X + GGG(1) * k0 * Y));
			B(i, j) = accu(A);
		}
	}

	cx_mat D = diagmat(B.diag());
	B = B.st() + B - D;
	B = B.st();
	B.col(0).print();
	return  B;
}

WVLmatrix solution_in_inHomogeneous_Layers(
	cx_mat& Kx, cx_mat& Ky, cx_mat& ERC, cx_mat& URC)
{
	cx_mat Pi = MatrixConnect(Kx * inv(ERC) * Ky, URC - Kx * inv(ERC) * Kx,
		Ky * inv(ERC) * Ky - URC, -Ky * inv(ERC) * Kx);
	cx_mat Qi = MatrixConnect(Kx * inv(URC) * Ky, ERC - Kx * inv(URC) * Kx,
		Ky * inv(URC) * Ky - ERC, -Ky * inv(URC) * Kx);
	cx_mat Oi = Pi * Qi;

	//cx_mat OI = full(Oi);
	//cx_mat [Wi, LAMi2] = eig(OI); // % [V, D] = eig(A) 返回特征值的对角矩阵 D 和矩阵 V，其列是对应的右特征向量，使得 A * V = V * D。
	cx_vec LAMi2;
	cx_mat Wi;
	eig_gen(LAMi2, Wi, Oi);
	cx_vec LAMi = sqrt(LAMi2);
	//% Vi = Qi * Wi / LAMi;

	//%%% 稳定性排列

	/*cx_mat DV = imag(LAMi) > 0;
	cx_mat Wi =  Wi(:, DV), Wi(:, ~DV) ;
	cx_mat LAMi = diag([LAMi(DV); LAMi(~DV)]);*/


	cx_mat Vi = Qi * Wi * inv ( diagmat(  LAMi ));

	return { Wi , Vi , diagmat(LAMi) };
}

cx_mat invD(const cx_mat  & In)
{
	cx_mat Out(In.n_rows, In.n_cols);
	for (size_t i = 0; i < In.n_rows; i++)
	{
		Out(i, i) = 1.0 / In(i, i);
	}
	return Out;
}

mat invD(const mat& In)
{
	mat Out(In.n_rows, In.n_cols);
	for (size_t i = 0; i < In.n_rows; i++)
	{
		Out(i, i) = 1.0 / In(i, i);
	}
	return Out;
}

 


 

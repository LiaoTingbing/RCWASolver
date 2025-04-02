#include "convulationMatrix.h"


cx_mat Convulation_Matrix(cx_mat& ER, vec& m, vec& n)
{

	size_t M = m.size();
	size_t N = n.size();

	size_t Nh = M * N;


	size_t Ncx = ER.n_rows;
	size_t Ncy = ER.n_cols;

	size_t m0 = 1 + Ncx / 2 ;
	size_t n0 = 1 + Ncy / 2; 

	cx_mat ft = fftMat(ER);
	cx_mat A = fftshift  ( ft ) / (Ncx * Ncy );

	cx_mat ERC(Nh, Nh);

	size_t row , col,mfft,nfft;
	for (size_t nrow = 0; nrow < N; nrow++)
	{
		for (size_t mrow = 0; mrow < M; mrow++)
		{
			row = nrow  * M + mrow;
			
			for (size_t ncol = 0; ncol < N; ncol++)
			{
				for (size_t mcol = 0; mcol < M; mcol++)
				{
					col = ncol * M + mcol;

					mfft = m(mrow) - m(mcol);
					nfft = n(nrow) - n(ncol);
					ERC(row, col) = A(m0 + mfft-1, n0 + nfft-1);
				}

			}
		}

	}

	return ERC;
}

void Convulation_Matrix_()
{
	size_t m = 20;
	size_t n = m;
	size_t N = m * n;
	cx_mat A(m, n);
	for (size_t i = 0;i < m * n;i++)
	{
		A(i) = double(i) + double(i) * iI;
	}
	A.print();

	vec x = linspace(-2, 2, 5);

	Convulation_Matrix(A, x, x).col(0).print();

}

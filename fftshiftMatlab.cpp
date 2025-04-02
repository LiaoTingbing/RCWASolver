#include "fftshiftMatlab.h"
 

cx_mat fftshift(cx_mat& In)
{
	size_t m = In.n_rows;
	size_t n = In.n_cols;
 
	size_t md = round (m / 2.0 );
	size_t nd =round ( n / 2.0);

	// 行交换
	cx_mat Temp = join_cols(In.rows(md, m - 1), In.rows(0, md - 1));
	//Temp.print();
	// 列交换
	cx_mat Out = join_rows(Temp.cols(nd, n - 1), Temp.cols(0, nd - 1));

	return Out;
}

void shiftfft_()
{
	int m = 3, n = 4;
	cx_mat A(m, n);
	for (size_t i = 0; i < m*n; i++)
	{
		A(i) = i;
	}
	A.print();

	cout << endl;

	fftshift(A).print();
	
}

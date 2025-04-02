#pragma once
#include "fftMatlab.h"


cx_mat fftMat(cx_mat& In)
{
	size_t m = In.n_rows;
	size_t n = In.n_cols;
	size_t N = m * n;

	fftw_complex* in = new fftw_complex[N];
	fftw_complex* out = new fftw_complex[N];
	fftw_plan p;

	for (size_t i = 0;i < N; i++)
	{
		in[i][0] = real(In(i));
		in[i][1] = imag(In(i));
	}

	p = fftw_plan_dft_2d(m, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	cx_mat Out(m, n);

	for (size_t i = 0;i < N; i++)
	{
		Out(i) = out[i][0] + iI * out[i][1];
	}

	fftw_destroy_plan(p);

	delete[] in, out;

	return Out;
}
 
void Test() //²âÊÔFFTW
{

	//using namespace fftwpp;
	using namespace arma;
	size_t m = 4;
	size_t n = m;
	size_t N = m * n;
	cx_mat A(m, n);
	for (size_t i = 0;i < m * n;i++)
	{
		A(i) =  double(i) + double(i) * iI;
	}
	A.print();

 
	fftw_complex* in = new fftw_complex[m * n];
	fftw_complex* out = new fftw_complex[m * n];
	fftw_plan p;
 
	for (size_t i = 0;i < m * n; i++)
	{
		in[i][0] = real(A(i));
		in[i][1] = imag(A(i));
	}
	

	p = fftw_plan_dft_2d(m,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);

	for (size_t i = 0;i < m * n; i++)
	{
		A(i) = out[i][0] + iI * out[i][1];
	}

	cout << endl;
	A.print();
	
	fftw_destroy_plan(p);

	delete[] in, out;

}

void fft_()   //²âÊÔ
{
	size_t m = 4;
	size_t n = m;
	size_t N = m * n;
	cx_mat A(m, n, fill::zeros);
	for (size_t i = 0;i < m * n;i++)
	{
		A(i) = double(i) + double(i) * iI;
	}
	A.print();

	cout << endl;
	fftMat(A).print();

	A.save("test.txt", arma::raw_ascii);

}

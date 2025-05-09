#pragma once
#include "rcwa.h"
  
cx_mat MatrixConnect(cx_mat  A, cx_mat B, cx_mat  C, cx_mat  D);

Smatrix SconnectRight(Smatrix& G, Smatrix& S);

//template<class T>
//T invS(const T& A) {
//	T uni = eye<T>(size(A));
//	return solve(A, uni);
//}




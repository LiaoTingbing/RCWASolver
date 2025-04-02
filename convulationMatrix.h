#pragma once

#include "common.h"
#include "fftMatlab.h"
#include "fftshiftMatlab.h"

cx_mat Convulation_Matrix(cx_mat& ER, vec& m, vec & n);

void Convulation_Matrix_();
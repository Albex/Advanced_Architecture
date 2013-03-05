//============================================================================
// Name        : SVD.cpp
// Author      : George Rokos
// Description : SVD solver for 2x2 linear systems - implementation
//============================================================================

#include <math.h>

#include "SVD2x2.hpp"

/*
 * Calculates the eigenvalues of a 2x2 matrix A
 * eigenvalues[0] contains the largest eigenvalue
 * eigenvalues[1] contains the smallest eigenvalue
 */

void svd_solve_2x2(const double A[4], double p[2], const double q[2]){
  double const det = A[ 0 ] * A[ 3 ] - A[ 1 ] * A[ 2 ];
  P[ 0 ] = ( q[ 0 ] * A[ 3 ] - q[ 1 ] * A[ 2 ] ) / det;
  P[ 1 ] = ( A[ 0 ] * q[ 1 ] - A[ 2 ] * q[ 0 ] ) / det;
}

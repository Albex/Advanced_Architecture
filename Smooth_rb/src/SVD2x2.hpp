//============================================================================
// Name        : SVD.cpp
// Author      : George Rokos
// Description : SVD solver for 2x2 linear systems - definitions
//============================================================================

#ifndef SVD2X2_HPP_
#define SVD2X2_HPP_

static inline void
svd_solve_2x2( double const A[ 4 ], double P[ 2 ], double const q[ 2 ] ){

  double const det = A[ 0 ] * A[ 3 ] - A[ 1 ] * A[ 2 ];
  P[ 0 ] = ( q[ 0 ] * A[ 3 ] - q[ 1 ] * A[ 2 ] ) / det;
  P[ 1 ] = ( A[ 0 ] * q[ 1 ] - A[ 2 ] * q[ 0 ] ) / det;
}

#endif /* SVD_HPP_ */

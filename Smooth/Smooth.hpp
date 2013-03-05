//============================================================================
// Name        : Smooth.hpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel prototype
//============================================================================

#ifndef SMOOTH_HPP_
#define SMOOTH_HPP_

#include "Mesh.hpp"

void smooth(Mesh * __restrict__ mesh, size_t niter);

#endif /* SMOOTH_HPP_ */

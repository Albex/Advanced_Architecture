#ifndef COLORGRAPH_HPP_INCLUDED
#define COLORGRAPH_HPP_INCLUDED

#include <SmoothConfig.hpp>
#include <vector>
#include <blitz/array.h>

void
Coloring(
	std::vector< std::vector< uint32_t > > const & graph,
	blitz::Array< int, 1 > & colors
) noexcept;

#endif
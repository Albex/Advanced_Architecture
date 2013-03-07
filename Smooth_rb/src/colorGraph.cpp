#include <colorGraph.hpp>

INLINE static int
get_colors( std::vector< uint32_t > const & graph_node )
{
    uint32_t size = graph_node.size( );
    uint32_t C = 0;
    while ( true ) {
        uint32_t i;
        for ( i = 0; i < size; ++i ) {
            if ( C == graph_node[ i ] ) {
                ++C;
                break;
            }
        }
        if ( i == size ) {
            return C;
        }
    }
}

std::vector< int >
Coloring( std::vector< std::vector< uint32_t > > const & graph )
{
    uint32_t size = graph.size( );
    std::vector< int > colors( size, -1 );
    for ( uint32_t i = 0; i < size; ++i ) {
        colors[ i ] = get_colors( graph[ i ] );
    }
    return colors;
}
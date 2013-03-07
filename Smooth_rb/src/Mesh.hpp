//============================================================================
// Name        : Mesh.hpp
// Author      : George Rokos
// Description : Mesh description
//============================================================================

#ifndef MESH_HPP_
#define MESH_HPP_

#include <SmoothConfig.hpp>
#include <string>
#include <vector>
#include <blitz/array.h>

struct Quality {
    real mean;
    real min;
    real rms;
};

class Mesh{
public:
    // Constructor
     Mesh( void ) = delete;
     Mesh( std::string const & filename) noexcept;
    ~Mesh( void ) noexcept;

    INLINE bool isSurfaceNode   ( uint32_t ) const noexcept;
    INLINE bool isCornerNode    ( uint32_t ) const noexcept;
    INLINE real element_area    ( uint32_t ) const noexcept;
    INLINE real element_quality ( int, int ) const noexcept;
    INLINE real element_quality ( uint32_t ) const noexcept;
    Quality     get_mesh_quality( void     ) const noexcept;

    void smooth( uint32_t niter ) noexcept;

private:

    void create_adjacency( void ) noexcept;
    void find_surface    ( void ) noexcept;
    void set_orientation ( void ) noexcept;

    uint32_t NNodes;    // Number of mesh vertices.
    uint32_t NElements; // Number of mesh elements.
    int max_color;

    blitz::Array< int, 1 > colors;

    // Element eid is comprised of the vertices
    // ENList[3*eid], ENList[3*eid+1] and ENList[3*eid+2].
    blitz::Array< uint32_t, 1 > ENList;
    /* If vid is on the surface, the normal vector
    * (normals[2*vid],normals[2*vid+1] =
    *                            = (0.0,1.0) if vid is on the top surface
    *                            = (0.0,-1.0) if vid is on the bottom surface
    *                            = (1.0,0.0) if vid is on the right surface
    *                            = (-1.0,0.0) if vid is on the left surface
    * For all other vertices, the normal vector is (0.0,0.0).
    */
    blitz::Array< real, 1 > normals;


    // Vertex vid has coordinates x=coords[2*vid] and y=coords[2*vid+1].
    std::vector< real > coords;

    // The metric tensor at vertex vid is M_00 = metric[3*vid],
    //                                    M_01 = M_10 = metric[3*vid+1] and
    //                                    M_11 = metric[3*vid+2].
    std::vector< real > metric;

    // For every vertex i, NNList[i] contains the IDs of all adjacent vertices.
    std::vector< std::vector< uint32_t > > NNList;

    // For every vertex i, NEList[i] contains the IDs of all adjacent elements.
    std::vector< std::vector< uint32_t > > NEList;

    int orientation;
};

bool
Mesh::isSurfaceNode( uint32_t vid ) const noexcept {
  return NEList[ vid ].size( ) < NNList[ vid ].size( );
}

bool
Mesh::isCornerNode( uint32_t vid ) const noexcept {
  return
    std::abs( normals( 2 * vid     ) ) == 1.0 &&
    std::abs( normals( 2 * vid + 1 ) ) == 1.0;
}

/* This function evaluates the quality of an element, based on the 2D quality
 * functional proposed by Lipnikov et. al.. The description for the functional
 * is taken from: Yu. V. Vasileskii and K. N. Lipnikov, An Adaptive Algorithm
 * for Quasioptimal Mesh Generation, Computational Mathematics and Mathematical
 * Physics, Vol. 39, No. 9, 1999, pp. 1468 - 1486.
 */
real Mesh::element_quality( int vid, int it ) const noexcept {

    uint32_t eid = NEList[ vid ][ it ];
    return element_quality( eid );
}

real Mesh::element_quality( uint32_t eid ) const noexcept {

    uint32_t const *n = &ENList( 3 * eid );

    // Pointers to the coordinates of each vertex
    real const *c0 = &coords[ 2 * n[ 0 ] ];
    real const *c1 = &coords[ 2 * n[ 1 ] ];
    real const *c2 = &coords[ 2 * n[ 2 ] ];

    // Pointers to the metric tensor at each vertex
    real const *m0 = &metric[ 3 * n[ 0 ] ];
    real const *m1 = &metric[ 3 * n[ 1 ] ];
    real const *m2 = &metric[ 3 * n[ 2 ] ];

    // Metric tensor averaged over the element
    real m00 = ( m0[ 0 ] + m1[ 0 ] + m2[ 0 ] ) / 3;
    real m01 = ( m0[ 1 ] + m1[ 1 ] + m2[ 1 ] ) / 3;
    real m11 = ( m0[ 2 ] + m1[ 2 ] + m2[ 2 ] ) / 3;

    // l is the length of the perimeter, measured in metric space
    real s1l1 = ( c0[ 1 ] - c1[ 1 ] );
    real s2l1 = ( c0[ 0 ] - c1[ 0 ] );
    real s1l2 = ( c0[ 1 ] - c2[ 1 ] );
    real s2l2 = ( c0[ 0 ] - c2[ 0 ] );
    real s1l3 = ( c2[ 1 ] - c1[ 1 ] );
    real s2l3 = ( c2[ 0 ] - c1[ 0 ] );

    real l =
        std::sqrt(
            s1l1 * ( s1l1 * m11 + s2l1 * m01 ) +
            s2l1 * ( s1l1 * m01 + s2l1 * m00 )
        ) +
        std::sqrt(
            s1l2 * ( s1l2 * m11 + s2l2 * m01 ) +
            s2l2 * ( s1l2 * m01 + s2l2 * m00 )
        ) +
        std::sqrt(
            s1l3 * ( s1l3 * m11 + s2l3 * m01 ) +
            s2l3 * ( s1l3 * m01 + s2l3 * m00 )
        );

    // Area in physical space
    real a = element_area( eid );

    // Area in metric space
    real a_m = a * sqrt( m00 * m11 - m01 * m01 );

    // Function
    real f = std::min( l / 3.0, 3.0 / l );
    real F = f * ( 2.0 - f );

    // This is the 2D Lipnikov functional.
    real constexpr cte = 12.0 * sqrt( 3.0 );
    real quality = cte * a_m / ( l * l );

    return quality * F * F * F;
}

/* Element area in physical (Euclidean) space. Recall that the area of a
 * triangle ABC is calculated as area=0.5*(AB⋅AC), i.e. half the inner product
 * of two of the element's edges (e.g. AB and AC). The result is corrected by
 * the orientation factor ±1.0, so that the area is always a positive number.
 */
real
Mesh::element_area( uint32_t eid ) const noexcept {

    uint32_t const eid_off = 3 * eid;
    uint32_t const c0_off = 2 * ENList( eid_off     );
    uint32_t const c1_off = 2 * ENList( eid_off + 1 );
    uint32_t const c2_off = 2 * ENList( eid_off + 2 );

    return orientation * static_cast< real >( 0.5 ) * (
        ( coords[ c0_off + 1 ] - coords[ c2_off + 1 ] ) *
        ( coords[ c0_off     ] - coords[ c1_off     ] ) -
        ( coords[ c0_off + 1 ] - coords[ c1_off + 1 ] ) *
        ( coords[ c0_off     ] - coords[ c2_off     ] )
    );
}

#endif /* MESH_HPP_ */

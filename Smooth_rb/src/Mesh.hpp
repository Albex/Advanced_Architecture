//============================================================================
// Name        : Mesh.hpp
// Author      : George Rokos
// Description : Mesh description
//============================================================================

#ifndef MESH_HPP_
#define MESH_HPP_

#include <SmoothConfig.hpp>

struct Quality {
    real mean;
    real min;
    real rms;
};

class Mesh{
public:
    // Constructor
     Mesh( void ) = delete;
     Mesh( std::string const & filename);
    ~Mesh( void );

    INLINE bool isSurfaceNode   ( uint32_t ) const;
    INLINE bool isCornerNode    ( uint32_t ) const;
    INLINE real element_area    ( uint32_t ) const;
    INLINE real element_quality ( int, int ) const;
    INLINE real element_quality ( uint32_t ) const;
    Quality     get_mesh_quality( void     ) const;

    void smooth( uint32_t niter );

private:

    void create_adjacency( void );
    void find_surface    ( void );
    void set_orientation ( void );

    uint32_t NNodes;    // Number of mesh vertices.
    uint32_t NElements; // Number of mesh elements.

    // Element eid is comprised of the vertices
    // ENList[3*eid], ENList[3*eid+1] and ENList[3*eid+2].
    std::vector< uint32_t > ENList;

    // Vertex vid has coordinates x=coords[2*vid] and y=coords[2*vid+1].
    std::vector< real > coords;

    // The metric tensor at vertex vid is M_00 = metric[3*vid],
    //                                    M_01 = M_10 = metric[3*vid+1] and
    //                                    M_11 = metric[3*vid+2].
    std::vector< real > metric;

    /* If vid is on the surface, the normal vector
    * (normals[2*vid],normals[2*vid+1] =
    *                            = (0.0,1.0) if vid is on the top surface
    *                            = (0.0,-1.0) if vid is on the bottom surface
    *                            = (1.0,0.0) if vid is on the right surface
    *                            = (-1.0,0.0) if vid is on the left surface
    * For all other vertices, the normal vector is (0.0,0.0).
    */
    blitz::Array< real, 1 > normals;

    // For every vertex i, NNList[i] contains the IDs of all adjacent vertices.
    std::vector< std::vector< uint32_t > > NNList;

    // For every vertex i, NEList[i] contains the IDs of all adjacent elements.
    std::vector< std::vector< uint32_t > > NEList;

    int orientation;
};

bool
Mesh::isSurfaceNode( uint32_t vid ) const {
  return NEList[ vid ].size( ) < NNList[ vid ].size( );
}

bool
Mesh::isCornerNode( uint32_t vid ) const{
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
real Mesh::element_quality( int vid, int it ) const {

    uint32_t eid = NEList[ vid ][ it ];
    return element_quality( eid );
}

real Mesh::element_quality( uint32_t eid ) const{

    uint32_t const eid_off = 3 * eid;

    // Pointers to the coordinates of each vertex
    const real * __restrict__ c[ 3 ] = {
        &coords[ 2 * ENList[ eid_off     ] ],
        &coords[ 2 * ENList[ eid_off + 1 ] ],
        &coords[ 2 * ENList[ eid_off + 2 ] ]
    };

    // uint32_t const met_off[ 3 ] = {
    //     3 * ENList[ eid_off     ],
    //     3 * ENList[ eid_off + 1 ],
    //     3 * ENList[ eid_off + 2 ]
    // };

    // // Metric tensor averaged over the element
    real m00;
    real m01;
    real m11;
    // for ( int i = 0; i < 3; ++i ) {
    //     m00 += metric( met_off[ i ]     );
    //     m01 += metric( met_off[ i ] + 1 );
    //     m11 += metric( met_off[ i ] + 2 );
    // }

    const real * __restrict__ m[ 3 ] = {
        &metric[ 3 * ENList[ eid_off     ] ],
        &metric[ 3 * ENList[ eid_off + 1 ] ],
        &metric[ 3 * ENList[ eid_off + 2 ] ]
    };

    for ( int i = 0; i < 3; ++i ) {
        m00 += m[ i ][ 0 ];
        m01 += m[ i ][ 1 ];
        m11 += m[ i ][ 2 ];
    }

    // l is the length of the perimeter, measured in metric space
    real s1l1 = ( c[ 0 ][ 1 ] - c[ 1 ][ 1 ] );
    real s2l1 = ( c[ 0 ][ 0 ] - c[ 1 ][ 0 ] );
    real s1l2 = ( c[ 0 ][ 1 ] - c[ 2 ][ 1 ] );
    real s2l2 = ( c[ 0 ][ 0 ] - c[ 2 ][ 0 ] );
    real s1l3 = ( c[ 2 ][ 1 ] - c[ 1 ][ 1 ] );
    real s2l3 = ( c[ 2 ][ 0 ] - c[ 1 ][ 0 ] );

    real s2l3b = s2l3 * ( s1l3 * m01 + s2l3 * m00 );
    real s2l2b = s2l2 * ( s1l2 * m01 + s2l2 * m00 );
    real s2l1b = s2l1 * ( s1l1 * m01 + s2l1 * m00 );
    real s1l3b = s1l3 * ( s1l3 * m11 + s2l3 * m01 );
    real s1l2b = s1l2 * ( s1l2 * m11 + s2l2 * m01 );
    real s1l1b = s1l1 * ( s1l1 * m11 + s2l1 * m01 );

    real l1 = static_cast< real >( std::sqrt( ( s1l1b + s2l1b ) ) / 3 );
    real l2 = static_cast< real >( std::sqrt( ( s1l2b + s2l2b ) ) / 3 );
    real l3 = static_cast< real >( std::sqrt( ( s1l3b + s2l3b ) ) / 3 );

    real l = l1 + l2 + l3;

    // Area in physical space
    real a = element_area( eid );

    real f = static_cast< real >( std::min( l / 3, 3 / l ) );

    // Area in metric space
    real a_m = a * static_cast< real >( std::sqrt( m00 * m11 - m01 * m01 ) );

    // Function
    real F = f * ( 2 - f );

    // This is the 2D Lipnikov functional.
    real constexpr sqrt3 = static_cast< real >( std::sqrt( 3 ) );
    real quality = 12 * sqrt3 * a_m / ( l * l);

    return quality * F * F * F;
}

/* Element area in physical (Euclidean) space. Recall that the area of a
 * triangle ABC is calculated as area=0.5*(AB⋅AC), i.e. half the inner product
 * of two of the element's edges (e.g. AB and AC). The result is corrected by
 * the orientation factor ±1.0, so that the area is always a positive number.
 */
real
Mesh::element_area( uint32_t eid ) const{

    uint32_t const eid_off = 3 * eid;
    uint32_t const c0_off = 2 * ENList[ eid_off     ];
    uint32_t const c1_off = 2 * ENList[ eid_off + 1 ];
    uint32_t const c2_off = 2 * ENList[ eid_off + 2 ];

    return orientation * static_cast< real >( 0.5 ) * (
        ( coords[ c0_off + 1 ] - coords[ c2_off + 1 ] ) *
        ( coords[ c0_off     ] - coords[ c1_off     ] ) -
        ( coords[ c0_off + 1 ] - coords[ c1_off + 1 ] ) *
        ( coords[ c0_off     ] - coords[ c2_off     ] )
    );
}

#endif /* MESH_HPP_ */

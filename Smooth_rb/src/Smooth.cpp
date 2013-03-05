//============================================================================
// Name        : Smooth.cpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel - Smart Laplacian variant
//============================================================================

#include <algorithm>
#include <cmath>

#include "SVD2x2.hpp"
#include "Smooth.hpp"

void smooth( Mesh & mesh, size_t niter ) {
    // For the specified number of iterations, loop over all mesh vertices.
    uint32_t max = 0;
    for ( uint32_t i = 0; i < mesh.NEList.size( ); ++i ) {
        if ( mesh.NEList[ i ].size( ) > max ) {
            max = mesh.NEList[ i ].size( );
        }
    }
    double * buffer = new double[ max ];
    if ( buffer == nullptr ) {
        return;
    }
    for( size_t iter = 0; iter < niter; ++iter ) {
        for( size_t vid = 0; vid < mesh.NNodes; ++vid ) {

            // If this is a corner node, it cannot be moved.
            if( mesh.isCornerNode( vid ) ) continue;

            // Find the quality of the worst element adjacent to vid
            uint32_t const size = mesh.NEList[ vid ].size( );
            for( uint32_t it = 0; it < size; ++it ) {
                buffer[ it ] = mesh.element_quality( vid, it );
            }
            double worst_q = 1.0;
            for( uint32_t it = 0; it < size; ++it ) {
                worst_q = std::min( worst_q, buffer[ it ] );
            }

            /* Find the barycentre (centre of mass) of the cavity. A cavity is
            * defined as the set containing vid and all its adjacent vertices and
            * triangles. Since we work on metric space, all lengths have to measured
            * using the metric. The metric tensor is a 2x2 symmetric matrix which
            * encodes the ideal length and orientation of an edge containing vid. As
            * an edge is defined by two vertices, we calculate the edge length using
            * the value of the metric in the middle of the edge, i.e. the average of
            * the two metric tensors of the vertices defining the edge.
            */

            double x0 = mesh.coords[ 2 * vid     ];
            double y0 = mesh.coords[ 2 * vid + 1 ];

            double A[ 4 ]; std::memset( A, 0, 4 * sizeof( double ) );
            double q[ 2 ]; std::memset( A, 0, 2 * sizeof( double ) );

            // Iterate over all edges and assemble matrices A and q.
            uint32_t const vid_offset  = 3 * vid;
            uint32_t const size_nnlist = mesh.NNList[ vid ].size( );
            for( uint32_t it = 0; it < size_nnlist; ++it ) {

                size_t   const il         = mesh.NNList[ vid ][ it ];
                uint32_t const met_offset = 3 * il;

                // Find the metric in the middle of the edge.
                double ml00 = 0.5 * ( mesh.metric[ vid_offset     ] +
                    mesh.metric[ met_offset + 0 ] );
                double ml01 = 0.5 * ( mesh.metric[ vid_offset + 1 ] +
                    mesh.metric[ met_offset + 1 ] );
                double ml11 = 0.5 * ( mesh.metric[ vid_offset + 2 ] +
                    mesh.metric[ met_offset + 2 ] );

                double x = mesh.coords[ 2 * il     ] - x0;
                double y = mesh.coords[ 2 * il + 1 ] - y0;

                // Calculate and accumulate the contribution of
                // this vertex to the barycentre of the cavity.
                q[ 0 ] += ( ml00 * x + ml01 * y );
                q[ 1 ] += ( ml01 * x + ml11 * y );

                A[ 0 ] += ml00;
                A[ 1 ] += ml01;
                A[ 3 ] += ml11;
            }

            // The metric tensor is symmetric, i.e. ml01=ml10, so A[2]=A[1].
            A[ 2 ] = A[ 1 ];

            // Displacement vector for vid
            double p[ 2 ];

            /* The displacement p for vid is found by solving the linear system:
            * ┌─       ─┐   ┌    ┐   ┌    ┐
            * │A[0] A[1]│   │p[0]│   │q[0]│
            * │         │ x │    │ = │    │
            * │A[2] A[3]│   │p[1]│   │q[0]│
            * └─       ─┘   └    ┘   └    ┘
            */
            svd_solve_2x2( A, p, q );

            /* If this is a surface vertex, restrict the displacement
            * to the surface. The new displacement is the projection
            * of the old displacement on the surface.
            */
            if( mesh.isSurfaceNode( vid ) ){
                p[ 0 ] -= p[ 0 ] * std::abs( mesh.normals[ 2 * vid    ] );
                p[ 1 ] -= p[ 1 ] * std::abs( mesh.normals[ 2 * vid + 1] );
            }

            // Update the coordinates
            mesh.coords[ 2 * vid     ] += p[ 0 ];
            mesh.coords[ 2 * vid + 1 ] += p[ 1 ];

            /************************************************************************
            * At this point we must also interpolate the metric tensors from all   *
            * neighbouring vertices in order to calculate the new value of vid's   *
            * metric tensor at the new location. This is a quite complex procedure *
            * and has been omitted for simplicity of the exercise. A vertex will   *
            * always use its original metric tensor, no matter whether it has been *
            * relocated or not.                                                    *
            ************************************************************************/

            /* Find the quality of the worst element after smoothing. If an element
            * of the cavity was inverted, i.e. if vid was relocated outside the
            * interior convex hull of the cavity, then the calculated area of that
            * element will be negative and mesh.element_quality() will return a
            * negative number. In such a case, the smoothing operation has to be
            * rejected.
            */
            for( uint32_t it = 0; it < size; ++it ) {
                buffer[ it ] = mesh.element_quality( vid, it );
            }
            double new_worst_q = 1.0;
            for( uint32_t it = 0; it < size; ++it ) {
                new_worst_q = std::min( new_worst_q, buffer[ it ] );
            }
            /* If quality is worse than before, either because of element inversion
            * or just because relocating vid to the barycentre of the cavity does
            * not improve quality, revert the changes.
            */
            if( new_worst_q < worst_q ){
                mesh.coords[ 2 * vid     ] -= p[ 0 ];
                mesh.coords[ 2 * vid + 1 ] -= p[ 1 ];
            }
        }
    }
}

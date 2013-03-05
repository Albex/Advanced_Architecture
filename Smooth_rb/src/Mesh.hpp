//============================================================================
// Name        : Mesh.hpp
// Author      : George Rokos
// Description : Mesh description
//============================================================================

#ifndef MESH_HPP_
#define MESH_HPP_

#include <cstddef>
#include <set>
#include <vector>
#include <cmath>

struct Quality{
  double mean;
  double min;
  double rms;
};

class Mesh{
public:
  // Constructor
  Mesh(const char * __restrict__ filename);

  size_t NNodes;    // Number of mesh vertices.
  size_t NElements; // Number of mesh elements.

  // Element eid is comprised of the vertices
  // ENList[3*eid], ENList[3*eid+1] and ENList[3*eid+2].
  std::vector<size_t> ENList;

  // Vertex vid has coordinates x=coords[2*vid] and y=coords[2*vid+1].
  std::vector<double> coords;

  // The metric tensor at vertex vid is M_00 = metric[3*vid],
  //                                    M_01 = M_10 = metric[3*vid+1] and
  //                                    M_11 = metric[3*vid+2].
  std::vector<double> metric;

  /* If vid is on the surface, the normal vector
   * (normals[2*vid],normals[2*vid+1] =
   *                            = (0.0,1.0) if vid is on the top surface
   *                            = (0.0,-1.0) if vid is on the bottom surface
   *                            = (1.0,0.0) if vid is on the right surface
   *                            = (-1.0,0.0) if vid is on the left surface
   * For all other vertices, the normal vector is (0.0,0.0).
   */
  std::vector< double > normals;

  // For every vertex i, NNList[i] contains the IDs of all adjacent vertices.
  std::vector< std::vector< size_t > > NNList;

  // For every vertex i, NEList[i] contains the IDs of all adjacent elements.
  std::vector< std::vector< size_t > > NEList;

  inline bool isSurfaceNode(size_t vid) const;
  inline bool isCornerNode(size_t vid) const;
  double element_area(size_t eid) const;
  inline double element_quality(size_t eid) const;
  Quality get_mesh_quality() const;

private:
  void create_adjacency( );
  void find_surface( );
  void set_orientation( );

  int orientation;
};

bool
Mesh::isSurfaceNode( size_t vid ) const {
  return NEList[ vid ].size( ) < NNList[ vid ].size( );
}

bool
Mesh::isCornerNode( size_t vid ) const{
  return std::abs(
    normals[ 2 * vid] ) == 1.0 && std::abs( normals[ 2 * vid + 1 ] == 1.0
  );
}

/* This function evaluates the quality of an element, based on the 2D quality
 * functional proposed by Lipnikov et. al.. The description for the functional
 * is taken from: Yu. V. Vasileskii and K. N. Lipnikov, An Adaptive Algorithm
 * for Quasioptimal Mesh Generation, Computational Mathematics and Mathematical
 * Physics, Vol. 39, No. 9, 1999, pp. 1468 - 1486.
 */
double Mesh::element_quality(size_t eid) const{
  const int eid_off = 3*eid;

  // Pointers to the coordinates of each vertex
  const double * __restrict__ c[ 3 ] = {
    &coords[2*ENList[eid_off  ]],
    &coords[2*ENList[eid_off+1]],
    &coords[2*ENList[eid_off+2]]
  };

  // Pointers to the metric tensor at each vertex
  const double * __restrict__ m[ 3 ] = {
    &metric[3*ENList[eid_off  ]],
    &metric[3*ENList[eid_off+1]],
    &metric[3*ENList[eid_off+2]]
  };

  // Metric tensor averaged over the element
  double m00;
  double m01;
  double m11;
  for ( int i = 0; i < 3; ++i ) {
    m00 += m[i][0];
    m01 += m[i][1];
    m11 += m[i][2];
  }
  // l is the length of the perimeter, measured in metric space

  double s1l1 = ( c[ 0 ][ 1 ] - c[ 1 ][ 1 ] );
  double s2l1 = ( c[ 0 ][ 0 ] - c[ 1 ][ 0 ] );
  double s1l2 = ( c[ 0 ][ 1 ] - c[ 2 ][ 1 ] );
  double s2l2 = ( c[ 0 ][ 0 ] - c[ 2 ][ 0 ] );
  double s1l3 = ( c[ 2 ][ 1 ] - c[ 1 ][ 1 ] );
  double s2l3 = ( c[ 2 ][ 0 ] - c[ 1 ][ 0 ] );

  double s2l3b = s2l3 * ( s1l3 * m01 + s2l3 * m00 );
  double s2l2b = s2l2 * ( s1l2 * m01 + s2l2 * m00 );
  double s2l1b = s2l1 * ( s1l1 * m01 + s2l1 * m00 );
  double s1l3b = s1l3 * ( s1l3 * m11 + s2l3 * m01 );
  double s1l2b = s1l2 * ( s1l2 * m11 + s2l2 * m01 );
  double s1l1b = s1l1 * ( s1l1 * m11 + s2l1 * m01 );

  double l1 = sqrt( ( s1l1b + s2l1b ) / 3.0 );
  double l2 = sqrt( ( s1l2b + s2l2b ) / 3.0 );
  double l3 = sqrt( ( s1l3b + s2l3b ) / 3.0 );

  double l = l1 + l2 + l3;

  // Area in physical space
  double a = element_area( eid );

  double f = std::min( l / 3.0, 3.0 / l );

  // Area in metric space
  double a_m = a * sqrt( m00 * m11 - m01 * m01 );

  // Function
  double F = f * ( 2.0 - f );

  // This is the 2D Lipnikov functional.
  double quality = 12.0 * sqrt( 3.0 ) * a_m / ( l * l);

  return quality * F * F * F;
}

#endif /* MESH_HPP_ */

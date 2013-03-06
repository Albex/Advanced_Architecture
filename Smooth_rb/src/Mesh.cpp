//============================================================================
// Name        : Mesh.cpp
// Author      : George Rokos
// Description : Mesh implementation
//============================================================================

#include <Mesh.hpp>

Mesh::Mesh( std::string const & filename){
  // Check whether the provided file exists.
  ifstream ifile(filename);
  if(!ifile){
    std::cerr << "File " << filename << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }

  vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName( filename.c_str( ) );
  reader->Update();

  vtkUnstructuredGrid * __restrict__ ug = reader->GetOutput();

  NNodes = ug->GetNumberOfPoints();
  NElements = ug->GetNumberOfCells();

  // Get the coordinates of each mesh vertex. There is no z coordinate in 2D,
  // but VTK treats 2D and 3D meshes uniformly, so we have to provide memory
  // for z as well (r[2] will always be zero and we ignore it).
  for(uint32_t i=0;i<NNodes;++i){
    double r[ 3 ];
    ug->GetPoints()->GetPoint(i, r);
    coords.push_back(r[0]);
    coords.push_back(r[1]);
  }
  assert(coords.size() == 2*NNodes);

  // Get the metric at each vertex.
  metric.resize( 3 * NNodes );
  for(uint32_t i=0;i<NNodes;++i){
    double * __restrict__ tensor = ug->GetPointData()->GetArray("Metric")->GetTuple4(i);
    metric[ 3 * i     ] = tensor[0];
    metric[ 3 * i + 1 ] = tensor[1];
    assert(tensor[1] == tensor[2]);
    metric[ 3 * i + 2 ] = tensor[3];
  }

  // Get the 3 vertices comprising each element.
  for(uint32_t i=0;i<NElements;++i){
    vtkCell * __restrict__ cell = ug->GetCell(i);
    for(int j=0;j<3;j++){
      ENList.push_back(cell->GetPointId(j));
    }
  }
  assert(ENList.size() == 3*NElements);

  reader->Delete();

  create_adjacency();
  find_surface();
  set_orientation();
}

Mesh::~Mesh( void ) {

}

void Mesh::create_adjacency(){
  NNList.resize(NNodes);
  NEList.resize(NNodes);

  std::vector< std::set< uint32_t > > temp_sets;
  temp_sets.resize( NNodes );

  for(uint32_t eid=0; eid<NElements; ++eid){
    // Get a pointer to the three vertices comprising element eid.
    const uint32_t * __restrict__ n = &ENList[3*eid];

    // For each vertex, add the other two vertices to its node-node adjacency
    // list and element eid to its node-element adjacency list.
    for(uint32_t i=0; i<3; ++i){
      NNList[n[i]].push_back(n[(i+1)%3]);
      NNList[n[i]].push_back(n[(i+2)%3]);

      temp_sets[n[i]].insert(eid);
    }
  }
  for(uint32_t i=0; i<NNodes; ++i){
    NEList[ i ].assign( temp_sets[ i ].begin( ), temp_sets[ i ].end( ) );
  }
}

void Mesh::find_surface(){
  // Initialise all normal vectors to (0.0,0.0).
  normals.resize(2*NNodes, 0.0);

  // If an edge is on the surface, then it belongs to only 1 element. We
  // traverse all edges (vid0,vid1) and for each edge we find the intersection
  // of NEList[vid0] and NEList[vid1]. If the intersection size is 1, then
  // this edge belongs to only one element, so it lies on the mesh surface.
  for(uint32_t vid=0; vid<NNodes; ++vid){
    for(std::vector<uint32_t>::const_iterator it=NNList[vid].begin();
      it!=NNList[vid].end(); ++it){
      // In order to avoid processing an edge twice, one in the
      // form of (vid0,vid1) and one in the form of (vid1,vid0),
      // an edge is processed only of vid0 < vid1.
      if(vid > *it)
        continue;

      std::set<uint32_t> intersection;
      std::set_intersection(NEList[vid].begin(), NEList[vid].end(),
          NEList[*it].begin(), NEList[*it].end(),
          std::inserter(intersection, intersection.begin()));

      if(intersection.size()==1){ // We have found a surface edge
        real x=coords[2*vid], y=coords[2*vid+1];

        // Find which surface vid and *it belong to and set the corresponding
        // coordinate of the normal vector to ±1.0. The other coordinate is
        // intentionally left intact. This way, the normal vector for corner
        // vertices will be at the end (±1.0,±1.0), which enables us to detect
        // that they are corner vertices and are not allowed to be smoothed.
        if(fabs(y-1.0) < 1E-12){// vid is on the top surface
          normals( 2 * vid + 1 ) = 1.0;
          normals( 2 * ( *it ) + 1 ) = 1.0;
        }
        else if(fabs(y) < 1E-12){// vid is on the bottom surface
          normals( 2 * vid + 1 ) = -1.0;
          normals( 2 * ( *it ) + 1 ) = -1.0;
        }
        else if(fabs(x-1.0) < 1E-12){// vid is on the right surface
          normals( 2 * vid ) = 1.0;
          normals( 2 * ( *it ) ) = 1.0;
        }
        else if(fabs(x) < 1E-12){// vid is on the left surface
          normals( 2 * vid ) = -1.0;
          normals( 2 * ( *it ) ) = -1.0;
        }
        else{
          std::cerr << "Invalid surface vertex coordinates" << std::endl;
        }
      }
    }
  }
}

/* Computing the area of an element as the inner product of two element edges
 * depends on the order in which the three vertices comprising the element have
 * been stored in ENList. Using the right-hand rule for the cross product of
 * two 2D vectors, we can find whether the three vertices define a positive or
 * negative element. If the order in ENList suggests a clockwise traversal of
 * the element, the cross product faces the negative z-axis, so the element is
 * negative and we will have to correct the calculation of its area by
 * multiplying the result by -1.0. During smoothing, after a vertex is
 * relocated, if the area of an adjacent element is found to be negative it
 * means that the element has been inverted and the new location of the vertex
 * should be discarded. The orientation is the same for all mesh elements, so
 * it is enough to calculate it for one mesh element only.
 */
void Mesh::set_orientation(){
  // Find the orientation for the first element
  const uint32_t * __restrict__ n = &ENList[0];

  // Pointers to the coordinates of each vertex
  const real *c0 = &coords[2*n[0]];
  const real *c1 = &coords[2*n[1]];
  const real *c2 = &coords[2*n[2]];

  real x1 = (c0[0] - c1[0]);
  real y1 = (c0[1] - c1[1]);

  real x2 = (c0[0] - c2[0]);
  real y2 = (c0[1] - c2[1]);

  real A = x1*y2 - x2*y1;

  if(A<0)
    orientation = -1;
  else
    orientation = 1;
}



// Finds the mean quality, averaged over all mesh elements,
// and the quality of the worst element.
Quality Mesh::get_mesh_quality() const{
  Quality q;

  real mean_q = 0.0;
  real min_q = 1.0;

  for(uint32_t i=0;i<NElements;++i){
    real ele_q = element_quality(i);

    mean_q += ele_q;
    min_q = std::min(min_q, ele_q);
  }

  q.mean = mean_q/NElements;
  q.min = min_q;

  return q;
}

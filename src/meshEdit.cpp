#include <float.h>
#include <assert.h>
#include "meshEdit.h"
#include "vsa/vsa.h"
#include "mutablePriorityQueue.h"

namespace CMU462 {

  VertexIter HalfedgeMesh::splitEdge(EdgeIter e0, bool markNew) {

    // This method should split the given edge and return an iterator to the
    // newly inserted vertex. The halfedge of this vertex should point along
    // the edge that was split, rather than the new edges.
    // Only runs on triangle meshes
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->twin();
    if (h0->face()->degree() != 3 || h1->face()->degree() != 3) {
      return h0->vertex();
    }
    // PHASE I: Collect Elements
    HalfedgeIter h2 = h0->next(), h3 = h1->next(),
        h4 = h2->next(), h5 = h3->next();
    VertexIter v0 = h0->vertex(), v1 = h1->vertex(),
        v2 = h4->vertex(), v3 = h5->vertex();
    FaceIter     f0 = h0->face(), f1 = h1->face();

    //PHASE II: Allocate new elements
    VertexIter vn = newVertex();
    vn->position = ( v0->position + v1->position ) / 2;
    HalfedgeIter hn0 = newHalfedge(), hn1 = newHalfedge(),
        hn2 = newHalfedge(), hn3 = newHalfedge(),
        hn4 = newHalfedge(), hn5 = newHalfedge();
    EdgeIter en0 = newEdge(), en1 = newEdge(), en2 = newEdge();
    FaceIter fn0 = newFace(), fn1 = newFace();

    // PHASE III: Reassign Elements
    h0->next() = hn2; h1->vertex() = vn; h2->next() = hn3; h2->face() = fn0; h3->next() = hn4;
    h5->next() = hn1; h5->face() = fn1;

    hn0->next() =  h2; hn0->twin() = hn1; hn0->vertex() = vn; hn0->edge() = en0; hn0->face() = fn0;
    hn1->next() = hn5; hn1->twin() = hn0; hn1->vertex() = v1; hn1->edge() = en0; hn1->face() = fn1;
    hn2->next() =  h4; hn2->twin() = hn3; hn2->vertex() = vn; hn2->edge() = en1; hn2->face() =  f0;
    hn3->next() = hn0; hn3->twin() = hn2; hn3->vertex() = v2; hn3->edge() = en1; hn3->face() = fn0;
    hn4->next() =  h1; hn4->twin() = hn5; hn4->vertex() = v3; hn4->edge() = en2; hn4->face() =  f1;
    hn5->next() =  h5; hn5->twin() = hn4; hn5->vertex() = vn; hn5->edge() = en2; hn5->face() = fn1;

    v1 ->halfedge() =  h2; vn ->halfedge() = hn0;
    en0->halfedge() = hn0; en1->halfedge() = hn2; en2->halfedge() = hn4;
    f0 ->halfedge() =  h0; f1 ->halfedge() =  h1; fn0->halfedge() = hn0; fn1->halfedge() = hn1;
    // Mark New edges
    if (markNew){
      /*e0 ->isNew = true; en0->isNew = true; */
      en1->isNew = true; en2->isNew = true;
    }
    return vn;
  }

  VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {

    // This method should collapse the given edge and return an iterator to
    // the new vertex created by the collapse.


    // PHASE I: Collect Elements
    // HALFEDGES
    HalfedgeIter h0 = e->halfedge();
    HalfedgeIter h1 = h0->twin();
    HalfedgeIter h2 = h0->next();
    HalfedgeIter h3 = h1->next();
    HalfedgeIter h4 = h0; do{h4 = h4->next();} while(h4->next() != h0);
    HalfedgeIter h5 = h1; do{h5 = h5->next();} while(h5->next() != h1);

    // Corner Case (faces with two identical edges continuously):
    // doesn't do anything ^_^
    if (h2->twin() == h5 || h3->twin() == h4) return h0->vertex();

    HalfedgeIter h6 = h2->twin();
    HalfedgeIter h7 = h3->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h5->twin();

    // VERTICES
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h1->vertex();

    // EDGES
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h3->edge();
    EdgeIter e4 = h4->edge();
    EdgeIter e5 = h5->edge();

    // FACES
    FaceIter f0 = h0->face();
    FaceIter f1 = h1->face();

    //PHASE II: Allocate new elements
    VertexIter vn = newVertex();
    vn->position = (v0->position + v1->position) / 2.f;
    vn->halfedge() = h9;

    // PHASE III: Reassign Elements
    // HALFEDGES
    h4->next() = h2;
    h5->next() = h3;
    HalfedgeIter h = h3;
    do{
      h->vertex() = vn;
      h = h->twin()->next();
    } while(h != h3);

    // VERTICES
    v1->halfedge() = h2->twin()->next();

    // FACES done in PHASE IV

    //PHASE IV: Delete unused elements
    //delete faces with less than 3 edges
    if (h2->next()->next() == h2){
      h6->twin() = h8; h8->twin() = h6; h8->edge() = e2;
      e2->halfedge() = h6;
      h4->vertex()->halfedge() = h6;
      deleteHalfedge(h2); deleteHalfedge(h4);
      deleteEdge(e4);
      deleteFace(f0);
    }
    else{f0->halfedge() = h2;}
    if (h3->next()->next() == h3){
      h7->twin() = h9; h9->twin() = h7; h7->edge() = e5;
      e5->halfedge() = h9;
      h5->vertex()->halfedge() = h7;
      deleteHalfedge(h3); deleteHalfedge(h5);
      deleteEdge(e3);
      deleteFace(f1);
    }
    else{f1->halfedge() = h3;}

    //delete common-case elements
    deleteEdge(e);
    deleteVertex(v0); deleteVertex(v1);
    deleteHalfedge(h0); deleteHalfedge(h1);

    return vn;
  }

  VertexIter HalfedgeMesh::collapseFace(FaceIter f) {

    // TODO: (meshEdit)
    // This method should collapse the given face and return an iterator to
    // the new vertex created by the collapse.
    return VertexIter();
  }

  FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {

    // TODO: (meshEdit)
    // This method should replace the given vertex and all its neighboring
    // edges and faces with a single face, returning the new face.
    return FaceIter();
  }

  FaceIter HalfedgeMesh::eraseEdge( EdgeIter e ) {
    // This method should erase the given edge and return an iterator to the
    // merged face.

    // Don't erase edge on boundary
    if (e->halfedge()->isBoundary() || e->halfedge()->twin()->isBoundary()){
      return e->halfedge()->face();
    }
    // PHASE I: Collect Elements
    // HALFEDGES

    HalfedgeIter h0 = e->halfedge();
    HalfedgeIter h1 = h0->twin();
    HalfedgeIter h2 = h0->next();
    HalfedgeIter h3 = h1->next();
    HalfedgeIter h4 = h0; do{h4 = h4->next();} while(h4->next() != h0);
    HalfedgeIter h5 = h1; do{h5 = h5->next();} while(h5->next() != h1);

    // VERTICES
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h1->vertex();

    // EDGES
    // FACES
    FaceIter f0 = h0->face();
    FaceIter f1 = h1->face();

    //PHASE II: Allocate new elements
    FaceIter fn = newFace();

    // PHASE III: Reassign Elements
    // HALFEDGES
    h4->next() = h3;
    h5->next() = h2;

    HalfedgeIter h = h2;
    do{
      h->face() = fn;
      h = h->next();
    } while(h != h2);

    //VERTICES
    v0->halfedge() = h3;
    v1->halfedge() = h2;

    //EDGES
    //FACES
    fn->halfedge() = h2;

    //PHASE IV: Delete unused elements
    deleteEdge(e);
    deleteHalfedge(h0);
    deleteHalfedge(h1);
    deleteFace(f0);
    deleteFace(f1);

    return fn;
  }

  EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
    // This method should flip the given edge and return an iterator to the
    // flipped edge.

    // PHASE I: Collect Elements
    // HALFEDGES
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h0; do{h2 = h2->next();} while(h2->next() != h0);
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h3; do{h5 = h5->next();} while(h5->next() != h3);
    HalfedgeIter h6 = h1->next();
    HalfedgeIter h7 = h4->next();

    // VERTICES
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h1->twin()->vertex();
    VertexIter v3 = h4->twin()->vertex();
    // EDGES
    // FACES
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    // PHASE III: Reassign Elements
    // HALFEDGES
    h0->next() = h6; h0->vertex() = v3;
    h1->next() = h3; h1->face() = f1;
    h2->next() = h4;
    h3->next() = h7; h3->vertex() = v2;
    h4->next() = h0; h4->face() = f0;
    h5->next() = h1;

    // VERTICES
    v0->halfedge() = h4;
    v1->halfedge() = h1;

    // EDGES
    // FACES
    f0->halfedge() = h0;
    f1->halfedge() = h3;

    return e0;
  }

  void HalfedgeMesh::subdivideQuad( bool useCatmullClark )
  {
    // Unlike the local mesh operations (like bevel or edge flip), we will perform
    // subdivision by splitting *all* faces into quads "simultaneously."  Rather
    // than operating directly on the halfedge data structure (which as you've seen
    // is quite difficult to maintain!) we are going to do something a bit nicer:
    //
    //    1. Create a raw list of vertex positions and faces (rather than a full-
    //       blown halfedge mesh).
    //
    //    2. Build a new halfedge mesh from these lists, replacing the old one.
    //
    // Sometimes rebuilding a data structure from scratch is simpler (and even more
    // efficient) than incrementally modifying the existing one.  These steps are
    // detailed below.

    // Step I: Compute the vertex positions for the subdivided mesh.  Here we're
    // going to do something a little bit strange: since we will have one vertex in
    // the subdivided mesh for each vertex, edge, and face in the original mesh, we
    // can nicely store the new vertex *positions* as attributes on vertices, edges,
    // and faces of the original mesh.  These positions can then be conveniently
    // copied into the new, subdivided mesh.
    // [See subroutines for actual "TODO"s]
    if( useCatmullClark )
    {
      computeCatmullClarkPositions();
    }
    else
    {
      computeLinearSubdivisionPositions();
    }

    // Step II: Assign a unique index (starting at 0) to each vertex, edge, and
    // face in the original mesh.  These indices will be the indices of the vertices
    // in the new (subdivided mesh).  They do not have to be assigned in any particular
    // order, so long as no index is shared by more than one mesh element, and the
    // total number of indices is equal to V+E+F, i.e., the total number of vertices
    // plus edges plus faces in the original mesh.  Basically we just need a one-to-one
    // mapping between original mesh elements and subdivided mesh vertices.
    assignSubdivisionIndices();

    // Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
    // the element indices defined above.  In other words, each new quad should be of
    // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
    // original mesh elements.  Note that it is essential to get the orientation right
    // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
    // circulate in the same direction as old faces (think about the right-hand rule).
    vector< vector<Index> > subDFaces;
    vector< Vector3D > subDVertices;
    buildSubdivisionFaceList( subDFaces );
    buildSubdivisionVertexList( subDVertices );

    // Step IV: Pass the list of vertices and quads to a routine that clears the
    // internal data for this halfedge mesh, and builds new halfedge data from scratch,
    // using the two lists.
    rebuild( subDFaces, subDVertices );
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * simple linear interpolation, e.g., the edge midpoints and face
   * centroids.
   */
  void HalfedgeMesh::computeLinearSubdivisionPositions()
  {
    // For each vertex, assign Vertex::newPosition to
    // its original position, Vertex::position.
    VertexIter v = verticesBegin();
    while (v != verticesEnd()){
      VertexIter nextVertex = v;
      nextVertex++;
      v->newPosition = Vector3D (v->position.x, v->position.y, v->position.z);
      v = nextVertex;
    }

    // For each edge, assign the midpoint of the two original
    // positions to Edge::newPosition.
    EdgeIter e = edgesBegin();
    VertexCIter v0, v1;
    while (e != edgesEnd()){
      EdgeIter nextEdge = e;
      nextEdge++;
      v0 = e->halfedge()        ->vertex();
      v1 = e->halfedge()->twin()->vertex();
      e->newPosition = (v0->position + v1->position) / 2.;
      e = nextEdge;
    }

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::newPosition.  Note
    // that in general, NOT all faces will be triangles!
    FaceIter f = facesBegin();
    HalfedgeCIter h;
    while (f != facesEnd()){
      FaceIter nextface = f;
      nextface++;
      h = f->halfedge();
      double count = 0.;
      f->newPosition = Vector3D(0.,0.,0.);
      do{
        f->newPosition += h->vertex()->position;
        count += 1.;
        h = h->next();
      } while(h != f->halfedge());
      f->newPosition /= count;
      f = nextface;
    }
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * the Catmull-Clark rules for subdivision.
   */
  void HalfedgeMesh::computeCatmullClarkPositions()
  {
    // The implementation for this routine should be
    // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules.  (These rules are outlined in the Developer Manual.)

    // face
    FaceIter f = facesBegin();
    HalfedgeCIter h;
    while (f != facesEnd()){
      FaceIter nextface = f;
      nextface++;
      h = f->halfedge();
      double count = 0.;
      f->newPosition = Vector3D(0.,0.,0.);
      do{
        f->newPosition += h->vertex()->position;
        count += 1.;
        h = h->next();
      } while(h != f->halfedge());
      f->newPosition /= count;
      f = nextface;
    }
    // edges
    EdgeIter e = edgesBegin();
    VertexCIter v0, v1;
    FaceCIter f0, f1;
    while (e != edgesEnd()){
      EdgeIter nextEdge = e;
      nextEdge++;
      v0 = e->halfedge()        ->vertex();
      v1 = e->halfedge()->twin()->vertex();
      f0 = e->halfedge()		  ->face();
      f1 = e->halfedge()->twin()->face();
      e->newPosition = (v0->position + v1->position +
                        f0->newPosition + f1->newPosition) / 4.;
      e = nextEdge;
    }
    // vertices
    VertexIter v = verticesBegin();
    while (v != verticesEnd()){
      VertexIter nextVertex = v;
      HalfedgeCIter hv = v->halfedge();
      double n = (double) v->degree();
      v->newPosition = Vector3D(0.,0.,0.);
      do {
        v->newPosition += hv->face()->newPosition;
        v->newPosition += hv->edge()->newPosition * 2.;
        hv = hv->twin()->next();
      } while (hv != v->halfedge());
      v->newPosition /= n;
      v->newPosition += (n - 3) * v->position;
      v->newPosition /= n;
      nextVertex++;
      v = nextVertex;
    }
  }

  /**
   * Assign a unique integer index to each vertex, edge, and face in
   * the mesh, starting at 0 and incrementing by 1 for each element.
   * These indices will be used as the vertex indices for a mesh
   * subdivided using Catmull-Clark (or linear) subdivision.
   */
  void HalfedgeMesh::assignSubdivisionIndices()
  {
    // Start a counter at zero; if you like, you can use the
    // "Index" type (defined in halfedgeMesh.h)
    Index counter = 0;
    // Iterate over vertices, assigning values to Vertex::index
    VertexIter v = verticesBegin();
    while (v != verticesEnd()){
      VertexIter nextVertex = v;
      nextVertex++;
      v->index = counter;
      v = nextVertex;
      counter++;
    }
    // Iterate over edges, assigning values to Edge::index
    EdgeIter e = edgesBegin();
    while (e != edgesEnd()){
      EdgeIter nextEdge = e;
      nextEdge++;
      e->index = counter;
      e = nextEdge;
      counter++;
    }
    // Iterate over faces, assigning values to Face::index
    FaceIter f = facesBegin();
    while (f != facesEnd()){
      FaceIter nextface = f;
      nextface++;
      f->index = counter;
      f = nextface;
      counter++;
    }
  }

  /**
   * Build a flat list containing all the vertex positions for a
   * Catmull-Clark (or linear) subdivison of this mesh.  The order of
   * vertex positions in this list must be identical to the order
   * of indices assigned to Vertex::newPosition, Edge::newPosition,
   * and Face::newPosition.
   */
  void HalfedgeMesh::buildSubdivisionVertexList( vector<Vector3D>& subDVertices )
  {
    // Resize the vertex list so that it can hold all the vertices.
    size_t v_size = vertices.size(), e_size = edges.size(), f_size = faces.size();
    size_t size = v_size + e_size + f_size;
    subDVertices.resize(size);
    Index i = 0;

    // Iterate over vertices, assigning Vertex::newPosition to the appropriate
    // location in the new vertex list.
    VertexCIter v = verticesBegin();
    while (v != verticesEnd()){
      subDVertices[i] = v->newPosition;
      i++; v++;
    }

    // Iterate over edges, assigning Edge::newPosition to the appropriate
    // location in the new vertex list.
    EdgeCIter e = edgesBegin();
    while (e != edgesEnd()){
      subDVertices[i] = e->newPosition;
      i++; e++;
    }

    // Iterate over faces, assigning Face::newPosition to the appropriate
    // location in the new vertex list.
    FaceCIter f = facesBegin();
    while (f != facesEnd()){
      subDVertices[i] = f->newPosition;
      i++; f++;
    }
  }

  /**
   * Build a flat list containing all the quads in a Catmull-Clark
   * (or linear) subdivision of this mesh.  Each quad is specified
   * by a vector of four indices (i,j,k,l), which come from the
   * members Vertex::index, Edge::index, and Face::index.  Note that
   * the ordering of these indices is important because it determines
   * the orientation of the new quads; it is also important to avoid
   * "bowties."  For instance, (l,k,j,i) has the opposite orientation
   * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
   * will look like a bowtie.
   */
  void HalfedgeMesh::buildSubdivisionFaceList( vector< vector<Index> >& subDFaces )
  {
    // This routine is perhaps the most tricky step in the construction of
    // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
    // vertex positions).  Basically what you want to do is iterate over faces,
    // then for each for each face, append N quads to the list (where N is the
    // degree of the face).  For this routine, it may be more convenient to simply
    // append quads to the end of the list (rather than allocating it ahead of
    // time), though YMMV.  You can of course iterate around a face by starting
    // with its first halfedge and following the "next" pointer until you get
    // back to the beginning.  The tricky part is making sure you grab the right
    // indices in the right order---remember that there are indices on vertices,
    // edges, AND faces of the original mesh.  All of these should get used.  Also
    // remember that you must have FOUR indices per face, since you are making a
    // QUAD mesh!

    FaceCIter f = facesBegin();
    HalfedgeCIter h, next_h;
    // iterate over faces
    while(f != facesEnd()){
      h = f->halfedge();
      // loop around face
      do {
        next_h = h->next();
        // build lists of four indices for each sub-quad
        vector<Index> quad(4);
        quad[0] = h->edge()->index;
        quad[1] = next_h->vertex()->index;
        quad[2] = next_h->edge()->index;
        quad[3] = f->index;
        // append each list of four indices to face list
        subDFaces.push_back(quad);
        h = next_h;
      } while (h != f->halfedge());
      f++;
    }
  }

  void HalfedgeMesh::_bevel_fc_reposition_with_dist( vector<Vector3D>& orig, // list of vertex positions of the original face (before bevel)
                                                     vector<HalfedgeIter>& hs, // list of halfedges pointing from the vertices of the new, beveled face to the old, original face
                                                     double shift, // user-requested amount to shift the face in the normal direction
                                                     double inset ) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled face.
    // FIXME Not a good interaction design
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    int N = hs.size();
    FaceCIter f = hs[0]->twin()->next()->twin()->face();
    Vector3D normal = f->normal();
    for( int i = 0; i < N; i++ )
    {
      Vector3D pi = hs[i]->next()->vertex()->position; // get the original vertex position correponding to vertex i
      Vector3D p0 = hs[(i+N-1) % N]->next()->vertex()->position;
      Vector3D p1 = hs[(i+1) % N]->next()->vertex()->position;
      // get inset direction
      Vector3D vi0 = p0 - pi, vi1 = p1 - pi;
      vi0.normalize(); vi1.normalize();
      hs[i]->vertex()->position += inset * (vi0 + vi1) - shift * normal;
    }
  }

  void HalfedgeMesh::_bevel_vtx_reposition_with_dist( Vector3D orig, // original vertex position, before the bevel
                                                      vector<HalfedgeIter>& hs, // list of halfedges pointing from the vertices of the new, beveled face to the neighbors of the original vertex
                                                      double inset ) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled vertex.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
  }

  void HalfedgeMesh::_bevel_edge_reposition_with_dist( vector<Vector3D>& origs,  // list of vertex positions of the neighbors of the two endpoints of the edge, before the bevel
                                                       vector<HalfedgeIter>& hs,  // list of halfedges pointing from the vertices of the new, beveled face to the neighbors of the endpoints of the old, original edge
                                                       double inset) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled edge.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //
  }

  FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {

    // TODO This method should replace the vertex v with a face, corresponding to a bevel operation.
    // It should return the new face.  NOTE: This method is responsible for updating the *connectivity*
    // of the mesh only---it does not need to update the vertex positions.  These positions will be
    // updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to implement!)

    return facesBegin();
  }

  FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {

    // TODO This method should replace the edge e with a face, corresponding to a bevel operation.
    // It should return the new face.  NOTE: This method is responsible for updating the *connectivity*
    // of the mesh only---it does not need to update the vertex positions.  These positions will be
    // updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to implement!)

    return facesBegin();
  }

  FaceIter HalfedgeMesh::bevelFace(FaceIter f) {
    // This method should replace the face f with an additional, inset face (and ring of faces around it),
    // corresponding to a bevel operation. It should return the new face.  NOTE: This method is responsible for
    // updating the *connectivity* of the mesh only---it does not need to update the vertex positions.  These
    // positions will be updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to
    // implement!)

    // Some variables
    HalfedgeIter h = f->halfedge(); HalfedgeIter ht = h;  HalfedgeIter tmp;

    // New inner face
    FaceIter fn = newFace();

    // Make first new side-face
    HalfedgeIter hn0 = newHalfedge(); HalfedgeIter hn1 = newHalfedge();
    HalfedgeIter hn2 = newHalfedge(); HalfedgeIter hn3 = newHalfedge();
    VertexIter vn0 = newVertex();     VertexIter vn1 = newVertex();
    EdgeIter en0 = newEdge(); EdgeIter en1 = newEdge(); EdgeIter en2 = newEdge();
    FaceIter fn0;

    HalfedgeIter hne = hn2; EdgeIter ene = en2; VertexIter vne = vn1; HalfedgeIter hnei = hn3;
    ene->halfedge() = hne; vne->halfedge() = hne;
    hne->edge() = ene; hne->face() = fn0; hne->vertex() = vne; hne->next() = h;

    h = h->next();
    ht->next() = hn0; hn0->next() = hn1; hn1->next() = hn2; hn1->next() = hn2; hn2->next() = ht;
    hn0->vertex() = h->vertex(); hn1->vertex() = vn0; hn2->vertex() = vn1; hn3->vertex() = vn1;
    hn1->twin() = hn3; hn3->twin() = hn1;
    hn0->edge() = en0; hn1->edge() = en1; hn2->edge() = en2; hn3->edge() = en1;
    hn0->face() = f; hn1->face() = f; hn2->face() = f; hn3->face() = fn;

    en0->halfedge() = hn0; en1->halfedge() = hn1; en2->halfedge() = hn2;
    vn0->halfedge() = hn1; vn1->halfedge() = hn2;

    //temporary position to be changed
    vn0->position = hn0->vertex()->position; //vn0->position.z += 1.;

    // Make all new side-faces
    while(h != f->halfedge()){
      // Allocate new stuff
      ht = h; h = h->next();
      hn2 = newHalfedge(); hn2->twin() = hn0; hn0->twin() = hn2;
      hn0 = newHalfedge(); hn1 = newHalfedge(); tmp = hn3; hn3 = newHalfedge(); tmp->next() = hn3;
      vn1 = vn0; en2 = en0;
      // Last iteration
      if(h == f->halfedge()){vn0 = vne; en0 = ene; hn0->twin() = hne; hne->twin() = hn0; hn3->next() = hnei;}
      else{vn0 = newVertex(); en0 = newEdge();}

      en1 = newEdge(); fn0 = newFace();
      // Assign halfedges
      ht->next() = hn0; hn0->next() = hn1; hn1->next() = hn2; hn1->next() = hn2; hn2->next() = ht;
      hn0->vertex() = h->vertex(); hn1->vertex() = vn0; hn2->vertex() = vn1; hn3->vertex() = vn1;
      hn1->twin() = hn3; hn3->twin() = hn1;
      hn0->edge() = en0; hn1->edge() = en1; hn2->edge() = en2; hn3->edge() = en1;
      hn0->face() = fn0; hn1->face() = fn0; hn2->face() = fn0; hn3->face() = fn; ht->face() = fn0;
      // Assign vertices, edges and faces
      en0->halfedge() = hn0; en1->halfedge() = hn1;
      vn0->halfedge() = hn1;
      fn0->halfedge() = hn0;
      // temporary position to be changed
      vn0->position = hn0->vertex()->position; //vn0->position.z += 1.;
    }

    // Assign new inner face
    fn->halfedge() = hn3;

    return fn;
  }

  void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
    for (auto f : fcs) splitPolygon(f);
  }

  void HalfedgeMesh::splitPolygon(FaceIter f) {
    // keep original triangles
    if (f->degree() < 4) return;
    // 1 step triangulation
    HalfedgeIter h0 = f->halfedge(), h1 = h0->next(), h2 = h1->next();
    HalfedgeIter h3 = h0; do{h3 = h3->next();} while(h3->next() != h0);
    HalfedgeIter hn = newHalfedge(), hnt = newHalfedge(); hn->twin() = hnt; hnt->twin() = hn;
    EdgeIter e0 = h0->edge(), e1 = h1->edge(), en = newEdge();
    VertexIter v0 = h0->vertex(), v1 = h1->vertex(), v2 = h2->vertex();
    FaceIter f0 = h0->face(), fn = newFace();

    h1->next() = hn; hn->next() = h0; h3->next() = hnt; hnt->next() = h2;
    hn->edge() = en; hnt->edge() = en;
    hn->vertex() = v2; hnt->vertex() = v0;
    h0->face() = fn; h1->face() = fn; hn->face() = fn; hnt->face() = f0;

    en->halfedge() = hnt; f0->halfedge() = hnt; fn->halfedge() = h0;

    // recursively triangulate
    splitPolygon(f0);
  }

  EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {

    // TODO: (meshEdit)
    // Compute the combined quadric from the edge endpoints.
    // -> Build the 3x3 linear system whose solution minimizes the quadric error
    //    associated with these two endpoints.
    // -> Use this system to solve for the optimal position, and store it in
    //    EdgeRecord::optimalPoint.
    // -> Also store the cost associated with collapsing this edg in
    //    EdgeRecord::Cost.

  }

  void MeshResampler::upsample(HalfedgeMesh& mesh)
  // This routine should increase the number of triangles in the mesh using Loop subdivision.
  {
    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::newPosition.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::newPosition.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::isNew. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::position.

    // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using the Loop subdivision rule.
    VertexIter v = mesh.verticesBegin();
    while (v != mesh.verticesEnd()){
      VertexIter nextVertex = v;
      HalfedgeCIter hv = v->halfedge();
      size_t n = v->degree();
      double u = (n == 3) ? 0.1875 : 0.375/n ;
      v->newPosition = Vector3D(0.,0.,0.);
      do {
        v->newPosition += hv->twin()->vertex()->newPosition;
        hv = hv->twin()->next();
      } while (hv != v->halfedge());
      v->newPosition *= u;
      v->newPosition += (1. - n*u) * v->position;
      v->isNew = false;
      nextVertex++;
      v = nextVertex;
    }
    // Next, compute the updated vertex positions associated with edges.
    EdgeIter e = mesh.edgesBegin();
    VertexCIter v0, v1, v2, v3;
    while (e != mesh.edgesEnd()){
      EdgeIter nextEdge = e;
      nextEdge++;
      v0 = e->halfedge()->vertex();
      v1 = e->halfedge()->twin()->vertex();
      v2 = e->halfedge()->next()->next()->vertex();
      v3 = e->halfedge()->twin()->next()->next()->vertex();
      e->newPosition = (3.*v0->position + 3.*v1->position +
                        v2->position + v3->position) / 8.;
      e->isNew = false;
      e = nextEdge;
    }
    // Next, we're going to split every edge in the mesh, in any order.  For future
    // reference, we're also going to store some information about which subdivided
    // edges come from splitting an edge in the original mesh, and which edges are new.
    // In this loop, we only want to iterate over edges of the original mesh---otherwise,
    // we'll end up splitting edges that we just split (and the loop will never end!)
    vector<EdgeIter> copy_of_old_edges;
    for (e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++){
      copy_of_old_edges.push_back(e);
    }
    size_t edges_size = copy_of_old_edges.size();
    for (Index i = 0; i < edges_size; i++){
      e = copy_of_old_edges[i];
      Vector3D ePos = e->newPosition;
      // split and mark new edges
      v = mesh.splitEdge(e, true);
      v->newPosition = ePos;
      v->isNew = true;
    }
    // Finally, flip any new edge that connects an old and new vertex.
    e = mesh.edgesBegin();
    while (e != mesh.edgesEnd()){
      EdgeIter nextEdge = e;
      nextEdge++;
      // find new edges
      if (e->isNew){
        VertexIter end0 = e->halfedge()->vertex(),
            end1 = e->halfedge()->twin()->vertex();
        // flip any new edge that connects an old and new vertex
        if (end0->isNew ^ end1->isNew){
          mesh.flipEdge(e);
        }
      }
      e = nextEdge;
    }
    // Copy the updated vertex positions to the subdivided mesh.
    v = mesh.verticesBegin();
    while (v != mesh.verticesEnd()){
      VertexIter nextVertex = v;
      nextVertex++;
      v->position = v->newPosition;
      v = nextVertex;
    }
  }

  void MeshResampler::downsample(HalfedgeMesh& mesh)
  {

    // TODO: (meshEdit)
    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in Face::quadric
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in Vertex::quadric
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an EdgeRecord for each edge and sticking it in the
    //    queue.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

  }

  void MeshResampler::resample(HalfedgeMesh &mesh, size_t numProxy, size_t numIterations, double edgeSplitThreshold) {

    // TODO: (meshEdit)
    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // 15418 Final Project
    // We are using this method as an entry point to start VSA routines

    if (!mesh.displayVSA)
      vsaBuilder.build(mesh, numProxy, numIterations, edgeSplitThreshold);
    else
      vsaBuilder.clear();
  }

} // namespace CMU462

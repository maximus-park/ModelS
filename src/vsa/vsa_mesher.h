/**
 * Scotty3D - vsa_mesher.h
 */

#ifndef SCOTTY3D_VSA_MESHER_H
#define SCOTTY3D_VSA_MESHER_H

#include "vsa_proxy.h"
#include "vsa.h"

/** Hash functions */
namespace std {
  template<>
  struct hash<VertexCIter> {
    size_t operator()(VertexCIter const& v) const noexcept {
      size_t const hash_result (std::hash<size_t>{}((size_t) &*v));
      return hash_result;
    }
  };
}

/*
 * Meshing routines for Variational Shape Approximation
 */
class VSAMesher {
public:
  VSAMesher(HalfedgeMeshPtr oldMesh) { this->oldMeshPtr = oldMesh; }
  ~VSAMesher() {}

  void buildMesh(ProxyList &proxyList, double threshold);

private:
  /** Meshing related global variables */
  std::unordered_map<VertexCIter, vector<ProxyLabel>> anchorVertices;
  HalfedgeMeshPtr oldMeshPtr = nullptr; ///< The pointer to the old input mesh

  /**
   * Initialize global anchor vertex set. Attach a border halfedge to all proxies.
   * When this method returns, the global anchorVertices map should contain all vertices
   * that is surrounded by three or more proxies.
   * @param proxyList
   */
  void initAnchorVertices(ProxyList &proxyList);

  /**
   * Split edges and make sure each proxy has at least three anchors
   */
  void refineAnchorVertices(ProxyList &proxyList, double threshold);

  void newAnchor(ProxyList &proxyList, VertexCIter vertex);

  /**
   * Recursively split an edge define by its two end points v1 and v2 in the given proxy.
   * v1 must be the previous anchor point in the given proxy.
   * @param proxyList
   * @param proxy
   * @param v1
   * @param v2
   * @param threshold split threshold, a negative value means don't recurse
   * @return the newly added anchor vertex or
   *         nullptr if split criterion isn't satisfied or split cannot be performed
   */
  VertexCIter splitEdge(ProxyList &proxyList, Proxy &proxy, VertexCIter v1, VertexCIter v2, double threshold);

  void buildNewVertexList(std::vector<Vector3D> &newVertices, std::unordered_map<VertexCIter, Index> &newIndices);
  void buildNewFaceList(std::vector<vector<Index> > &newFaces, ProxyList &proxyList,
                        std::unordered_map<VertexCIter, Index> &newIndices);

}; // class VSAMesher


#endif //SCOTTY3D_VSA_MESHER_H

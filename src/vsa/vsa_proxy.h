/**
 * Scotty3D - vsa_proxy.h
 */

#ifndef SCOTTY3D_VSA_PROXY_H
#define SCOTTY3D_VSA_PROXY_H

#include "../halfEdgeMesh.h"

class Proxy {
public:
  /** The normal of the proxy */
  Vector3D normal;
  /** The centroid of the proxy */
  Vector3D centroid;
  /**
   * The seed of a proxy is the face from which a region can be grown.
   * It is updated every iteration using Lloyd algorithm.
   */
  FaceIter seed;
  /** Total area of this proxy */
  double totalArea;
  /** Global approximation error over the proxy */
  double totalDistorsion;
  /** The label of proxy */
  ProxyLabel label;

  /** Used by meshing routines */
  HalfedgeCIter borderHalfedge;
  size_t borderEdgeCount;
  std::list<VertexCIter> anchors; ///< anchors must be sorted by the same order of halfedge

  /**
   * For a given border halfedge on this proxy, find the next halfedge on border
   * @param he a given border halfedge on this proxy
   * @param invalidReturnValue
   * @return the next halfedge on border or nullptr if input is not on border
   */
  HalfedgeCIter nextHalfedgeOnBoarder(HalfedgeCIter he, HalfedgeCIter invalidReturnValue);

  /**
   * Given a boarder vertex of this proxy, find the outgoing halfedge on border
   * @param v
   * @param invalidReturnValue
   * @return the outgoing halfedge on border
   */
  HalfedgeCIter findHalfedgeOnBorder(VertexCIter v, HalfedgeCIter invalidReturnValue);

  /**
   * Determine whether a given halfedge is on the border of this proxy
   */
  bool isBorder(HalfedgeCIter he) const;

  /**
   * Add an anchor vertex to this proxy. After adding, anchor vector should remain sorted.
   * Requires that this vertex is on the border of this proxy.
   * @param newAnchor
   */
  void addAnchor(VertexCIter newAnchor);
};


#endif //SCOTTY3D_VSA_PROXY_H

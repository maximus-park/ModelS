/**
 * Scotty3D - vsa_mesher.cpp
 * Meshing routines for Variational Shape Approximation
 */

#include "vsa_mesher.h"

void VSAMesher::buildMesh(ProxyList &proxyList, double threshold) {
  initAnchorVertices(proxyList);
  refineAnchorVertices(proxyList, threshold);

  std::unordered_map<VertexCIter, Index> newIndices;
  std::vector< vector<Index> > newFaces;
  std::vector< Vector3D > newVertices;
  buildNewVertexList(newVertices, newIndices);
  buildNewFaceList(newFaces, proxyList, newIndices);
  oldMeshPtr->rebuild( newFaces, newVertices);
}

void VSAMesher::initAnchorVertices(ProxyList &proxyList) {
  // TODO more clever initialization rather than iterating through all vertices
  for (auto v = oldMeshPtr->verticesBegin(); v != oldMeshPtr->verticesEnd(); v++) {
    int proxyCount = 0;
    vector<ProxyLabel> labelMapping;
    auto he = v->halfedge(), h = he;
    do {
      if (h->isBoundary()) {
        proxyCount++;
        continue;
      }
      Proxy &p = proxyList[h->face()->label];
      if (p.isBorder(h)) {
        proxyCount++;
        labelMapping.push_back(p.label);
        if (p.borderHalfedge == oldMeshPtr->halfedgesEnd()) p.borderHalfedge = h;
      }
      h = h->twin()->next(); // next outgoing halfedge of this vertex
    } while (h != he);

    if (proxyCount >= 3) {
      anchorVertices.insert(std::make_pair(v, labelMapping));
    }
  }

  // initialize the sorted anchor list for all proxies
  for (auto &p : proxyList) {
    // visit all border vertices in order
    size_t count = 0;
    auto he = p.borderHalfedge;
    do {
      auto v = he->vertex();
      if (anchorVertices.find(v) != anchorVertices.end()) {
        p.anchors.push_back(v);
      }
      count++;
      he = p.nextHalfedgeOnBoarder(he, oldMeshPtr->halfedgesEnd());
    } while (he != p.borderHalfedge);
    p.borderEdgeCount = count;
  }
}

void VSAMesher::newAnchor(ProxyList &proxyList, VertexCIter vertex) {
  std::unordered_set<ProxyLabel> finishedLabels;
  std::vector<ProxyLabel> labelMapping;
  // find the border halfedges corresponding to the new anchor
  auto he = vertex->halfedge();
  do {
    if (he->isBoundary()) continue;
    ProxyLabel l = he->face()->label;
    if (finishedLabels.find(l) == finishedLabels.end()) {
      proxyList[l].addAnchor(vertex);
      labelMapping.push_back(l);
      finishedLabels.insert(l);
    }
    he = he->twin()->next(); // next outgoing halfedge of this vertex
  } while (he != vertex->halfedge());
  // add to global anchor map
  anchorVertices.insert(std::make_pair(vertex, labelMapping));
  // TODO add to anchor edge list
}

VertexCIter VSAMesher::splitEdge(ProxyList &proxyList, Proxy &proxy,
                                 VertexCIter v1, VertexCIter v2, double threshold) {
  // find the vertex on edge that has largest distance to the vector (v1, v2)
  double largestDistance = 0.0;
  VertexCIter newAnchorVertex;
  Vector3D v1v2 = v2->position - v1->position;
  double edgeLength = v1v2.norm();
  v1v2.normalize();

  auto he = proxy.findHalfedgeOnBorder(v1, oldMeshPtr->halfedgesEnd());
  he = proxy.nextHalfedgeOnBoarder(he, oldMeshPtr->halfedgesEnd());
  while (he->vertex() != v2) {
    Vector3D vec = he->vertex()->position - v1->position;
    double dist = cross(vec, v1v2).norm();
    if (dist > largestDistance) {
      largestDistance = dist;
      newAnchorVertex = he->vertex();
    }
    he = proxy.nextHalfedgeOnBoarder(he, oldMeshPtr->halfedgesEnd());
  }

  if (largestDistance == 0.0) return oldMeshPtr->verticesEnd();

  // non-recursive split:
  if (threshold < 0) {
    newAnchor(proxyList, newAnchorVertex);
    return newAnchorVertex;
  }

  // recursive spilt:
  double sinProxyNormals;
  // calculate split criterion
  auto boarderHalfedge = proxy.findHalfedgeOnBorder(newAnchorVertex, oldMeshPtr->halfedgesEnd());
  if (boarderHalfedge->twin()->isBoundary()) { sinProxyNormals = 1.0; } // we want an accurate edge
  else {
    Vector3D N1 = proxy.normal;
    Vector3D N2 = proxyList[boarderHalfedge->twin()->face()->label].normal;
    sinProxyNormals = cross(N1, N2).norm();
  }
  double splitCriterion = largestDistance * sinProxyNormals / edgeLength;
  if (splitCriterion > threshold) {
    newAnchor(proxyList, newAnchorVertex);
    // recursion
    splitEdge(proxyList, proxy, v1, newAnchorVertex, threshold);
    splitEdge(proxyList, proxy, newAnchorVertex, v2, threshold);
    return newAnchorVertex;
  } else {
    return oldMeshPtr->verticesEnd();
  }
}

void VSAMesher::refineAnchorVertices(ProxyList &proxyList, double threshold) {
  // first make sure every proxy has 3 or more anchors
  for (auto &p : proxyList) {
    if (p.anchors.size() == 0) {
      // let proxy at least have a anchor to start with
      newAnchor(proxyList, p.borderHalfedge->vertex());
    }
    if (p.anchors.size() == 1) {
      // add an anchor to the far side of the original anchor
      int steps = p.borderEdgeCount / 2;
      auto he = p.borderHalfedge;
      for (int i = 0; i < steps; i++) {
        he = p.nextHalfedgeOnBoarder(he, oldMeshPtr->halfedgesEnd());
      }
      newAnchor(proxyList, he->vertex());
    }
    if (p.anchors.size() == 2) {
      // a threshold value of 0.0 means a split must happen regardless of the split criterion
      // but don't recurse split
      splitEdge(proxyList, p, p.anchors.front(), p.anchors.back(), -1.0);
    }
  } // end for loop

  // next recursively split each edge
  // TODO each edge is actually checked twice for split, fix by keeping track of split edges
  for (auto &p : proxyList) {
    auto prev = --p.anchors.end();
    auto next = p.anchors.begin();
    while (next != p.anchors.end()) {
      splitEdge(proxyList, p, *prev, *next, threshold);
      prev = next;
      next++;
    }
  }
}

void VSAMesher::buildNewVertexList(std::vector<Vector3D> &newVertices,
                                   std::unordered_map<VertexCIter, Index> &newIndices) {
  Index i = 0;
  for (auto vpair : anchorVertices) {
    VertexCIter v = vpair.first;
    newIndices.insert(std::make_pair(v, i++));
    // TODO better approximation of vertex position
    newVertices.push_back(v->position);
  }
}

void VSAMesher::buildNewFaceList(std::vector<vector<Index> > &newFaces, ProxyList &proxyList,
                                 std::unordered_map<VertexCIter, Index> &newIndices) {
  for (auto &p : proxyList) {
    std::vector<Index> faceList;
    for (auto v : p.anchors) {
      faceList.push_back(newIndices.find(v)->second);
    }
    newFaces.push_back(faceList);
  }
}
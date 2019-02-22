/**
 * Scotty3D - vsa_proxy.cpp
 */

#include "vsa_proxy.h"

HalfedgeCIter Proxy::nextHalfedgeOnBoarder(HalfedgeCIter he, HalfedgeCIter invalidReturnValue) {
  if (!isBorder(he))
    return invalidReturnValue;
  auto nextVertex = he->next()->vertex();
  return findHalfedgeOnBorder(nextVertex, invalidReturnValue);
}

HalfedgeCIter Proxy::findHalfedgeOnBorder(VertexCIter v, HalfedgeCIter invalidReturnValue) {
  auto he = v->halfedge();
  bool found = false;
  do {
    if (isBorder(he)) {
      found = true;
      break;
    }
    he = he->twin()->next(); // next outgoing halfedge of this vertex
  } while (he != v->halfedge());
  if (!found) return invalidReturnValue;
  else return he;
}

bool Proxy::isBorder(HalfedgeCIter he) const {
  return !he->isBoundary() &&
         he->face()->label == this->label &&
         (he->twin()->isBoundary() || he->twin()->face()->label != this->label);
}

void Proxy::addAnchor(VertexCIter newAnchor) {
  // find the border halfedge corresponding to the new anchor
  auto he = findHalfedgeOnBorder(newAnchor, CMU462::HalfedgeCIter());
  // don't check return value, this shouldn't be nullptr
  // if (he == nullptr) return;

  // TODO duplicate anchor check, we may not need this
  for (auto it = anchors.begin(); it != anchors.end(); it++) {
    if (he->vertex() == *it) {
      return;
    }
  }
  // simple case
  if (anchors.size() < 2) {
    anchors.push_back(newAnchor);
    return;
  }
  // travel along the border to find the anchor after it
  auto h = he;
  do {
    for (auto it = anchors.begin(); it != anchors.end(); it++) {
      if (h->vertex() == *it) {
        anchors.insert(it, newAnchor);
        return;
      }
    }
    h = nextHalfedgeOnBoarder(h, CMU462::HalfedgeCIter());
  } while (h != he);
}
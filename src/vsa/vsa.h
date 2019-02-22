/**
 * Scotty3D - vsa.h
 * Implementation of Variational Shape Approximation
 *
 * Reference Implementation from Guillaume Lavoué (MEPP)
 * Credit: Variational Shape Approximation
 *	       David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004.
 */


#ifndef SCOTTY3D_VSA_H
#define SCOTTY3D_VSA_H

#include <vector>
#include <utility>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "../halfEdgeMesh.h"
#include "CMU462/timer.h"
#include "vsa_proxy.h"
#include "vsa_types.h"

#ifdef CUDA_VSA
#include "vsa_cuda.h"
#endif //CUDA_VSA

typedef HalfedgeMesh* HalfedgeMeshPtr;
typedef std::vector<Proxy> ProxyList;

typedef struct {
  FaceIter face; ///< One face of mesh
  double distance; ///< Distance used for Lloyd algorithm
  ProxyLabel possibleLabel; ///< Possible label of cluster
} MetricFace;

class VSABuilder {
public:

  VSABuilder(){ currentIteration = 0; }
  ~VSABuilder() {}

  /** Segment input mesh with VSA */
  void build(HalfedgeMesh &mesh, size_t numProxies, size_t numIterations, double edgeSplitThreshold);

  /** Clear segmentation data */
  void clear();

  std::vector<Proxy> proxyList; ///< The list of proxies

private:
  size_t numProxies; ///< The number of proxies
  size_t currentIteration; ///< Debug use, VSA iteration by iteration
  double edgeSplitThreshold;
  HalfedgeMeshPtr meshPtr = nullptr; ///< The pointer to the input mesh
  Timer timer;                   ///< performance test timer

  /** Initialize Proxy Randomly */
  void initProxy();
  /** Distortion minimizing flooding */
  void flood();
  /** Proxy Fitting using L_2,1 Error Metric */
  void fitProxy();

  /** Mesh related operation should be abstracted out for portability*/
  double calculateDistortionError(FaceCIter face, Proxy proxy);
  double faceArea(FaceCIter face);

#ifdef CUDA_VSA
  /** Host data structures for CUDA VSA */
  std::vector<FaceCu> hostFaceCu;
  std::vector<ProxyCu> hostProxyCu;
  /** IO functions for CUDA VSA */
  void getCudaResult();
  void buildFaceCu();
  void setFaceNeighbors(FaceCIter f, std::vector<FaceCu> &faceCu, std::map<FaceCIter, FaceIndex> &faceIndexMap);
  void writeFaceNeighbors(std::vector<FaceCu> &faceCu, char* fileName, int numberOfFaces);
  void loadFaceNeighbors(std::vector<FaceCu> &faceCu, char* fileName);
#endif //CUDA_VSA

}; // class VSABuilder



#endif //SCOTTY3D_VSA_H

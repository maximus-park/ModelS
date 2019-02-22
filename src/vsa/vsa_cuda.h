/**
 * Scotty3D - vsa_cuda.h
 * Header file for CUDA based vsa
 */

#ifndef SCOTTY3D_VSA_CUDA_H
#define SCOTTY3D_VSA_CUDA_H

#include <stdio.h>
#include <vector>
#include <string>
#include <set>
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include <device_launch_parameters.h>
#include "cycleTimer.h"
#include "vsa_types.h"

typedef struct {
  FaceIndex faceIndex; ///< One face of mesh
  double distance; ///< Distance used for Lloyd algorithm
  ProxyLabel possibleLabel; ///< Possible label of cluster
} MetricFaceCu;

typedef struct face_cu_s {
  double3 centroid;         ///< this is essentially a Vector3D
  double3 normal;           ///< this is essentially a Vector3D
  FaceIndex neighbors[3];   ///< if a face is on the boundary, its neighbor index should be -1
  double area;              ///< this value should be pre-calculated
  ProxyLabel label;         ///< defaults to -1
  bool isBoundary;          ///< dummy face that represents an open boundary
} FaceCu;

typedef struct proxy_cu_s {
  double3 centroid;         ///< this is essentially a Vector3D
  double3 normal;           ///< this is essentially a Vector3D
  double totalArea;
  FaceIndex seed;
} ProxyCu;

typedef struct config_s {
  unsigned long numFaces;
  unsigned long numProxies;
} CudaVSAConfig;

class CudaVSAPartitioner {
public:
  CudaVSAPartitioner(std::vector<FaceCu> *hostFaceCu, std::vector<ProxyCu> *hostProxyCu);
  ~CudaVSAPartitioner();

  void partition(size_t numIterations);

private:
  CudaVSAConfig vsaConfig;
  // data structures on host
  std::vector<FaceCu> *hostFaceCu;
  std::vector<ProxyCu> *hostProxyCu;
//  std::vector<FaceIndex> hostFaceResult;
  // data structures on device
  FaceCu *deviceFaceCu;
  ProxyCu *deviceProxyCu;
//  FaceIndex *deviceFaceResult;

  // TODO use shared memory for this
//  long *deviceFaceNearestProxy;
  double *deviceLowestNormalDiff;


  /** Build data structures for CUDA */
  void setup();
  /** Initialize Proxy Randomly */
  void initProxy();
  /** Distortion minimizing flooding */
  void flood();
  /** Proxy Fitting using L_2,1 Error Metric */
  void fitProxy(bool updateSeedFace);
  /** Distortion minimizing flooding with a global priority queue */
  void sequentialFlood();
};

struct MetricFaceComp {
  bool operator()(MetricFaceCu const &mf1,
                  MetricFaceCu const &mf2) const {
    return mf1.distance < mf2.distance;
  }
};

typedef std::multiset<MetricFaceCu, MetricFaceComp>::iterator FaceQueueIter;


#endif //SCOTTY3D_VSA_CUDA_H

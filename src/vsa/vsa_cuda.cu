/**
 * Scotty3D - vsa_cuda.cu
 */
#include "vsa_cuda.h"

#define cudaCheckError(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
  {
    fprintf(stderr, "CUDA Error: %s at %s:%d\n",
            cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}
/****************************** Device Parameters *****************************/
#define THREADS_PER_BLOCK 512
#define NEAREST_PROXY_NUM 8

/*************************** Device Constant Memory ***************************/
__constant__ CudaVSAConfig cudaVsaConfig;

/********************************* Device Code ********************************/
__device__ __host__ __inline__ double distance2(double3 v1, double3 v2) {
  double dx = v1.x - v2.x, dy = v1.y - v2.y, dz = v1.z - v2.z;
  return dx*dx + dy*dy + dz*dz;
}
__device__ __host__ __inline__ double3 vecAdd(double3 v1, double3 v2) {
  return make_double3(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}

__device__ __host__ __inline__ double3 vecMul(double3 v1, double scaler) {
  return make_double3(v1.x*scaler, v1.y*scaler, v1.z*scaler);
}

__device__ double atomicMinDouble(double *addr, double val) {
  double old = *addr, assumed;
  if(old <= val) return old;
  do {
    assumed = old;
    old = atomicCAS((unsigned long long int*)addr, __double_as_longlong(assumed), __double_as_longlong(val));
  } while(old != assumed);
  return old;
}

#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do { assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif // #if __CUDA_ARCH__ < 600
#endif // #ifdef __CUDA_ARCH__

__global__ void kernelInitProxy(FaceCu *faces, ProxyCu *proxies) {
  // TODO: k-means++ initialization or some other smart initialization

  // a naive fixed initialization
  size_t offset = cudaVsaConfig.numFaces / cudaVsaConfig.numProxies;
  FaceIndex index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= cudaVsaConfig.numFaces) return;
  if (index % offset == 0) {
    ProxyLabel proxyLabel = index / offset;
    faces[index].label = proxyLabel;
    proxies[proxyLabel].seed      = index;
    proxies[proxyLabel].normal    = faces[index].normal;
    proxies[proxyLabel].centroid  = faces[index].centroid;
  }
}

__global__ void kernelFlood(FaceCu *faces, ProxyCu *proxies) {

  FaceIndex index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= cudaVsaConfig.numFaces) return;
  double nearestProxyDist[NEAREST_PROXY_NUM];
  ProxyLabel nearestProxyLabel[NEAREST_PROXY_NUM];
  double maxDist = 0;
  size_t maxDistIdx = 0;

  // TODO use shared memeory for proxy centroids and normals
  // Find the nearest <NEAREST_PROXY_NUM> proxies to this face
  for (ProxyLabel l = 0; l < cudaVsaConfig.numProxies; l++) {
    double dist2 = distance2(faces[index].centroid, proxies[l].centroid);
    if (l < NEAREST_PROXY_NUM) {
      nearestProxyDist[l] = dist2;
      nearestProxyLabel[l] = l;
      if (dist2 > maxDist) {
        maxDist = dist2;
        maxDistIdx = (size_t)l;
      }
    } else { // l >= NEAREST_PROXY_NUM
      // TODO consider reduce conditional branches
      // replace proxy with largest distance with the current proxy
      if (dist2 < maxDist) {
        nearestProxyDist[maxDistIdx] = dist2;
        nearestProxyLabel[maxDistIdx] = l;
        maxDist = 0;
        for (size_t i = 0; i < NEAREST_PROXY_NUM; i++) {
          if (nearestProxyDist[i] > maxDist) {
            maxDist = nearestProxyDist[i];
            maxDistIdx = i;
          }
        }
      }
    }
  }

  // Among those nearest proxies find one that has the minimum distortion
  // error and assign its label to this host
  ProxyLabel labelToAssign = -1;
  double minError = 4.0; // this is the maximum distance2 of two normalized vector
  for (size_t i = 0; i < NEAREST_PROXY_NUM; i++) {
    double error = distance2(faces[index].normal,
                             proxies[nearestProxyLabel[i]].normal);
    if (error < minError) {
      minError = error;
      labelToAssign = nearestProxyLabel[i];
    }
  }

  faces[index].label = labelToAssign;
}

__global__ void kernelAddUpProxy(FaceCu *faces, ProxyCu *proxies) {
  FaceIndex index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= cudaVsaConfig.numFaces) return;
  ProxyLabel l = faces[index].label;
  double3 areaNormal = vecMul(faces[index].normal, faces[index].area);
  double3 areaCentroid = vecMul(faces[index].centroid, faces[index].area);

  // update total area of proxy
  atomicAdd(&(proxies[l].totalArea), faces[index].area);

  // update normal of proxy
  atomicAdd(&(proxies[l].normal.x), areaNormal.x);
  atomicAdd(&(proxies[l].normal.y), areaNormal.y);
  atomicAdd(&(proxies[l].normal.z), areaNormal.z);

  // update centroid of proxy
  atomicAdd(&(proxies[l].centroid.x), areaCentroid.x);
  atomicAdd(&(proxies[l].centroid.y), areaCentroid.y);
  atomicAdd(&(proxies[l].centroid.z), areaCentroid.z);
}

__global__ void kernelUpdateProxy(ProxyCu *proxies) {
  ProxyLabel label = blockIdx.x * blockDim.x + threadIdx.x;
  if (label >= cudaVsaConfig.numProxies) return;
  proxies[label].centroid = vecMul(proxies[label].centroid, 1.0/proxies[label].totalArea);
  proxies[label].normal = vecMul(proxies[label].normal, 1.0/proxies[label].totalArea);
}

__global__ void kernelUpdateSeedFace(FaceCu *faces, ProxyCu *proxies, double *lowestNormalDiff) {
  FaceIndex index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < cudaVsaConfig.numProxies) {
    lowestNormalDiff[index] = 4.0;
  }
  __syncthreads();

  if (index >= cudaVsaConfig.numFaces || faces[index].isBoundary) return;
  // update seed face for final round of flooding
  ProxyLabel label = faces[index].label;
  double diff = distance2(faces[index].normal, proxies[label].normal);

  // TODO Fix race condition here
  double oldDiff = atomicMinDouble(&lowestNormalDiff[label], diff);
  if (oldDiff != diff) {
    proxies[label].seed = index;
  }
  // clear proxy label of a face
  faces[index].label = -1;
}

/********************************** Host Code *********************************/
CudaVSAPartitioner::CudaVSAPartitioner(std::vector<FaceCu> *hostFaceCu,
                                       std::vector<ProxyCu> *hostProxyCu) {
  this->hostFaceCu = hostFaceCu;
  this->hostProxyCu = hostProxyCu;
  this->vsaConfig.numFaces = hostFaceCu->size();
  this->vsaConfig.numProxies = hostProxyCu->size();
  deviceFaceCu = NULL;
  deviceProxyCu = NULL;
  setup();
}

CudaVSAPartitioner::~CudaVSAPartitioner() {
  if(deviceProxyCu) {
    cudaFree(deviceProxyCu);
    cudaFree(deviceFaceCu);
    cudaFree(deviceLowestNormalDiff);
    cudaCheckError(cudaThreadSynchronize());
//    cudaFree(deviceFaceNearestProxy);
  }
}

void CudaVSAPartitioner::setup() {

  int deviceCount = 0;
  std::string name;
  cudaError_t err = cudaGetDeviceCount(&deviceCount);

  printf("---------------------------------------------------------\n");
  printf("Initializing CUDA for VSA Partitioner\n");
  printf("Found %d CUDA devices\n", deviceCount);

  for (int i = 0; i < deviceCount; i++) {
    cudaDeviceProp deviceProps;
    cudaGetDeviceProperties(&deviceProps, i);
    name = deviceProps.name;

    printf("Device %d: %s\n", i, deviceProps.name);
    printf("   SMs:        %d\n", deviceProps.multiProcessorCount);
    printf("   Global mem: %.0f MB\n", static_cast<float>(deviceProps.totalGlobalMem) / (1024 * 1024));
    printf("   CUDA Cap:   %d.%d\n", deviceProps.major, deviceProps.minor);
  }
  printf("---------------------------------------------------------\n");

  // allocate device memory
  double startTime, endTime;
  startTime = CycleTimer::currentSeconds();
  cudaCheckError(cudaMalloc(&deviceFaceCu, sizeof(FaceCu) * vsaConfig.numFaces));
  cudaCheckError(cudaMalloc(&deviceProxyCu, sizeof(ProxyCu) * vsaConfig.numProxies));
  cudaCheckError(cudaMalloc(&deviceLowestNormalDiff, sizeof(double) * vsaConfig.numProxies));
  cudaCheckError(cudaThreadSynchronize());
  endTime = CycleTimer::currentSeconds();
  fprintf(stdout, "[VSA CUDA] cudaMalloc()        (%.4f sec)\n", endTime - startTime);
//  cudaCheckError(cudaMalloc(&deviceFaceNearestProxy, sizeof(ProxyLabel) * vsaConfig.numFaces));
}

void CudaVSAPartitioner::partition(size_t numIterations) {
  // TODO allocate result proxy mapping data structure which is essentially map of (face -> proxy label)
  double startTime, endTime;
  // Copy data from host to device
  startTime = CycleTimer::currentSeconds();
  cudaCheckError(cudaMemcpy(deviceFaceCu, &(*hostFaceCu)[0],
                            sizeof(FaceCu) * hostFaceCu->size(), cudaMemcpyHostToDevice));
  cudaCheckError(cudaMemcpy(deviceProxyCu, &(*hostProxyCu)[0],
                            sizeof(ProxyCu) * hostProxyCu->size(), cudaMemcpyHostToDevice));
  cudaCheckError(cudaMemcpyToSymbol(cudaVsaConfig, &vsaConfig, sizeof(CudaVSAConfig)));
  cudaCheckError(cudaThreadSynchronize());
  endTime = CycleTimer::currentSeconds();
  fprintf(stdout, "[VSA CUDA] Host to Device I/O  (%.4f sec)\n", endTime - startTime);

  // Actual Compute
  startTime = CycleTimer::currentSeconds();
  initProxy();
  endTime = CycleTimer::currentSeconds();
  fprintf(stdout, "[VSA CUDA] Initialization Time (%.4f sec)\n", endTime - startTime);
  double flooding_time = 0.0, fitting_time = 0.0, sequential_flooding_time = 0.0;

  startTime = CycleTimer::currentSeconds();
  flood();
  endTime = CycleTimer::currentSeconds();
  flooding_time += endTime - startTime;

  // Lloyd Iterations
  for (size_t i = 1; i < numIterations; i++) {
    startTime = CycleTimer::currentSeconds();
    fitProxy(false);
    endTime = CycleTimer::currentSeconds();
    fitting_time += endTime - startTime;

    startTime = CycleTimer::currentSeconds();
    flood();
    endTime = CycleTimer::currentSeconds();
    flooding_time += endTime - startTime;
  }

  // fit proxy with seed face update
  startTime = CycleTimer::currentSeconds();
  fitProxy(true);
  endTime = CycleTimer::currentSeconds();
  fitting_time += endTime - startTime;


  fprintf(stdout, "[VSA CUDA] Flooding Time       (%.4f sec)\n", flooding_time);
  fprintf(stdout, "[VSA CUDA] Proxy fitting Time  (%.4f sec)\n", fitting_time);

  // Copy data back to host
  startTime = CycleTimer::currentSeconds();
  cudaCheckError(cudaMemcpy(&(*hostFaceCu)[0], deviceFaceCu,
                            sizeof(FaceCu) * vsaConfig.numFaces, cudaMemcpyDeviceToHost));
  cudaCheckError(cudaMemcpy(&(*hostProxyCu)[0], deviceProxyCu,
                            sizeof(ProxyCu) * vsaConfig.numProxies, cudaMemcpyDeviceToHost));
  cudaCheckError(cudaThreadSynchronize());
  endTime = CycleTimer::currentSeconds();
  fprintf(stdout, "[VSA CUDA] Device to Host I/O  (%.4f sec)\n", endTime - startTime);

  startTime = CycleTimer::currentSeconds();
  sequentialFlood();
  endTime = CycleTimer::currentSeconds();
  sequential_flooding_time += endTime - startTime;
  fprintf(stdout, "[VSA CUDA] Seq Flooding Time   (%.4f sec)\n", sequential_flooding_time);
}

/** Initialize Proxy Randomly */
void CudaVSAPartitioner::initProxy() {
  const size_t blocks = (vsaConfig.numFaces + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  kernelInitProxy<<<blocks, THREADS_PER_BLOCK>>>(deviceFaceCu, deviceProxyCu);
  cudaCheckError(cudaThreadSynchronize());
}
/**
 * A data parallel approach of distortion minimizing flooding:
 *
 * 1. We first compute the euclidean distance (as approximation to geodesic distance)
 *    of each face to each proxy centroid.
 * 2. Pick the nearest (maybe 4) proxies of each face and compute the differences of
 *    their normals. Assign proxy label of least normal diff to face.
 * 3. (Maybe we can omit this step) Clean up disconnected clusters with same label
 *    using BFS, keep the cluster with largest total area, unlabel the rest of faces
 * 4. Re-flood using the sequential method to ensure connectivity
 *
 * Inspired by: Fan, Fengtao, et al. "Mesh clustering by approximating centroidal Voronoi tessellation."
 */
void CudaVSAPartitioner::flood() {
  const size_t blocks = (vsaConfig.numFaces + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  kernelFlood<<<blocks, THREADS_PER_BLOCK>>>(deviceFaceCu, deviceProxyCu);

  cudaCheckError(cudaThreadSynchronize());
}
/** Proxy Fitting using L_2,1 Error Metric */
void CudaVSAPartitioner::fitProxy(bool updateSeedFace) {
  // clear old proxy values
  cudaMemset(deviceProxyCu, 0, sizeof(ProxyCu) * vsaConfig.numProxies);
  cudaCheckError(cudaThreadSynchronize());
  // add up related values per proxy
  const size_t faceBlocks = (vsaConfig.numFaces + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  kernelAddUpProxy<<<faceBlocks, THREADS_PER_BLOCK>>>(deviceFaceCu, deviceProxyCu);
  cudaCheckError(cudaThreadSynchronize());
  // update normal and seed of all proxies
  const size_t proxyBlocks = (vsaConfig.numProxies + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  kernelUpdateProxy<<<proxyBlocks, THREADS_PER_BLOCK>>>(deviceProxyCu);
  cudaCheckError(cudaThreadSynchronize());
  // We pick the face with the lowest difference of normal as the new seed
  if (updateSeedFace) {
    kernelUpdateSeedFace <<<faceBlocks, THREADS_PER_BLOCK>>>(deviceFaceCu, deviceProxyCu, deviceLowestNormalDiff);
    cudaCheckError(cudaThreadSynchronize());
  }
}

void CudaVSAPartitioner::sequentialFlood() {
  // Global priority queue of faces
  std::multiset<MetricFaceCu, MetricFaceComp> faceQueue;

  // For each proxy, first add its seed face to priority queue
  for (ProxyLabel label = 0; label < vsaConfig.numProxies; label++) {
    ProxyCu proxy = (*hostProxyCu)[label];
    FaceIndex seedFace = proxy.seed;

    MetricFaceCu m;
    m.faceIndex = seedFace;
    m.distance = 0.0;
    m.possibleLabel = label;

    faceQueue.insert(m);
  }

  // For each face in the priority queue, perform flooding by adding its adjacent faces into the priority queue
  for (FaceQueueIter queueIter = faceQueue.begin(); queueIter != faceQueue.end(); ) {
    FaceIndex face = queueIter->faceIndex;
    FaceCu *facePtr = &(*hostFaceCu)[face];
    if (facePtr->label == -1) { // this means the face hasn't been labeled before
      // label this face, because it is supposed to have the least distortion error
      facePtr->label = queueIter->possibleLabel;
      // add all unlabeled adjacent faces to the priority queue
      for (int i = 0; i < 3; i++) {
        FaceIndex neighborIndex = facePtr->neighbors[i];
        if (neighborIndex < 0) continue;
        FaceCu *neighborPtr = &(*hostFaceCu)[neighborIndex];
        if (!neighborPtr->isBoundary && neighborPtr->label == -1) {
          MetricFaceCu m;
          m.faceIndex = neighborIndex;
          m.distance = neighborPtr->area *
                       distance2(neighborPtr->normal, (*hostProxyCu)[queueIter->possibleLabel].normal);
          m.possibleLabel = queueIter->possibleLabel;
          faceQueue.insert(m);
        }
      }
    }
    // pop the top of priority queue
    faceQueue.erase(queueIter);
    queueIter = faceQueue.begin();
  }
}
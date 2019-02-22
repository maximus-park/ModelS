/**
 * Scotty3D - vsa.cpp
 * Implementation of Variational Shape Approximation
 *
 * Reference Implementation from Guillaume Lavoué (MEPP)
 * Credit: Variational Shape Approximation
 *	       David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004.
 */

#include "vsa.h"
#include "vsa_mesher.h"

#include <iostream>
#include <fstream>
#include <string>

void VSABuilder::build(HalfedgeMesh &mesh, size_t numProxies, size_t numIterations, double edgeSplitThreshold) {
  this->meshPtr = &mesh;
  this->numProxies = numProxies;
  this->edgeSplitThreshold = edgeSplitThreshold;
  meshPtr->displayVSA = true;

  // Debug routines to check convergence
//  if (currentIteration == 0) {
//    init();
//    flood();
//    currentIteration++;
//  } else if (currentIteration < numIterations) {
//    fitProxy();
//    flood();
//    currentIteration++;
//  } else {
//    clear();
//  }
#ifdef CUDA_VSA
  /** GPU segmentation routine - parallel */
  fprintf(stdout, "[VSA CUDA] Input Face Count = %ld, Proxy Count = %ld, Iterations = %ld\n",
          meshPtr->nFaces(), numProxies, numIterations);
  // first gather mesh faces to a random access adjacency list data structure
  timer.start();
  buildFaceCu();
  timer.stop();
  fprintf(stdout, "[VSA CUDA] buildFaceCu()       (%.4f sec)\n", timer.duration());
  // allocate result proxy list data structure of type ProxyList
  proxyList.resize(numProxies);
  // allocate result proxy mapping data structure which is essentially map of (face -> proxy label)
  hostProxyCu.resize(numProxies);
  // run CUDA VSA segmentation routine, wait for it to finish
  CudaVSAPartitioner p = CudaVSAPartitioner(&hostFaceCu, &hostProxyCu);
  p.partition(numIterations);
  // copy back result data structure and assign labels to faces in the original halfedge mesh
  timer.start();
  getCudaResult();
  timer.stop();
  fprintf(stdout, "[VSA CUDA] getCudaResult()     (%.4f sec)\n", timer.duration());
  return;
#endif

  /** CPU segmentation routine - sequential */
  fprintf(stdout, "[VSA] Input Face Count = %ld, Proxy Count = %ld, Iterations = %ld\n",
          this->meshPtr->nFaces(), numProxies, numIterations);
  timer.start();
  initProxy();
  timer.stop();
  fprintf(stdout, "[VSA] Initialization Time (%.4f sec)\n", timer.duration());
  double flooding_time = 0.0, fitting_time = 0.0;

  timer.start();
  flood();
  timer.stop();
  flooding_time += timer.duration();

  // Lloyd Iterations
  for (size_t i = 0; i < numIterations; i++) {
    timer.start();
    fitProxy();
    timer.stop();
    fitting_time += timer.duration();

    timer.start();
    flood();
    timer.stop();
    flooding_time += timer.duration();
  }
  fprintf(stdout, "[VSA] Flooding Time       (%.4f sec)\n", flooding_time);
  fprintf(stdout, "[VSA] Proxy fitting Time  (%.4f sec)\n", fitting_time);
}

void VSABuilder::clear() {
  if (meshPtr == nullptr) return;
  timer.start();
  VSAMesher mesher = VSAMesher(meshPtr);
  mesher.buildMesh(proxyList, edgeSplitThreshold);
  timer.stop();
  fprintf(stdout, "[VSA] Meshing Time        (%.4f sec)\n", timer.duration());

  for (auto f = meshPtr->facesBegin(); f != meshPtr->facesEnd(); f++) {
    f->label = -1;
  }
  proxyList.clear();
  currentIteration = 0;
  meshPtr->displayVSA = false;
}

/** Initialize Proxy Randomly */
void VSABuilder::initProxy() {
  size_t numFaces = meshPtr->nFaces();
  size_t offset = numFaces / numProxies;

  // TODO: Optimization with random access array, maybe with KMeans++
  // Iterate the whole set of faces and initialize proxy at equal distance
  // mimicking a procedure of initializing at random
  size_t counter = 0;
  ProxyLabel label = 0;
  for (auto f = meshPtr->facesBegin(); f != meshPtr->facesEnd(); f++) {
    if (counter % offset == 0) {
      Proxy newProxy;
      newProxy.normal = f->normal();
      newProxy.centroid = f->centroid();
      newProxy.seed = f;
      newProxy.totalArea = 0.0;
      newProxy.label = label;
      newProxy.borderHalfedge = meshPtr->halfedgesEnd();
      proxyList.push_back(newProxy);
      if (++label >= numProxies) break;
    }
    counter++;
  }
}

/** Distortion minimizing flooding */
void VSABuilder::flood() {

  // Global priority queue of faces
  auto metricComp = [](MetricFace mf1, MetricFace mf2) { return mf1.distance < mf2.distance; };
  auto faceQueue = std::multiset<MetricFace, decltype(metricComp)> (metricComp);

  ProxyLabel label;
  ProxyLabel labelEnd = proxyList.size();

  // For each proxy, first add its seed face to priority queue
  for (label = 0; label < labelEnd; label++) {
    Proxy proxy = proxyList[label];
    FaceIter seedFace = proxy.seed;

    MetricFace m;
    m.face = seedFace;
    m.distance = 0.0;
    m.possibleLabel = label;

    faceQueue.insert(m);
  }

  // For each face in the priority queue, perform flooding by adding its adjacent faces into the priority queue
  for (auto queueIter = faceQueue.begin(); queueIter != faceQueue.end(); ) {
    FaceIter face = queueIter->face;
    if (face->label == -1) { // this means the face hasn't been labeled before
      // label this face, because it is supposed to have the least distortion error
      face->label = queueIter->possibleLabel;

      // add all unlabeled adjacent faces to the priority queue
      HalfedgeIter he = face->halfedge();
      do {
        HalfedgeIter twin = he->twin();
        if (!twin->isBoundary() && twin->face()->label == -1) {
          MetricFace m;
          m.face = twin->face();
          m.distance = calculateDistortionError(m.face, proxyList[queueIter->possibleLabel]);
          m.possibleLabel = queueIter->possibleLabel;
          faceQueue.insert(m);
        }
        he = he->next();
      } while (he != face->halfedge());
    }
    // pop the top of priority queue
    faceQueue.erase(queueIter);
    queueIter = faceQueue.begin();
  }
}

void VSABuilder::fitProxy() {

  // add up related values per proxy
  double *totalArea = new double[numProxies](); ///< allocate and initialize to zero
  double *totalDistortionError = new double[numProxies]();
  double *lowestNormalDiff = new double[numProxies]();
  Vector3D *areaWeightedNormal = new Vector3D[numProxies]();
  Vector3D *areaWeightedCentroid = new Vector3D[numProxies]();
  for (auto f = meshPtr->facesBegin(); f != meshPtr->facesEnd(); f++) {
    ProxyLabel label = f->label;
    double area = faceArea(f);
    totalArea[label]              += area;
    totalDistortionError[label]   += calculateDistortionError(f, proxyList[label]);
    areaWeightedNormal[label]     += f->normal() * area;
    areaWeightedCentroid[label]   += f->centroid() * area;
  }

  // update normal and seed of all proxies
  for (ProxyLabel label = 0; label < numProxies; label++) {
    proxyList[label].totalDistorsion = totalDistortionError[label];
    proxyList[label].totalArea = totalArea[label];
    proxyList[label].normal = areaWeightedNormal[label] / totalArea[label];
    proxyList[label].centroid = areaWeightedCentroid[label] / totalArea[label];
    // initialize with double max value
    lowestNormalDiff[label] = 4.0;
  }

  // We pick the face with the lowest difference of normal as the new seed
  for (auto f = meshPtr->facesBegin(); f != meshPtr->facesEnd(); f++) {
    ProxyLabel label = f->label;
    double normalDiff = (f->normal() - proxyList[label].normal).norm2();
    if (normalDiff < lowestNormalDiff[label]) {
      proxyList[label].seed = f;
      lowestNormalDiff[label] = normalDiff;
    }
    // clear face label
    f->label = -1;
  }

  delete[] totalArea;
  delete[] totalDistortionError;
  delete[] lowestNormalDiff;
  delete[] areaWeightedNormal;
  delete[] areaWeightedCentroid;
}

double VSABuilder::calculateDistortionError(FaceCIter face, Proxy proxy) {
  Vector3D diffNormal = face->normal() - proxy.normal;
  double norm2 = diffNormal.norm2();
  double area = faceArea(face);
  return norm2 * area;
}

double VSABuilder::faceArea(FaceCIter face) {
  // Assuming face is triangular
  auto he = face->halfedge();
  Vector3D P = he->vertex()->position;
  Vector3D Q = he->next()->vertex()->position;
  Vector3D R = he->next()->next()->vertex()->position;

  Vector3D PQ = Q - P;
  Vector3D PR = R - P;
  double area = 0.5 * cross(PQ, PR).norm();

  return area;
}

#ifdef CUDA_VSA
/** IO functions for CUDA VSA */
void VSABuilder::getCudaResult() {
  FaceIndex faceIndex;
  FaceIter f;
  std::vector<FaceIter> faceIndexMap = std::vector<FaceIter>(meshPtr->nFaces());
  for (f = meshPtr->facesBegin(), faceIndex = 0;
       f != meshPtr->facesEnd();
       f++, faceIndex++) {
    faceIndexMap[faceIndex] = f;
    f->label = hostFaceCu[faceIndex].label;
  }

  for (ProxyLabel p = 0; p < numProxies; p++) {
    proxyList[p].label = p;
    proxyList[p].seed = faceIndexMap[hostProxyCu[p].seed];
    proxyList[p].centroid = Vector3D(hostProxyCu[p].centroid.x,
                                     hostProxyCu[p].centroid.y,
                                     hostProxyCu[p].centroid.z);
    proxyList[p].normal = Vector3D(hostProxyCu[p].normal.x,
                                   hostProxyCu[p].normal.y,
                                   hostProxyCu[p].normal.z);
    proxyList[p].borderHalfedge = meshPtr->halfedgesEnd();
  }
}
void VSABuilder::buildFaceCu(){
  hostFaceCu.resize(meshPtr->nFaces());
  std::map<FaceCIter, FaceIndex> faceIndexMap;
  FaceIndex count = 0;
  for (FaceCIter f = meshPtr->facesBegin(); f != meshPtr->facesEnd(); f++) {
    faceIndexMap[f] = count;
    count++;
  }
  for (FaceCIter f = meshPtr->facesBegin(); f != meshPtr->facesEnd(); f++) {
    setFaceNeighbors(f, hostFaceCu, faceIndexMap);
  }
}
void VSABuilder::setFaceNeighbors(FaceCIter f, std::vector<FaceCu> &faceCu,
                                  std::map<FaceCIter, FaceIndex> &faceIndexMap) {
  FaceIndex currentFaceIndex = faceIndexMap[f];
  Vector3D normal = f->normal();
  faceCu[currentFaceIndex].normal = make_double3(normal.x, normal.y, normal.z);
  Vector3D centroid = f->centroid();
  faceCu[currentFaceIndex].centroid = make_double3(centroid.x, centroid.y, centroid.z);
  faceCu[currentFaceIndex].isBoundary = f->isBoundary();
  faceCu[currentFaceIndex].label = f->label;
  faceCu[currentFaceIndex].area = faceArea(f);

  auto h = f->halfedge();
  for (int i = 0; i < 3; i++) {
    auto h_twin = h->twin();
    // check for boundary faces
    if (f->isBoundary() || h_twin->isBoundary()) {
      faceCu[currentFaceIndex].neighbors[i] = -1;
    } else {
      auto neighbor = h_twin->face();
      FaceIndex neighborIndex = faceIndexMap[neighbor];
      faceCu[currentFaceIndex].neighbors[i] = neighborIndex;
    }
    h = h->next();
  }
  if (h != f->halfedge()) fprintf(stderr, "[FaceNeighbors] didn't go through all neighbors.\n");
}
void VSABuilder::writeFaceNeighbors(std::vector<FaceCu> &faceCu, char* fileName, int numberOfFaces) {
	ofstream faceFile;
	faceFile.open(fileName);
	faceFile << numberOfFaces;
	faceFile << "\n";
	for (std::vector<FaceCu>::iterator it = faceCu.begin(); it != faceCu.end(); ++it) {
		faceFile << (*it).normal.x;
		faceFile << "\n";
		faceFile << (*it).normal.y;
		faceFile << "\n";
		faceFile << (*it).normal.z;
		faceFile << "\n";
		faceFile << (*it).centroid.x;
		faceFile << "\n";
		faceFile << (*it).centroid.y;
		faceFile << "\n";
		faceFile << (*it).centroid.z;
		faceFile << "\n";
		faceFile << (*it).isBoundary;
		faceFile << "\n";
		faceFile << (*it).label;
		faceFile << "\n";
		faceFile << (*it).area;
		faceFile << "\n";
		for (int i = 0; i < 3; i++) {
			faceFile << (*it).neighbors[i];
			faceFile << "\n";
		}
		//faceFile << "\n";
	}
	faceFile.close();
}

void VSABuilder::loadFaceNeighbors(std::vector<FaceCu> &faceCu, char* fileName) {
	ifstream faceFile;
	faceFile.open(fileName);
	if (!faceFile) {
		cerr << "Unable to open file datafile.txt";
		exit(1);   // call system to stop
	}
	int numFaces;
	faceFile >> numFaces;
	cout << numFaces;
	cout << "\n";
	hostFaceCu.resize(numFaces);
	//printf("numFaces: %d\n");
	char* line;
	int index = 0;
	int numPropertiesPerFace = 12; //vector3d counted as 3
	double elem;
	double oneFaceProps[12];
	while (!faceFile.eof()) {
		for (int i = 0; i < numPropertiesPerFace; i++) {
			faceFile >> elem;
			oneFaceProps[i] = elem;
		}
		faceCu[index].normal = make_double3(oneFaceProps[0], oneFaceProps[1], oneFaceProps[2]);
		faceCu[index].centroid = make_double3(oneFaceProps[3], oneFaceProps[4], oneFaceProps[5]);
		faceCu[index].isBoundary = oneFaceProps[6];
		faceCu[index].label = oneFaceProps[7];
		faceCu[index].area = oneFaceProps[8];
		faceCu[index].neighbors[0] = oneFaceProps[9];
		faceCu[index].neighbors[1] = oneFaceProps[10];
		faceCu[index].neighbors[2] = oneFaceProps[11];
		index++;
	}
	faceFile.close();
}
#endif //CUDA_VSA

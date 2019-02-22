#ifndef CMU462_MESHEDIT_H
#define CMU462_MESHEDIT_H

#include "halfEdgeMesh.h"
#include "vsa/vsa.h"

using namespace std;

namespace CMU462 {

  class MeshResampler{
    public:

      MeshResampler(){ vsaBuilder = VSABuilder(); };
      ~MeshResampler(){}

      void upsample  ( HalfedgeMesh& mesh );
      void downsample( HalfedgeMesh& mesh );
      void resample(HalfedgeMesh &mesh, size_t numProxy, size_t numIterations, double edgeSplitThreshold);

      VSABuilder vsaBuilder;
  };

} // namespace CMU462

#endif // CMU462_MESHEDIT_H

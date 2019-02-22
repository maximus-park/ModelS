#ifndef CMU462_DYNAMICSCENE_MESH_H
#define CMU462_DYNAMICSCENE_MESH_H

#include "scene.h"

#include "../collada/polymesh_info.h"
#include "../halfEdgeMesh.h"
#include "../meshEdit.h"
#include "skeleton.h"

#include <map>

namespace CMU462 { namespace DynamicScene {

// A structure for holding linear blend skinning information,
// One each vertex of mesh
class LBSInfo {
public:
  // You should compute each vertex's position with respect to each
  // joint j in the skeleton in j's coordinate frame when no
  // transformations have been applied to the skeleton (bind pose).
  // Then, you should find where this vertex would end up with respect to
  // joint j after all of the joints in the skeletons have been
  // transformed if it were to have the same relative position in j's coordinate frame.
  Vector3D blendPos;

  // for this vertex, you should find the closest point on joint j's
  // bone segment (axis) and store the distance to the closest point.
  double distance;
};

class Mesh : public SceneObject {
 public:

  Mesh( Collada::PolymeshInfo& polyMesh,
        const Matrix4x4& transform );

  ~Mesh();

  void set_draw_styles(DrawStyle *defaultStyle, DrawStyle *hoveredStyle,
                       DrawStyle *selectedStyle) override;
  virtual void draw() override;
  virtual void drawGhost() override;

  void draw_pretty() override;

  void draw_normal() override;
  void draw_shading() override;

  StaticScene::SceneObject *get_transformed_static_object(double t) override;

  BBox get_bbox() override;

  virtual Info getInfo() override;

  void _bevel_selection( double inset, double shift );

  virtual void drag( double x, double y, double dx, double dy, const Matrix4x4& modelViewProj ) override;

  BSDF *get_bsdf();
  StaticScene::SceneObject *get_static_object() override;

  void collapse_selected_element();
  void flip_selected_edge();
  void split_selected_edge();
  void erase_selected_element();
  void bevel_selected_element();
  void upsample();
  void downsample();
  void resample(size_t numProxy, size_t numIterations, double edgeSplitThreshold);
  void triangulate();

  HalfedgeMesh mesh;

  Skeleton* skeleton; // skeleton for mesh
  void linearBlendSkinning(bool useCapsuleRadius, double time);
  void forward_euler(float timestep, float damping_factor);
  void symplectic_euler(float timestep, float damping_factor);
  void resetWave();
  void keyframe(double t);
  void unkeyframe(double t);
  void restoreMesh();
  /**
   * Rather than drawing the object geometry for display, this method draws the
   * object with unique colors that can be used to determine which object was
   * selected or "picked" by the cursor.  The parameter pickID is the lowest
   * consecutive integer that has so far not been used by any other object as
   * a picking ID.  (Draw colors are then derived from these IDs.)  This data
   * will be used by Scene::update_selection to make the final determination
   * of which object (and possibly element within that object) was picked.
   */
  virtual void draw_pick( int& pickID, bool transformed = false ) override;

  /** Assigns attributes of the selection based on the ID of the
   * object that was picked.  Can assume that pickID was one of
   * the IDs generated during this object's call to draw_pick().
   */
  virtual void setSelection( int pickID, Selection& selection ) override;

 private:

  // Helpers for draw().
  void draw_faces_normal() const;
  void draw_faces(bool smooth=false) const;
  void draw_faces_lighting() const;
  void draw_edges() const;
  void draw_feature_if_needed( Selection* s ) const;
  void draw_vertex(const Vertex *v) const;
  void draw_halfedge_arrow(const Halfedge *h) const;
  DrawStyle *get_draw_style(const HalfedgeElement *element) const;

  // a vector of halfedges whose vertices are newly created with bevel
  // on scroll, reposition vertices referenced from these halfedges
  vector<HalfedgeIter> bevelVertices;
  // original position of beveled vertex
  Vector3D beveledVertexPos;
  // original positions of beveled edge, corresponding to bevelVertices
  vector<Vector3D> beveledEdgePos;
  // original vertex positions for face currently being beveled
  vector<Vector3D> beveledFacePos;

  DrawStyle *defaultStyle, *hoveredStyle, *selectedStyle;

  MeshResampler resampler;
  
  // map from picking IDs to mesh elements, generated during draw_pick
  // and used by setSelection
  std::vector<HalfedgeElement*> idToElement;
  int idOffset;


  // Assigns the next consecutive pickID to the given element, and
  // sets the GL color accordingly.  Also increments pickID.
  void newPickElement( int& pickID, HalfedgeElement* e );

  // material
  BSDF* bsdf;

  // Compute the distance from a skin vertex to a given joint axis
  double distance_to_bone(Vector3D v_local, Vector3D joint_axis) const;
};

} // namespace DynamicScene
} // namespace CMU462

#endif // CMU462_DYNAMICSCENE_MESH_H

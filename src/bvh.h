#ifndef CMU462_BVH_H
#define CMU462_BVH_H

#include "static_scene/scene.h"
#include "static_scene/aggregate.h"

#include <vector>

namespace CMU462 { namespace StaticScene {


  /**
   * A node in the BVH accelerator aggregate.
   * The accelerator uses a "flat tree" structure where all the primitives are
   * stored in one vector. A node in the data structure stores only the starting
   * index and the number of primitives in the node and uses this information to
   * index into the primitive vector for actual data. In this implementation all
   * primitives (index + range) are stored on leaf nodes. A leaf node has no child
   * node and its range should be no greater than the maximum leaf size used when
   * constructing the BVH.
   */
  struct BVHNode {

    BVHNode(BBox bb, size_t start, size_t range)
      : bb(bb), start(start), range(range), l(NULL), r(NULL) { }

    inline bool isLeaf() const { return l == NULL && r == NULL; }

    BBox bb;        ///< bounding box of the node
    size_t start;   ///< start index into the primitive list
    size_t range;   ///< range of index into the primitive list
    BVHNode* l;     ///< left child node
    BVHNode* r;     ///< right child node
  };

  /**
   * Bounding Volume Hierarchy for fast Ray - Primitive intersection.
   * Note that the BVHAccel is an Aggregate (A Primitive itself) that contains
   * all the primitives it was built from. Therefore once a BVHAccel Aggregate
   * is created, the original input primitives can be ignored from the scene
   * during ray intersection tests as they are contained in the aggregate.
   */
  class BVHAccel : public Aggregate {
    public:

      BVHAccel () { }

      /**
       * Parameterized Constructor.
       * Create BVH from a list of primitives. Note that the BVHAccel Aggregate
       * stores pointers to the primitives and thus the primitives need be kept
       * in memory for the aggregate to function properly.
       * \param primitives primitives to build from
       * \param max_leaf_size maximum number of primitives to be stored in leaves
       */
      BVHAccel(const std::vector<Primitive*>& primitives, size_t max_leaf_size = 4);

      /**
       * Destructor.
       * The destructor only destroys the Aggregate itself, the primitives that
       * it contains are left untouched.
       */
      ~BVHAccel();

      /**
       * Get the world space bounding box of the aggregate.
       * \return world space bounding box of the aggregate
       */
      BBox get_bbox() const;

      /**
       * Ray - Aggregate intersection.
       * Check if the given ray intersects with the aggregate (any primitive in
       * the aggregate), no intersection information is stored.
       * \param r ray to test intersection with
       * \return true if the given ray intersects with the aggregate,
       false otherwise
       */
      bool intersect(const Ray& r) const;

      /**
       * Ray - Aggregate intersection 2.
       * Check if the given ray intersects with the aggregate (any primitive in
       * the aggregate). If so, the input intersection data is updated to contain
       * intersection information for the point of intersection. Note that the
       * intersected primitive entry in the intersection should be updated to
       * the actual primitive in the aggregate that the ray intersected with and
       * not the aggregate itself.
       * \param r ray to test intersection with
       * \param i address to store intersection info
       * \return true if the given ray intersects with the aggregate,
       false otherwise
       */
      bool intersect(const Ray& r, Intersection* i) const;

      /**
       * Get BSDF of the surface material
       * Note that this does not make sense for the BVHAccel aggregate
       * because it does not have a surface material. Therefore this
       * should always return a null pointer.
       */
      BSDF* get_bsdf() const { return NULL; }

      /**
       * Get entry point (root) - used in visualizer
       */
      BVHNode* get_root() const { return root; }

      /**
       * Draw the BVH with OpenGL - used in visualizer
       */
      void draw(const Color& c) const { }

      /**
       * Draw the BVH outline with OpenGL - used in visualizer
       */
      void drawOutline(const Color& c) const { }
      /**
       * Sort the sub-primitive-list of a particular BVHNode
       * along the longest edge of node's BBox.
       */
      void sortPrimitives(BVHNode* node) {

	size_t start, range;
	start = node->start; range = node->range;

	// create sub-vector for sorting
	std::vector<Primitive*> sub(primitives.begin() + start,
				    primitives.begin() + start + range);
	// sort using helper compare functions
	double x = node->bb.extent.x,
	       y = node->bb.extent.y,
	       z = node->bb.extent.z;
	if (x > y && x > z)
	  std::sort(sub.begin(), sub.end(), compareX);
	else if (y > x && y > z)
	  std::sort(sub.begin(), sub.end(), compareY);
	else
	  std::sort(sub.begin(), sub.end(), compareZ);
	// copy sub back to primitives
	for (int i = 0; i < range; ++i){
	    primitives[start + i] = sub[i];
	}
      }

    

    private:
      BVHNode* root; ///< root node of the BVH
      std::vector<BVHNode*> nodes; ///< all nodes except root of the BVH

      /**
       * Recursively split a BVHNode, until the primitive count
       * of every node is less than or equal to max_leaf_size
       */
      void splitNode(BVHNode* node, size_t max_leaf_size) {

        if (node->range <= max_leaf_size) return;

        sortPrimitives(node);

        size_t splt = findSplitPoint(node);
        // left child node BBox
        BBox bbl;
        for (size_t i = node->start; i < splt; ++i) {
            bbl.expand(primitives[i]->get_bbox());
        }
        // right child node BBox
        BBox bbr;
        for (size_t i = splt; i < node->start + node->range; ++i) {
            bbr.expand(primitives[i]->get_bbox());
        }
        // build child nodes
        node->l = new BVHNode(bbl, node->start, splt - node->start);
        node->r = new BVHNode(bbr, splt, node->range - (splt - node->start));
        // record new nodes
        nodes.push_back(node->l);
        nodes.push_back(node->r);

        // Recursively split left child and right child
        splitNode(node->l, max_leaf_size);
        splitNode(node->r, max_leaf_size);

      }
      /**
       * Find split point in a node that has least cost, using
       * Surface Area Heuristic method.
       * Assuming that cost is a convex function, so it may return
       * a local minimum cost.
       */
      size_t findSplitPoint(BVHNode* node) const{
	size_t splt = node->start + node->range / 2;
	size_t left_range  = splt - node->start;
	size_t right_range = node->range - left_range;
	double sah, sah_min;
	BBox bbl, bbr;
	// left child node BBox
	for (size_t i = node->start; i < splt; ++i) {
	    bbl.expand(primitives[i]->get_bbox());
	}
	// right child node BBox
	for (size_t i = splt; i < node->start + node->range; ++i) {
	    bbr.expand(primitives[i]->get_bbox());
	}
	sah_min = left_range  * bbl.surface_area()
		+ right_range * bbr.surface_area();

	// splt moves left to see if sah decreases
	while(left_range > 1) {
	    splt--; left_range--; right_range++;
	    // left BBox shrinks
	    bbl = BBox();
	    for (size_t i = node->start; i < splt; ++i) {
		  bbl.expand(primitives[i]->get_bbox());
	    }
	    // right BBox expands
	    bbr.expand(primitives[splt]->get_bbox());
	    sah = left_range  * bbl.surface_area()
		+ right_range * bbr.surface_area();
	    if (sah > sah_min) break;
	    else sah_min = sah;
	}

	if (splt != (node->start + node->range / 2) - 1){
	    return (splt + 1);
	}
	// splt moves right to see if sah decreases
	while(right_range > 1) {
	    left_range++; right_range--;
	    // left BBox expands
	    bbl.expand(primitives[splt]->get_bbox());
	    splt++;
	    // right BBox shrinks
	    bbr = BBox();
	    for (size_t i = splt; i < node->start + node->range; ++i) {
		bbr.expand(primitives[i]->get_bbox());
	    }
	    sah = left_range  * bbl.surface_area()
		+ right_range * bbr.surface_area();
	    if (sah > sah_min){
		return (splt - 1);
	    }
	    else sah_min = sah;
	}
	return splt;
      }

      /**
       * Helper compare functions for std::sort()
       */
      static bool compareX(const Primitive* p1, const Primitive* p2) {
	  return (p1->get_bbox().centroid().x < p2->get_bbox().centroid().x);
      }
      static bool compareY(const Primitive* p1, const Primitive* p2){
	  return (p1->get_bbox().centroid().y < p2->get_bbox().centroid().y);
      }
      static bool compareZ(const Primitive* p1, const Primitive* p2){
	  return (p1->get_bbox().centroid().z < p2->get_bbox().centroid().z);
      }

      /**
       * Helper function that uses "front-to-back" traversal to find the
       * closest hit primitive in a BVH tree. The first one only does the
       * query, the second one updates Intersection struct.
       */
      bool findClosestHit(const Ray &ray, BVHNode* node) const {

	// In leaf node, check every primitive for intersection
	if (node->isLeaf()) {
	    size_t start = node->start,
		   range = node->range;
	    for (int i = 0; i < range; ++i){
		if(primitives[i + start]->intersect(ray)) return true;
	    }
	    return false;
	}
	double t0, t1;
	if (!node->bb.intersect(ray, t0, t1)) return false;
	// Normal case
	else{
	    double t0l, t0r, t1l, t1r;
	    bool hit_l = node->l->bb.intersect(ray, t0l, t1l),
		 hit_r = node->r->bb.intersect(ray, t0l, t1l);
	    if (hit_l && hit_r) {
		if (t0l < t0r) {
		    if (findClosestHit(ray, node->l)) return true;
		    if (findClosestHit(ray, node->r)) return true;
		}
		else{
		    if (findClosestHit(ray, node->r)) return true;
		    if (findClosestHit(ray, node->l)) return true;
		}
	    }
	    else {
		if      (hit_l){
		    if (findClosestHit(ray, node->l)) return true;
		}
		else if (hit_r){
		    if (findClosestHit(ray, node->r)) return true;
		}
	    }
	    return false;
	}
      }

      bool findClosestHit(const Ray &ray, BVHNode* node, Intersection* isect) const {

	// In leaf node, check every primitive for intersection
	if (node->isLeaf()) {
	    size_t start = node->start,
		   range = node->range;
	    bool hit = false;
	    for (int i = 0; i < range; ++i){
		if(primitives[i + start]->intersect(ray, isect)) hit = true;
	    }
	    return hit;
	}

	// Normal case
	else{
	    double t0l, t0r, t1l, t1r;
	    bool found_l, found_r;
	    bool hit_l = node->l->bb.intersect(ray, t0l, t1l),
		 hit_r = node->r->bb.intersect(ray, t0l, t1l);
	    if (hit_l){
		found_l = findClosestHit(ray, node->l, isect);
	    }
	    if (hit_r){
		found_r = findClosestHit(ray, node->r, isect);
	    }
	    return found_l || found_r;
	}
      }
  };

} // namespace StaticScene
} // namespace CMU462

#endif // CMU462_BVH_H

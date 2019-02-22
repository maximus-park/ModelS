#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  // r.d should be normalized direction

  Vector3D o_rel = r.o - this->o;
  double o_dot_d = dot(o_rel, r.d);
  double delta = o_dot_d * o_dot_d - o_rel.norm2() + this->r2;
  if (delta < 0) return false;
  else {
      double sqrt_delta = sqrt(delta);
      t1 = -o_dot_d - sqrt_delta;
      t2 = -o_dot_d + sqrt_delta;
      return true;
  }
}

bool Sphere::intersect(const Ray& r) const {

  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  if (test(r, t1, t2) && ((r.min_t <= t1 && t1 <= r.max_t) ||
	 (r.min_t <= t2 && t2 <= r.max_t))) {
      return true;
  }
  return false;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1, t2;
  if (test(r, t1, t2)) {
      // Check whether it is the closest hit, do we have to do that???
      if (r.min_t <= t1 && t1 <= r.max_t){
	  if (t1 < i->t) {
	      i->t = t1;
	      Vector3D normal = (r.o + t1*r.d - this->o) / this->r;
	      i->n = normal;
	      i->primitive = this;
	      i->bsdf = object->get_bsdf();
	      return true;
	  }
      }
      // Check whether it is the closest hit, do we have to do that???
      else if (r.min_t <= t2 && t2 <= r.max_t) {
	  if (t2 < i->t) {
	      i->t = t2;
	      Vector3D normal = (r.o + t2*r.d - this->o) / this->r;
	      i->n = normal; // Mustn't be negetive, because of refraction purposes
	      i->primitive = this;
	      i->bsdf = object->get_bsdf();
	      i->is_back_hit = true;
	      return true;
	  }
      }
  }
  return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CMU462

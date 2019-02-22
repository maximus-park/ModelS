#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, vector<size_t>& v) :
    mesh(mesh), v(v) { }
Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {
  
  // compute the bounding box of the triangle
  Vector3D p0 = mesh->positions[v1],
  	   p1 = mesh->positions[v2],
  	   p2 = mesh->positions[v3];
  double minX = min(p0.x, min(p1.x, p2.x)),
	 minY = min(p0.y, min(p1.y, p2.y)),
	 minZ = min(p0.z, min(p1.z, p2.z)),
	 maxX = max(p0.x, max(p1.x, p2.x)),
	 maxY = max(p0.y, max(p1.y, p2.y)),
	 maxZ = max(p0.z, max(p1.z, p2.z));
  return BBox(minX, minY, minZ, maxX, maxY, maxZ);
}

bool Triangle::intersect(const Ray& r) const {
  
  // ray-triangle intersection
  // Algorithm from Assignment 3 Notes
  Vector3D p0 = mesh->positions[v1],
	   p1 = mesh->positions[v2],
	   p2 = mesh->positions[v3];
  Vector3D e1 = p1 - p0, e2 = p2 - p0, s = r.o - p0;
  double detM = dot(cross(e1, r.d), e2);
  if (detM == 0) return false;
  double detMu = -dot(cross(s,  e2 ), r.d),
         detMv =  dot(cross(e1, r.d),  s ),
	 detMt = -dot(cross(s,  e2 ),  e1);
  double u = detMu / detM, v = detMv / detM, t = detMt / detM;
  return (0 <= u && u <= 1) && (0 <= v && v <= 1) && (0 <= u+v && u+v <= 1)
      && (r.min_t <= t && t <= r.max_t);
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  Vector3D p0 = mesh->positions[v1],
	   p1 = mesh->positions[v2],
	   p2 = mesh->positions[v3];
  Vector3D e1 = p1 - p0, e2 = p2 - p0, s = r.o - p0;
  double detM = dot(cross(e1, r.d), e2);
  if (detM == 0) return false;
  double detMu = -dot(cross(s,  e2 ), r.d),
         detMv =  dot(cross(e1, r.d),  s ),
	 detMt = -dot(cross(s,  e2 ),  e1);
  double u = detMu / detM, v = detMv / detM, t = detMt / detM;
  bool intersect = (0 <= u && u <= 1) && (0 <= v && v <= 1)
		&& (0 <= u+v && u+v <= 1) && (r.min_t <= t && t <= r.max_t);

  // Check whether it is the closest hit
  if (intersect && t < isect->t){
      isect->t = t;
      isect->n = mesh->normals[v1] * (1-u-v) +
		 mesh->normals[v2] * (u) +
		 mesh->normals[v3] * (v);
      if (dot(r.d, isect->n) > 0) {
	  isect->is_back_hit = true;
	  isect->n = -isect->n;
      }
      isect->primitive = this;
      isect->bsdf = mesh->get_bsdf();
      return true;
  }
  else{
      return false;
  }
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CMU462

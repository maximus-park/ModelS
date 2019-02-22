#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CMU462 {

  bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

    // Implement ray - bounding box intersection test
    // If the ray intersected the bouding box within the range given by
    // t0, t1, update t0 and t1 with the new intersection times.
    double minX = this->min.x, minY = this->min.y, minZ = this->min.z,
	   maxX = this->max.x, maxY = this->max.y, maxZ = this->max.z;
    bool zeroX = (r.d.x == 0), zeroY = (r.d.y == 0), zeroZ = (r.d.z == 0);
    if (zeroX) {
	if(maxX < r.o.x || r.o.x < minX) return false;
    }
    else if (zeroY) {
	if(maxY < r.o.y || r.o.y < minY) return false;
    }
    else if (zeroZ) {
	if(maxZ < r.o.z || r.o.z < minZ) return false;
    }

    double t_minX = (zeroX) ? -INF_D : (minX - r.o.x) / r.d.x;
    double t_maxX = (zeroX) ?  INF_D : (maxX - r.o.x) / r.d.x;
    if (t_minX > t_maxX) std::swap(t_minX, t_maxX);
    double t_minY = (zeroY) ? -INF_D : (minY - r.o.y) / r.d.y;
    double t_maxY = (zeroY) ?  INF_D : (maxY - r.o.y) / r.d.y;
    if (t_minY > t_maxY) std::swap(t_minY, t_maxY);
    double t_minZ = (zeroZ) ? -INF_D : (minZ - r.o.z) / r.d.z;
    double t_maxZ = (zeroZ) ?  INF_D : (maxZ - r.o.z) / r.d.z;
    if (t_minZ > t_maxZ) std::swap(t_minZ, t_maxZ);

    t0 = std::max(t_minX, std::max(t_minY, t_minZ));
    t1 = std::min(t_maxX, std::min(t_maxY, t_maxZ));

    if (t0 <= t1) return true;
    return false;
  }

  void BBox::draw(Color c) const {

    glColor4f(c.r, c.g, c.b, c.a);

    // top
    glBegin(GL_LINE_STRIP);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(max.x, max.y, max.z);
    glEnd();

    // bottom
    glBegin(GL_LINE_STRIP);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, min.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glEnd();

    // side
    glBegin(GL_LINES);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(min.x, min.y, max.z);
    glEnd();

  }

  std::ostream& operator<<(std::ostream& os, const BBox& b) {
    return os << "BBOX(" << b.min << ", " << b.max << ")";
  }

} // namespace CMU462

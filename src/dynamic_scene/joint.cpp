/*
 * Implementations for Joint Based Skeletons.
 *
 * Started on October 29th, 2015 by Bryce Summers.
 */

#include "joint.h"
#include "skeleton.h"
#include "mesh.h"

#include "GL/glew.h"

#include <iostream>

namespace CMU462 { namespace DynamicScene {

   BBox Joint::get_bbox() {
      BBox bbox;
      Vector3D p1 = position;
      Vector3D p2 = p1 + axis;
      bbox.expand(p1);
      bbox.expand(p2);
      return bbox;
   }

   Info Joint::getInfo()
   {
      Info info;

      if (!scene || !scene->selected.element)
      {
         info.push_back("JOINT");
      }
      else
      {
         info = scene->selected.element->getInfo();
      }

      return info;
   }

   void Joint::drag(double x, double y, double dx, double dy, const Matrix4x4& modelViewProj)
   {
      Vector4D q(position, 1.);

      // Transform into clip space
      q = modelViewProj * q;
      double w = q.w;
      q /= w;

      // Shift by (dx, dy).
      q.x += dx;
      q.y += dy;

      // Transform back into model space
      q *= w;
      q = modelViewProj.inv() * q;

      if (skeleton->root == this)
         position = q.to3D();
      else
         axis += q.to3D();
   }

   StaticScene::SceneObject *Joint::get_static_object() {
      return nullptr;
   }

   // The real calculation.
   void Joint :: calculateAngleGradient( Joint* goalJoint, Vector3D q )
   {
     // Implement Me! (task 2B)
     Vector3D p_theta = goalJoint->getEndPosInWorld();
     Vector3D q_minus_p = q - p_theta; // World Coordinate q-p

     Vector3D p = p_theta - this->getBasePosInWorld();
     vector<Vector3D> axes;
     this->getAxes(axes);
     // Rotation axes in world
     Vector3D axis_X_World = axes[0],
	      axis_Y_World = axes[1],
	      axis_Z_World = axes[2];
     // Jacobian Matrix for x, y, z rotations
     Vector3D Jx = cross(axis_X_World, p),
	      Jy = cross(axis_Y_World, p),
	      Jz = cross(axis_Z_World, p);

     // Update ikAngleGradient, should use += for multi-target purposes
     this->ikAngleGradient[0] += dot(Jx, q_minus_p);
     this->ikAngleGradient[1] += dot(Jy, q_minus_p);
     this->ikAngleGradient[2] += dot(Jz, q_minus_p);
   }


   // The constructor sets the dynamic angle and velocity of
   // the joint to zero (at a perfect vertical with no motion)
   Joint :: Joint(Skeleton * s)
   : skeleton( s ), capsuleRadius( 0.05 ), renderScale( 1.0 )
   {
     scale = Vector3D(1., 1., 1.);
     scales.setValue(0, scale);
   }

   Vector3D Joint::getAngle( double time )
   {
      return rotations(time);
   }

   void Joint::setAngle( double time, Vector3D value )
   {
      rotations.setValue( time, value );
   }

   bool Joint::removeAngle(double time)
   {
      return rotations.removeKnot(time, .1);
   }

   void Joint::keyframe(double t) {
     positions.setValue(t, position);
     rotations.setValue(t, rotation);
     scales.setValue(t, scale);
     for (Joint *j : kids) j->keyframe(t);
   }

   void Joint::unkeyframe(double t) {
     positions.removeKnot(t, 0.1);
     rotations.removeKnot(t, 0.1);
     scales.removeKnot(t, 0.1);
     for (Joint *j : kids) j->unkeyframe(t);
   }

   void Joint::removeJoint(Scene* scene)
   {
     if (this == skeleton->root)
       return;

     for (auto childJoint : kids)
     {
       childJoint->removeJoint(scene);
     }

     scene->removeObject(this);

     auto & kids = parent->kids;
     kids.erase(std::remove(kids.begin(), kids.end(), this), kids.end());

     auto & joints = skeleton->joints;
     joints.erase(std::remove(joints.begin(), joints.end(), this), joints.end());

     delete this;
   }

   void Joint::getAxes(vector<Vector3D>& axes)
   {
     Matrix4x4 T = Matrix4x4::identity();
     for (Joint* j = parent; j != nullptr; j = j->parent)
     {
       T = j->getRotation() * T;
     }
     T = skeleton->mesh->getRotation() * T;
     axes.resize(3);
     axes[0] = T * Vector3D(1., 0., 0.);
     axes[1] = T * Vector3D(0., 1., 0.);
     axes[2] = T * Vector3D(0., 0., 1.);
   }

   Matrix4x4 Joint::getTransformation()
   {
     /* Implement Me! (Task 2a)
     Initialize a 4x4 identity transformation matrix. Traverse the hierarchy starting from the parent of 
     this joint all the way up to the root (root has parent of nullptr) and accumulate their transformations 
     on the left side of your transformation matrix. Don't forget to apply a translation which extends along 
     the axis of those joints. Finally, apply the mesh's transformation at the end.
     */
     Matrix4x4 T = Matrix4x4::identity();
     for (Joint* j = parent; j != nullptr; j = j->parent)
     {
	 T = Matrix4x4::translation(j->axis) * T;
         T = j->SceneObject::getTransformation() * T;
     }
     T = skeleton->mesh->getTransformation() * T;
     return T;
   }

   Matrix4x4 Joint::getBindTransformation()
   {
     Matrix4x4 T = Matrix4x4::identity();
     for (Joint* j = parent; j != nullptr; j = j->parent)
     {
       T = Matrix4x4::translation(j->axis) * T;
     }

     // Allow skeleton translation by taking root's position into account
     T = Matrix4x4::translation(skeleton->root->init_position) * T;

     return T;
   }

   Vector3D Joint::getBasePosInWorld()
   {
     /* Implement Me! (Task 2a)
     This should be fairly simple once you implement Joint::getTransform().
     You can utilize the transformation returned by Joint::getTransform() to 
     compute the base position in world coordinate frame.
     */
     Vector3D ret = this->getTransformation() * Vector3D();
     return this->getTransformation() * Vector3D();
   }

   Vector3D Joint::getEndPosInWorld()
   {
     /* Implement Me! (Task 2a)
     In addition to what you did for getBasePosInWorld(), you need to apply this joint's 
     transformation and translate along this joint's axis to get the end position in world 
     coordinate frame.
     */
     return this->getTransformation() *
	    this->SceneObject::getTransformation() *
	    Matrix4x4::translation(this->axis) * Vector3D();
   }
} // namespace DynamicScene
} // namespace CMU462

#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CMU462 {

  void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
  }

  // Diffuse BSDF //

  Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    return albedo * (1.0 / PI);
  }

  Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // In Lambertian Diffusion, each sample direction has equal probability

    // *wi is in local coordinate, should be transformed by caller function
    *wi = sampler.get_sample(pdf);

    return Spectrum();
  }

  // Mirror BSDF //

  Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {

    Vector3D r;
    reflect(wo, &r);
    if (r == wi) return Spectrum(1.f, 1.f, 1.f);
    return Spectrum();
  }

  Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // Implement MirrorBSDF
    reflect(wo, wi);
    *pdf = 1.f;
    return Spectrum();
  }

  // Glossy BSDF //

  /*
     Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
     return Spectrum();
     }

     Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
   *pdf = 1.0f;
   return reflect(wo, wi, reflectance);
   }
   */

  // Refraction BSDF //

  Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    Vector3D r;
    refract(wo, &r, ior);
    if (r == wi) return Spectrum(1.f, 1.f, 1.f);
    return Spectrum();
  }

  Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // Implement RefractionBSDF
    if(!refract(wo, wi, ior)) reflect(wo, wi);
    *pdf = 1.f;
    return Spectrum();
  }

  // Glass BSDF //

  Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {

    Vector3D refr;
    refract(wo, &refr, ior);
    Vector3D refl;
    reflect(wo, &refl);

    if (refr == wi || refl == wi) {
	return Spectrum(1.f, 1.f, 1.f);
    }
    return Spectrum();
  }

  Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // Compute Fresnel coefficient and either reflect or refract based on it.
    bool refrc = refract(wo, wi, ior);
    // Total internal reflection
    if(!refrc){
	*pdf = 1.f;
	reflect(wo, wi);
    }
    // Using Fresnel Equations
    else{
	double rp, rv;
	double cos_theta_i = fabs(wo.z);
	double cos_theta_t = fabs(wi->z);
	double ior_i, ior_t;
	if (wo.z > 0) ior_i = 1., ior_t = (double) ior;
	else ior_i = (double) ior, ior_t = 1.;
	// Fresnel equations
	rp = (ior_t*cos_theta_i - ior_i*cos_theta_t)
	   / (ior_t*cos_theta_i + ior_i*cos_theta_t);
	rv = (ior_i*cos_theta_i - ior_t*cos_theta_t)
	   / (ior_i*cos_theta_i + ior_t*cos_theta_t);
	// the fraction of light reflected
	double Fr = (rp*rp + rv*rv) / 2;

	// Random generator [0, 1]
	double rand = (double) (std::rand()) / RAND_MAX;
	if (rand < Fr) reflect(wo, wi);
	*pdf = cos_theta_i / ((ior_t*ior_t)/(ior_i*ior_i) * (1-Fr));
    }
    return Spectrum();
  }

  void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

    // Implement reflection of wo about normal (0,0,1) and store result in wi.
    *wi = Vector3D(-wo.x, -wo.y, wo.z);
  }

  bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

    // Use Snell's Law to refract wo surface and store result ray in wi.
    // Return false if refraction does not occur due to total internal reflection
    // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
    // ray entering the surface through vacuum.

    Vector2D theta1 = Vector2D(wo.x, wo.y);
    double sin_theta1 = theta1.norm();

    // perpendicular light
    if (sin_theta1 == 0) {
	*wi = -wo;
	return true;
    }
    // total internal reflection
    if (wo.z < 0) {
	ior = 1.f / ior;
	if (sin_theta1 > ior) return false;
    }

    double sin2_theta1 = theta1.norm2();
    double n1_over_n2_sq = (1.f/ior) * (1.f/ior);
    double combine = sin2_theta1 * n1_over_n2_sq;
    double zi = -wo.z * sqrt((1-combine) / (n1_over_n2_sq - combine));
    *wi = Vector3D(-wo.x, -wo.y, zi);
    wi->normalize();
    return true;

  }

  // Emission BSDF //

  Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    return Spectrum();
  }

  Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    *wi = sampler.get_sample(pdf);
    return Spectrum();
  }

} // namespace CMU462

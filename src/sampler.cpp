#include "sampler.h"

namespace CMU462 {

  // Uniform Sampler2D Implementation //

  Vector2D UniformGridSampler2D::get_sample() const {

    // Implement uniform 2D grid sampler

    double Xi1 = (double)(std::rand()) / RAND_MAX;
    double Xi2 = (double)(std::rand()) / RAND_MAX;

    return Vector2D(Xi1, Xi2);

  }

  // Uniform Hemisphere Sampler3D Implementation //

  Vector3D UniformHemisphereSampler3D::get_sample() const {

    double Xi1 = (double)(std::rand()) / RAND_MAX;
    double Xi2 = (double)(std::rand()) / RAND_MAX;

    double theta = acos(Xi1);
    double phi = 2.0 * PI * Xi2;

    double xs = sinf(theta) * cosf(phi);
    double ys = sinf(theta) * sinf(phi);
    double zs = cosf(theta);

    return Vector3D(xs, ys, zs);

  }

  Vector3D CosineWeightedHemisphereSampler3D::get_sample() const {
    float f;
    return get_sample(&f);
  }

  Vector3D CosineWeightedHemisphereSampler3D::get_sample(float *pdf) const {

    // Concentric Disk Sampler
    double Xi1 = (double)(std::rand()) / RAND_MAX;
    double Xi2 = (double)(std::rand()) / RAND_MAX;
    // Map Xi1 Xi2 to [-1, 1]
    double o1 = 2. * Xi1 - 1., o2 = 2. * Xi2 - 1.;
    // Handle point at the origin
    if (o1 == 0 && o2 == 0) {
	*pdf = (float) (1.f / PI);
	return Vector3D(0, 0, 1);
    }

    double theta, r;
    if (fabs(o1) > fabs(o2)) {
	r = o1;
	theta = PI/4.* (o1 / o2);
    }
    else {
	r = o2;
	theta = PI/2.- PI/4.* (o1 / o2);
    }

    // Cosine Weighted Hemisphere Sampler
    double xs = r * cosf(theta);
    double ys = r * sinf(theta);
    double zs = sqrt(std::max(0., 1. - xs*xs - ys*ys));

    *pdf = (float) (zs / PI);
    return Vector3D(xs, ys, zs);
  }

} // namespace CMU462

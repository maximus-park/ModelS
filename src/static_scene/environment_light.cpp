#include "environment_light.h"

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {

    // initialize things here as needed
    std::vector<double> illu;
    double total_illu = 0;
    for (size_t i = 0; i < envMap->h; i++){
	double theta = (double) i / (double) envMap->h;
	for (size_t j = 0; j < envMap->w; j++){
	    double il = envMap->data[j + i*envMap->w].illum()*sin(theta);
	    illu.push_back(il);
	    total_illu += il;
	}
    }
    prob.resize(illu.size());
    prob_d.resize(illu.size());

    // prob is a transpose matrix of illu
    double sum_prob = 0;
    for (size_t i = 0; i < envMap->h; i++){
	for (size_t j = 0; j < envMap->w; j++){
	    double curr_prob = illu[j*envMap->h + i] / total_illu;
	    prob_d[j + i*envMap->w] = curr_prob;
	    sum_prob += curr_prob;
	    prob[j + i*envMap->w] = sum_prob;
	}
    }
}

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // env light comes from infinity
  *distToLight = std::numeric_limits<float>::infinity();

  // importance sampling
  double Xi = (double)(std::rand()) / RAND_MAX;
  size_t idx = search_closest_prob(0, prob.size(), Xi);

  // pixel idx to Vec3D
  size_t sx = idx / envMap->h;
  size_t sy = idx % envMap->h;

  double theta =     PI * ((double) sy / (double) envMap->h);
  double phi = 2.0 * PI * ((double) sx / (double) envMap->w);

  double xs = sinf(theta) * cosf(phi);
  double ys = sinf(theta) * sinf(phi);
  double zs = cosf(theta);
  *wi = Vector3D(xs, ys, zs);

  // prob to pdf
  *pdf = prob_d[idx] * (float) (envMap->w * envMap->h);
  if (*pdf < 0.0001f) *pdf = 0.0001f; // pdf can't be equal to 0.f

  return envMap->data[sx + sy*envMap->w];

}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {

  // Get theta and phi
  double theta = acos(r.d.y);
  Vector3D xz = Vector3D (r.d.x, r.d.z, 0);
  xz.normalize();
  double phi = xz.y > 0 ? acos(xz.x) : 2*PI - acos(xz.x);
  // get the sample point on texture
  double x = (double) envMap->w * phi / (2*PI),
	 y = (double) envMap->h * theta / PI;
  int sx1 = (int) floor(x - 0.5),
      sy1 = (int) floor(y - 0.5),
      sx2 = sx1++,
      sy2 = sy1++;
  float  tx1 = (float) sx2 - x + 0.5,
	 ty1 = (float) sy2 - y + 0.5,
	 tx2 = x - 0.5 - (float) sx1,
	 ty2 = y - 0.5 - (float) sy1;
  // get 4 adjacent texel samples
  if (x - 0.5 < 0 || x + 0.5 > (double)(envMap->w)) {sx1 = envMap->w - 1; sx2 = 0;}
  if (y - 0.5 < 0 || y + 0.5 > (double)(envMap->h)) {sy1 = envMap->h - 1; sy2 = 0;}

  // Bilinear interpolation
  Spectrum cuv = envMap->data[sx1+sy1*envMap->w]* tx1 * ty1 +
		 envMap->data[sx1+sy2*envMap->w]* tx1 * ty2 +
		 envMap->data[sx2+sy1*envMap->w]* tx2 * ty1 +
		 envMap->data[sx2+sy2*envMap->w]* tx2 * ty2;
  return cuv;
}

} // namespace StaticScene
} // namespace CMU462

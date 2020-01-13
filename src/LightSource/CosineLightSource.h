/* 
* Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom
* the Software is furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _COSINE_LIGHT_SOURCE_H
#define _COSINE_LIGHT_SOURCE_H

#include "bunnykiller.h"

#include <cmath>

#include "LightSource/PointLightSource.h"
#include "RayTracing/World.h"

/** CosineLightSource class represents a point light emiting 
	light with a cosine distribution. */
template<unsigned int D, class Radiance>
class CosineLightSource : public PointLightSource<D, Radiance>
{
protected:
	typedef PointLightSource<D, Radiance> PLS;

	VectorN<D> direction;
public:
	CosineLightSource(World<D, Radiance>* _world, const VectorN<D>& pos, const VectorN<D>& dir, const Radiance& ints, Real time = 0.) :
		PointLightSource<D, Radiance>(_world, pos, ints, time), direction(dir) {}
	CosineLightSource(World<D, Radiance>* _world, const Intersection<D> &it, const Radiance& ints) :
		PointLightSource<D, Radiance>(_world, it.get_position(),
						ints*it.material()->get_absorption(it.get_uv())*dot_clamped(it.get_normal(), -it.get_ray().get_direction()), 
						it.get_ray().get_parameter()/it.get_ray().get_refraction_index()), direction(it.get_normal()) {}

	CosineLightSource(World<D, Radiance>* _world, const Intersection<D> &it, const Radiance& ints, Real time) :
		PointLightSource<D, Radiance>(_world, it.get_position(),
				ints*it.material()->get_absorption(it.get_uv())*dot_clamped(it.get_normal(), -it.get_ray().get_direction()),
				time), direction(it.get_normal()) {}

	CosineLightSource(World<D, Radiance>* _world, const Ray<D> &ray, const Radiance& ints) :
		PointLightSource<D, Radiance>(_world, ray.get_origin(), ints)
	{
		Ray<D> r(ray);
		Intersection<D> it;
		this->world->first_intersection(r, it);

		if (!it.did_hit())
			throw std::runtime_error("Impossible to create virtual light!");

		PLS::position = it.get_position() + it.get_normal()*1.e-4;
		direction = it.get_normal();

		PLS::intensities = ints * (1.0 - it.material()->get_absorption(it));

		Medium<D> *m = r.get_medium();
		if (m) {
			PLS::intensities *= m->get_transmittance(r);
		}

		PLS::time = PLS::world->time_of_flight(r.get_parameter()) / r.get_ior();
	}

	virtual ~CosineLightSource() {}
  
	Radiance get_incoming_light(const VectorN<D> &point_lighted) const;
	bool sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const;
	void sample(LightSample<D, Radiance> &light_sample, Real &pdf)const;
	bool sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const;

	virtual void print(FILE *_f_log)const;
}; //CosineLightSource

template<unsigned int D, class Radiance>
Radiance CosineLightSource<D, Radiance>::get_incoming_light(const VectorN<D> &point_lighted) const
{
	VectorN<D> dir = point_lighted - PLS::position;
	Real inv_distance2 = 1. / dir.length2();
	Real inv_distance = sqrt(inv_distance2);

	dir *= inv_distance;

	Real att(0.);
	if (D == 2)
		att = inv_distance;
	
	if (D == 3)
		att = inv_distance2;

	return PLS::intensities*att*dot_clamped(direction, dir);
}

template<unsigned int D, class Radiance>
bool CosineLightSource<D, Radiance>::sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLS::position;
	light_sample.dir = it.get_position() - PLS::position;
	Real distance2 = light_sample.dir.length2();
	light_sample.dist = sqrt(distance2);
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = PLS::time;

	pdf = 1.;

	Real dot_light = dot_clamped(light_sample.dir, direction);

	if (it.get_normal().dot(-light_sample.dir) > 0. && dot_light>0. && PLS::is_visible(it.get_position()))
	{
		Real att(0.);
		if (D == 2)
			att = 1. / light_sample.dist;

		if (D == 3)
			att = 1. / distance2;

		light_sample.irradiance = PLS::intensities * att *dot_light;
		return true;
	}

	return false;
}

template<unsigned int D, class Radiance>
bool CosineLightSource<D, Radiance>::sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLS::position;
	light_sample.dir = p - PLS::position;
	Real distance2 = light_sample.dir.length2();
	light_sample.dist = sqrt(distance2);
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = PLS::time;

	pdf = 1.;

	Real dot_light = dot_clamped(light_sample.dir, direction);

	if (dot_light>0. && PLS::is_visible(p))
	{
		Real att(0.);
		if (D == 2)
			att = 1. / light_sample.dist;

		if (D == 3)
			att = 1. / distance2;

		light_sample.irradiance = PLS::intensities * att *dot_light;
		return true;
	}

	return false;
}

template<unsigned int D, class Radiance>
void CosineLightSource<D, Radiance>::sample(LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLS::position;
	light_sample.dist = 0.;
	light_sample.irradiance = PLS::intensities;
	light_sample.instant = PLS::time;

	// Sample semicircle
	Sampling.cosine_weighted(light_sample.dir, pdf);

	light_sample.dir.transform_matrix_to(VectorN<D>(0, 1), direction);
}

template<unsigned int D, class Radiance>
void CosineLightSource<D, Radiance>::print(FILE* f_log)const
{
	fprintf(f_log, "Cosine Light Source:\n");
	fprintf(f_log, "- Position: [%f %f %f]\n", PLS::position[0], PLS::position[1], PLS::position[2]);
	fprintf(f_log, "- Normal: [%f %f %f]\n", direction[0], direction[1], direction[2]);
	fprintf(f_log, "- Time: %f\n", PLS::time);
	Vector3 rgb = PLS::intensities.to_rgb();
	fprintf(f_log, "- Intensity: [%f %f %f]\n", rgb[0], rgb[1], rgb[2]);
}

#endif //_COSINE_LIGHT_SOURCE_H

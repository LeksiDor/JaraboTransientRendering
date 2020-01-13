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

#ifndef _RECTANGULAR_AREA_LIGHT_SOURCE_H_
#define _RECTANGULAR_AREA_LIGHT_SOURCE_H_

#include <Image/Image.h>
#include "LightSource/LightSource.h"
#include "LinearAlgebra/Vector3.h"
#include "RayTracing/World.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Utils/RandomNumbers.h"
#include "Sampling/Sampling.h"
#include "Material/Texture/Texture.h"

//=============================================================
// Class modeling an area light as a rectangle in 3D. This is a 
// common type of light source, and allows efficient solid-angle 
// sampling. The implementation is based on the method of Urena 
// et al.[1].
//
// [1]	C. Urena, M. Fajardo, A. King 2013. An Area Preserving 
//		Parametrization for Spherical Rectangles. In Computer 
//		Graphics Forum, Vol. 32(4)
//		https://www.solidangle.com/research/egsr2013_spherical_rectangle.pdf
template<class Radiance>
class RectangularAreaLightSource : public LightSource<3, Radiance>
{
protected:
	typedef LightSource<3, Radiance> LS;

public:
	enum SamplingMethod { UNIFORM, SPHERICAL_RECTANGLE };

protected:
	Vector3 m_position;
	Vector3 m_normal, m_up, m_left;
	Vector2 m_size;

	Real m_area, m_inv_area;

	Texture* m_emission_map;

	SamplingMethod m_sampling;

	void setup()
	{
		m_inv_area = 1./m_area;

		if( fabs(dot(m_normal, Vector3(0,1,0))) < (0.9999) )
			m_up = Vector3(0,1,0);
		else
			m_up = Vector3(0,0,1);

		m_left = cross( m_normal, m_up );
		m_up = cross( m_normal, m_left );

		m_left *= m_size[0];
		m_up *= m_size[1];
	}

	// ------------------------------------------------------------
	// Surface sampling functions
	void sampling_uniform( Vector2 &uv, Real &pdf ) const;
	void sampling_spherical_rectangle( const Vector3 &target, Vector2 &uv, Real &pdf ) const;
	inline void sample_position( const Vector3 &target, Vector2 &uv, Real &pdf ) const;
public:
	RectangularAreaLightSource(World<3, Radiance>* _world, const Vector3& pos, const Vector3 &normal,
		const Vector2& size, const Radiance& ints, Texture* map = nullptr, Real time = 0.,
		SamplingMethod sampling = SamplingMethod::UNIFORM) :
			LightSource<3, Radiance>(_world, ints, time), m_position(pos), m_normal(normal), m_size(size), m_area(size[0] * size[1]),
			m_emission_map(map), m_sampling(sampling)
	{
		setup();
	}

	virtual ~RectangularAreaLightSource()
	{
		delete m_emission_map;
	}

	// ------------------------------------------------------------
	// Deprecated functions. 
	//	Not useful for area lights since assume that the light source 
	//	is a singularity in space. 
	Vector3 get_position() const { return m_position; }
	Vector3 get_incoming_direction(const Vector3 &point_lighted) const
	{
		return normalize(point_lighted-m_position);
	}

	bool is_visible(const Vector3 &point_lighted) const
	{
		fprintf(stderr, "Warning - RectangularAreaLightSource::is_visible. Note that this function works as a point-to-point visibility function\n");
		Vector3 dir = point_lighted-m_position;
		Real l = dir.length(); dir /= l;
		Ray<3> ray(m_position, dir, false);
		return( !this->world->intersects(ray, l, _SIGMA_VISIBILITY_)) ;
		
				
		
		Intersection<3> it;
		this->world->first_intersection(ray,it,l+100);

		//printf("Distance: %f vs %f [%d]", l, ray.get_parameter(), ray.did_hit());
		if( ray.did_hit() && ray.get_parameter() < l )
			return false;

		return true;
	}

	virtual Radiance get_incoming_light(const Vector3 &point_lighted) const
	{
		fprintf(stderr, "Warning - RectangularAreaLightSource::is_visible. Note that this function works as a point-to-point visibility function\n");
		Real distance2 = (point_lighted-m_position).length2();
		return this->intensities*(1.f/distance2);
	}

	// ------------------------------------------------------------
	// Sampling functions
	virtual bool sample(const Intersection<3> &it, LightSample<3, Radiance> &light_sample, Real &pdf)const;

	virtual void sample(LightSample<3, Radiance> &light_sample, Real &pdf)const;
	
	virtual bool sample(const Vector3 &p, LightSample<3, Radiance> &light_sample, Real &pdf)const;
}; //RectangularAreaLightSource

//=============================================================================
// Naively samples the area light by uniformly sampling its surface.
template<class Radiance>
void RectangularAreaLightSource<Radiance>::sampling_uniform(Vector2 &uv, Real &pdf) const
{
	uv[0] = Random::StdRNG.next_real();
	uv[1] = Random::StdRNG.next_real();

	pdf = m_inv_area;
}

//=============================================================================
// Samples the surface of the rectangle based on its projection on the sphere 
// over the point target [1].
//
// NOTE: Not implemented yet, right now it only does the naive approach.
template<class Radiance>
void RectangularAreaLightSource<Radiance>::sampling_spherical_rectangle(const Vector3 &target, Vector2 &uv, Real &pdf) const
{
	uv[0] = Random::StdRNG.next_real()-.5;
	uv[1] = Random::StdRNG.next_real()-.5;

	pdf = m_inv_area;
}


//=============================================================================
// Samples the rectangle, based on the method determined by m_sampling.
template<class Radiance>
void RectangularAreaLightSource<Radiance>::sample_position(const Vector3 &target, Vector2 &uv, Real &pdf) const
{
	switch( m_sampling )
	{
	case UNIFORM:
		sampling_uniform( uv, pdf );
		break;
	case SPHERICAL_RECTANGLE:
		sampling_spherical_rectangle( target, uv, pdf );
		break;
	}
}

//=============================================================================
template<class Radiance>
bool RectangularAreaLightSource<Radiance>::sample(const Intersection<3> &it, LightSample<3, Radiance> &light_sample, Real &pdf)const
{
	Vector2 uv;
	sample_position(it.get_position(), uv, pdf );

	light_sample.pos = m_position + m_left*uv[0] + m_up*uv[1];
	light_sample.dir = it.get_position()-light_sample.pos;
	light_sample.dist = light_sample.dir.length();
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = LS::time;
	
	Real dot_light = dot(light_sample.dir, m_normal);

	Ray3D shadow(light_sample.pos, light_sample.dir);
	if (dot_light > 0. && it.get_normal().dot(-light_sample.dir) > 0. && !LS::world->intersects(shadow, light_sample.dist)) {
		light_sample.irradiance = LS::intensities*dot_light / (light_sample.dist*light_sample.dist);

		if( m_emission_map )
		{
			light_sample.irradiance *= (*m_emission_map)(uv[0], uv[1]);
		}
		return true;
	}
	
	return false;
}

//=============================================================================
template<class Radiance>
bool RectangularAreaLightSource<Radiance>::sample(const VectorN<3> &p, LightSample<3, Radiance> &light_sample, Real &pdf)const
{
	Vector2 uv;
	sample_position( p, uv, pdf );

	light_sample.pos = m_position + m_left*uv[0] + m_up*uv[1];
	light_sample.dir = p-light_sample.pos;
	light_sample.dist = light_sample.dir.length();
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = LS::time;
	
	Real dot_light = dot(light_sample.dir, m_normal);

	Ray3D shadow(light_sample.pos, light_sample.dir);
	if (dot_light > 0. && !LS::world->intersects(shadow, light_sample.dist) ){
		light_sample.irradiance = LS::intensities*dot_light / (light_sample.dist*light_sample.dist);

		if( m_emission_map )
		{
			light_sample.irradiance *= (*m_emission_map)(uv[0], uv[1]);
		}
		return true;
	}
	
	return false;
}

//=============================================================================
template<class Radiance>
void RectangularAreaLightSource<Radiance>::sample(LightSample<3, Radiance> &light_sample, Real &pdf)const
{
	Vector2 uv;
	sampling_uniform( uv, pdf );

	light_sample.pos = m_position + m_left*uv[0] + m_up*uv[1];
	light_sample.dist = 0.;
	light_sample.instant = LS::time;
	
	if( m_emission_map )
	{
		light_sample.irradiance *= (*m_emission_map)(uv[0], uv[1]);
	}
	
	Vector3 outgoing_direction;
	Real pdf_direction;
	Sampling.cosine_weighted(outgoing_direction, pdf_direction);

	light_sample.dir = outgoing_direction.transform_matrix_to(Vector3(0,1,0), m_normal);

	Real dot_light = dot(light_sample.dir, m_normal);
	light_sample.irradiance = LS::intensities*dot_light;

	pdf *= pdf_direction;
}

#endif //_RECTANGULAR_AREA_LIGHT_SOURCE_H_

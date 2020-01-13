/* 
 * Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _LAMBERTIAN_TEXTURED_H_
#define _LAMBERTIAN_TEXTURED_H_

#include "bunnykiller.h"

#include <vector>

#include "Image/Image.h"
#include "Material/Material.h"
#include "Material/Texture/Texture.h"
#include "Material/Reflectance/Lambertian.h"

/**	Class defining the simplest emissive
	material. 
	TODO: Parametrize the light distribution,
		now it's only gaussian. */
template<unsigned D>
class LambertianTextured: public Material<D>
{
	Spectrum m_kd;
	Texture* m_texture;
private:
	inline Spectrum get_color(const Vector2& uv) const
	{
		return (m_kd * (*m_texture)(uv));
	}
public:
	LambertianTextured() :
		Material<D>(Reflectance::DIFFUSE),
		m_kd(1.),
		m_texture(nullptr)
	{}

	LambertianTextured(const char* name) :
		Material<D>(Reflectance::DIFFUSE),
		m_kd(1.),
		m_texture(new Texture(Imaging::load<Real>(name)))
	{}

	LambertianTextured(Spectrum kd, const char* name) :
		Material<D>(Reflectance::DIFFUSE),
		m_kd(kd),
		m_texture(new Texture(Imaging::load<Real>(name)))
	{}

	virtual ~LambertianTextured()
	{
		if (m_texture) {
			delete m_texture;
		}
	}
	
	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, Spectrum &R) const
	{
		Lambertian<D> reflectance(get_color(uv)); 

		reflectance.f(omega_i, omega_o, normal, uv, R);
	}

	Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv)const
	{
		Lambertian<D> reflectance(1.);

		return reflectance.p(omega_i,omega_o,normal,uv);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		Lambertian<D> reflectance(get_color(uv)); 

		reflectance.sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
			Real &pdf) const
	{
		Lambertian<D> reflectance(get_color(it.get_uv())); 

		reflectance.sample_outgoing_ray(it, new_ray, R, pdf);
	}

	Spectrum get_absorption(const Intersection<D> &it) const
	{
		return Spectrum(1.) - get_color(it.get_uv());
	}
}; //LambertianTextured

typedef LambertianTextured<3> LambertianTextured3D;
typedef LambertianTextured<2> LambertianTextured2D;

#endif //_LAMBERTIAN_TEXTURED_H_

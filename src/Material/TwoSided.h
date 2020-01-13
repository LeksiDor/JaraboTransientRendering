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

#ifndef _TWO_SIDED_MATERIAL_H_
#define _TWO_SIDED_MATERIAL_H_

#include "bunnykiller.h"

#include <exception>
#include <vector>

#include "Material/Material.h"

template<unsigned int D>
class TwoSidedMaterial : public Material<D>
{
	const Material<D>* m_base_material;
public:
	TwoSidedMaterial(const Material<D> *base_material) :
			Material<D>(Reflectance::merge(base_material->get_type(), Reflectance::TWO_SIDED)),
			m_base_material(base_material)
	{
		if (m_base_material->is_type(Reflectance::TRANSMISSION))
			throw std::runtime_error("Error - TwoSided material with transmissive base material!");
	}

	~TwoSidedMaterial()
	{
		if (m_base_material)
			delete m_base_material;
	}

	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Spectrum &R) const
	{
		if (dot(-omega_i, normal) < 0.)
			m_base_material->f(omega_i, omega_o, -normal, uv, R);
		else
			m_base_material->f(omega_i, omega_o, normal, uv, R);
	}

	Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv) const
	{
		if (dot(-omega_i, normal) < 0.)
			return m_base_material->p(omega_i, omega_o, -normal, uv);

		return m_base_material->p(omega_i, omega_o, normal, uv);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		if (dot(-omega_i, normal) < 0.)
			m_base_material->sample_direction(omega_i, omega_o, -normal, uv, R, pdf);
		else
			m_base_material->sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
		Real &pdf) const
	{
		Intersection<D> it1 = it;

		if (dot(-it.get_ray().get_direction(), it.get_normal()) < 0.)
			it1.invert_coordinate_system();

		m_base_material->sample_outgoing_ray(it1, new_ray, R, pdf);
	}

	Spectrum get_absorption(const Intersection<D> &it) const
	{
		return m_base_material->get_absorption(it);
	}
}; // TwoSidedMaterial

typedef TwoSidedMaterial<3> TwoSidedMaterial3D;
typedef TwoSidedMaterial<2> TwoSidedMaterial2D;

#endif //_TWO_SIDED_MATERIAL_H_

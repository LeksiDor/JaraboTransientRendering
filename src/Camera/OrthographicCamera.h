/*
 * Copyright (C) 2017, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef __ORTHOGRAPHIC_H
#define __ORTHOGRAPHIC_H

#include "bunnykiller.h"

#include <cstdio>

#include "LinearAlgebra/Vector2.h"
#include "LinearAlgebra/Vector3.h"
#include "Camera/Camera.h"
#include "RayTracing/Ray.h"
#include "Geometry/Object.h"

/** This class represents an orthographic camera  */
class OrthographicCamera: public Camera<3>
{
protected:
	/// Position of camera
	Vector3 m_position;
	/// Direction of m_up
	Vector3 m_up;
	/// Basis of camera coordinate system
	Vector3 m_nv, m_uv, m_vv;
	/// Dimensions of the camera window
	Real m_w, m_h;
public:
	/** Build the camera from viewing information.
	 */
	OrthographicCamera(const Vector3& _pos, const Vector3& _focus, const Vector3& _up,
			const Real _w=1.0, const Real _h=1.0) :
		m_position(_pos), m_w(_w), m_h(_h)
	{
		m_up = _up;
		m_up.normalize();

		Vector3 dir = (_focus - m_position).normalized();

		// Calculate ref system
		m_nv = dir;
		m_uv = cross(m_nv, m_up).normalized();
		m_vv = cross(m_uv, m_nv).normalized();
	}

	virtual ~OrthographicCamera() {}

	/// Calculates the ray
	virtual const Vector3 get_ray_direction(const Vector2& p) const
	{
		return m_nv;
	}
	
	virtual Ray3D get_ray(const Vector2& p) const
	{
		Vector3 pos = m_position + p[0]*m_w*m_uv - p[1]*m_h*m_vv;

		return Ray3D(pos, m_nv);
	}

	virtual Ray3D get_ray(const Vector2& p, Real &pdf, Real &time_delay) const
	{
		pdf = 1.; time_delay = 0.;

		Vector3 pos = m_position + p[0]*m_w*m_uv + p[1]*m_h*m_vv;

		return Ray3D(pos, m_nv);
	}

	virtual const PolarizationFrame<3> get_frame() const override
	{
		return PolarizationFrame<3>(m_up, m_nv);
	}

	virtual void print(FILE *stream) const
	{
		fprintf(stream, "Orthographic Camera:\n");
		fprintf(stream, "- Position: [%f, %f, %f]\n",  m_position[0], m_position[1], m_position[2]);
		fprintf(stream, "- Direction: [%f, %f, %f]\n", m_nv[0], m_nv[1], m_nv[2]);
		fprintf(stream, "- View plane: [%f %f]\n", m_w, m_h);
		fprintf(stream, "- U: [%f, %f, %f]\n", m_uv[0], m_uv[1], m_uv[2]);
		fprintf(stream, "- V: [%f, %f, %f]\n", m_vv[0], m_vv[1], m_vv[2]);
	}
};

#endif

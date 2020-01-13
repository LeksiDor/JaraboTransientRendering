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

#ifndef _PINHOLE_CAMERA_H_
#define _PINHOLE_CAMERA_H_

#include "bunnykiller.h"

#include <cstdio>

#include "Camera/Camera.h"
#include "LinearAlgebra/Vector3.h"
#include "LinearAlgebra/Vector2.h"
#include "RayTracing/Ray.h"

/** This class represents a pinhole camera.*/
class PinholeCamera : public Camera<3>
{
protected:
	/// Position of camera
	Vector3 m_position;
	/// Point we look at
	Vector3 m_focus;
	/// Direction of m_up
	Vector3 m_up;
	/// Direction in which we look
	Vector3 m_direction;
	/// Basis of camera coordinate system
	Vector3 m_nv, m_uv, m_vv;
	/// Distance from camera to focal plane
	Real m_focal_dist;
public:
	/** Build the camera from viewing information.
	 c is that camera position, f is the point we look at,
	 u is the m_up vector and fd is the distance to the image
	 plane. */
	PinholeCamera(const Vector3& c, const Vector3& f, const Vector3& u, const Real fd) :
			m_position(c),
			m_focus(f),
			m_up(u),
			m_focal_dist(fd)
	{
		m_up.normalize();
		
		// line of sight is calculated on the 
		// basis of the camera position and 
		// focus point
		m_direction = m_focus - m_position;

		//----------------------------------------------
		// The basis of our coordinate system (nv,uv,vv)
		//----------------------------------------------
		m_nv = m_direction.normalized();
		m_uv = cross(m_nv, m_up).normalized();
		m_vv = cross(m_uv, m_nv).normalized();
	}

	virtual ~PinholeCamera()
	{
	}

	/// Get direction of viewing ray from image coords.
	virtual const Vector3 get_ray_direction(const Vector2& p) const
	{
		Vector3 dir(m_nv);
		dir *= m_focal_dist;
		dir += 2. * m_uv * p[0];
		dir -= 2. * m_vv * p[1];
		dir.normalize();

		return dir;
	}

	/// Returns position of camera.
	virtual const Vector3& get_position() const
	{
		return m_position;
	}

	Vector3& get_position()
	{
		return m_position;
	}
	
	/// Calculates the ray
	virtual Ray3D get_ray(const Vector2& p) const
	{
		return Ray3D(m_position, get_ray_direction(p));
	}

	virtual Ray3D get_ray(const Vector2& p, Real &pdf, Real &time_delay) const
	{
		pdf = 1.;
		time_delay = 0.;
		return Ray3D(m_position, get_ray_direction(p));
	}

	virtual const PolarizationFrame<3> get_frame() const override
	{
		return PolarizationFrame<3>(m_up, m_nv);
	}

	virtual void print(FILE *stream) const
	{
		fprintf(stream, "Pinhole Camera:\n");
		fprintf(stream, "- Position: [%f, %f, %f]\n", m_position[0], m_position[1], m_position[2]);
		fprintf(stream, "- Direction: [%f, %f, %f]\n", m_nv[0], m_nv[1], m_nv[2]);
		fprintf(stream, "- Focal Distance: %f\n", m_focal_dist);
		fprintf(stream, "- U: [%f, %f, %f]\n", m_uv[0], m_uv[1], m_uv[2]);
		fprintf(stream, "- V: [%f, %f, %f]\n", m_vv[0], m_vv[1], m_vv[2]);
	}
};
// PinholeCamera

#endif // _PINHOLE_CAMERA_H_

/*
 * Copyright (C) 2018, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _WARD_H_
#define _WARD_H_

#include "bunnykiller.h"

#include <cmath>

#include "Material/Reflectance/BSDF.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Utils/RandomNumbers.h"

/**
 * Ward BRDF
 * Based on the corrected version with bounded albedo [1].
 *
 * [1] D. Geisler-Moroder, A. Dür, 2010
 *     A New Ward BRDF Model with Bounded Albedo
 *     Computer Graphics Forum 29(4). Blackwell Publishing Ltd.
 *
 * TODO: Develop 2D version of the BRDF, right now only 3D.
 */
class Ward : public BSDF<3>
{
	/* Specular coefficients */
	Spectrum m_ps;
	/* Anisotropic Gloss Lobe */
	Real m_ax, m_ay;
public:
	Ward(const Spectrum& ps = Spectrum(1.0)) :
			BSDF<3>(Reflectance::GLOSSY),
			m_ps(ps),
			m_ax(0.),
			m_ay(0.)
	{
	}

	Ward(const Spectrum& ps, const Real a) :
			BSDF<3>(Reflectance::GLOSSY),
			m_ps(ps),
			m_ax(a),
			m_ay(a)

	{
	}

	Ward(const Spectrum& ps, const Real ax, const Real ay) :
			BSDF<3>(Reflectance::GLOSSY),
			m_ps(ps),
			m_ax(ax),
			m_ay(ay)

	{
	}
	
	void f(const Vector3 &omega_i, const Vector3 &omega_o, const Vector3 &normal, const Vector2 &uv,
		Spectrum &R) const override
	{
		const Vector3& i = -omega_i;
		const Vector3& o = omega_o;
		const Vector3& n = normal;

		if (dot(n, o) < 0. || dot(n, i) < 0.) {
			R = Spectrum(0.);
			return;
		}

		/* Build an orthonormal base by Hughes-Moeller method */
		const Vector3 tx = normalize(
				(std::abs(n[0]) > std::abs(n[2])) ?
						Vector3(-n[1], n[0], n[2]) : Vector3(0.0, -n[2], n[1]));
		const Vector3 ty = cross(n, tx);

		/* Calculate half vector */
		Vector3 h = (i + o);

		Real exponent = -(std::pow(dot(h, tx) / m_ax, 2.f) + std::pow(dot(h, ty) / m_ay, 2.f))
				/ std::pow(dot(h, n), 2.f);

		Real power = (std::exp(exponent) * dot(h, h))
				/ (M_PI * m_ax * m_ay * std::pow(dot(h, n), 4.f));

		R = m_ps * power;
	}

	Real p(const Vector3 &omega_i, const Vector3 &omega_o, const Vector3 &normal,
		const Vector2 &uv) const override
	{
		const Vector3& i = -omega_i;
		const Vector3& o = omega_o;
		const Vector3& n = normal;

		/* Build an orthonormal base by Hughes-Moeller method */
		Vector3 tx = normalize(
				(std::abs(n[0]) > std::abs(n[2])) ?
						Vector3(-n[1], n[0], n[2]) : Vector3(0.0, -n[2], n[1]));
		Vector3 ty = cross(n, tx);

		/* Calculate half vector */
		Vector3 h = normalize(i + o);

		Real cos_theta_h = dot(h, n);
		Real sin_theta_h = std::sqrt(1.f - cos_theta_h * cos_theta_h);
		Real tan_theta_h = sin_theta_h / cos_theta_h;

		Real cos_phi_h = dot(h, tx) / sin_theta_h;
		Real sin_phi_h = dot(h, ty) / sin_theta_h;

		Real cos_sin_phi_factor = (cos_phi_h * cos_phi_h) / (m_ax * m_ax)
				+ (sin_phi_h * sin_phi_h) / (m_ay * m_ay);

		/* Calculate pdf and BRDF value */
		return std::exp(-tan_theta_h * tan_theta_h * cos_sin_phi_factor)
				/ (4.f * M_PI * m_ax * m_ay * dot(i, h) * std::pow(cos_theta_h, 3.f));
	}

	void sample_direction(const Vector3 &omega_i, Vector3 &omega_o, const Vector3 &normal,
		const Vector2 &uv, Spectrum &R, Real &pdf) const override
	{
		const Vector3& i = -omega_i;
		const Vector3& n = normal;

		Real e1 = Random::StdRNG.next_real(), e2 = Random::StdRNG.next_real();

		/* Sample zenith angle */
		Real phi_h = std::atan2(m_ay * std::sin(2.f * M_PI * e2), m_ax * std::cos(2.f * M_PI * e2));

		Real cos_phi_h = std::cos(phi_h);
		Real sin_phi_h = std::sin(phi_h);

		Real cos_sin_phi_factor = (cos_phi_h * cos_phi_h) / (m_ax * m_ax)
				+ (sin_phi_h * sin_phi_h) / (m_ay * m_ay);

		/* Sample azimuth angle */
		Real theta_h = std::atan2(std::sqrt(-std::log(1 - e1)), std::sqrt(cos_sin_phi_factor));

		Real cos_theta_h = std::cos(theta_h);
		Real sin_theta_h = std::sin(theta_h);
		Real tan_theta_h = sin_theta_h / cos_theta_h;

		/* Build an orthonormal base by Hughes-Moeller method */
		const Vector3 tx = normalize(
				(std::abs(n[0]) > std::abs(n[2])) ?
						Vector3(-n[1], n[0], n[2]) : Vector3(0.0, -n[2], n[1]));
		const Vector3 ty = cross(n, tx);

		/* Calculate half vector */
		Vector3 h(sin_theta_h * cos_phi_h, sin_theta_h * sin_phi_h, cos_theta_h);

		/* Transform sampled half vector to normal reference frame */
		h = h[0] * tx + h[1] * ty + h[2] * n;

		/* Calculate sampled direction */
		Vector3 o(2.f * dot(i, h) * h - i);

		omega_o = o;

		/* Calculate pdf and BRDF value */
		if (dot(n, o) < 0. || dot(n, i) < 0.) {
			pdf = 0.0;
			R = Spectrum(0.);
		} else {
			pdf = std::exp(-tan_theta_h * tan_theta_h * cos_sin_phi_factor)
					/ (4.f * M_PI * m_ax * m_ay * dot(i, h) * std::pow(cos_theta_h, 3.f));

			f(omega_i, omega_o, normal, uv, R);
		}
	}

	void sample_outgoing_ray(const Intersection<3> &it, Ray<3> &new_ray, Spectrum &R,
		Real &pdf) const override
	{
		const Ray<3>& old_ray = it.get_ray();

		Vector3 new_omega_o;
		sample_direction(old_ray.get_direction(), new_omega_o, it.get_normal(), it.get_uv(), R,
				pdf);

		new_ray = Ray<3>(it.get_position(), new_omega_o, true, old_ray.get_level() + 1,
				old_ray.get_ior(), old_ray.get_medium());
	}

	Spectrum get_absorption(const Vector2 &uv) const override
	{
		return Spectrum(1.) - m_ps;
	}
}; /* Ward */

#endif /* _WARD_H_ */

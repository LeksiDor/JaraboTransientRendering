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

#ifndef _PHONG_H_
#define _PHONG_H_

#include "Material/Reflectance/BSDF.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Sampling/Sampling.h"
#include <math.h>

/**	Physically-based Blinn-Phong BRDF:
	Blinn, J. 1977. "Models of light reflection for computer synthesized 
	pictures". In SIGGRAPH 1977.
	http://dl.acm.org/citation.cfm?doid=563858.563893	
	[TODO] Develop 2D versions of the BRDF, right now only 3D. */
class Phong: public BSDF3D
{
	//Specular and Diffuse coefficients
	Spectrum Kd, Ks;
	//Specular sharpness (glossiness)
	Real Ns;
public:
	Phong() :
		BSDF<3>(Reflectance::GLOSSY), Kd(.67), Ks(1.), Ns(0.)
	{}
	Phong(const Spectrum &_Kd, const Real _Ns) :
		BSDF<3>(Reflectance::GLOSSY), Kd(_Kd), Ks(1.), Ns(_Ns)
	{}
	Phong(const Spectrum &_Kd, const Spectrum &_Ks, const Real _Ns) :
		BSDF<3>(Reflectance::GLOSSY), Kd(_Kd), Ks(_Ks), Ns(_Ns)
	{}
	
	virtual ~Phong() {}

	void f(const Vector3 &omega_i, const Vector3 &omega_o, const Vector3 &normal, const Vector2 &uv, Spectrum &R) const
	{
		Vector3 half_vector = normalize(-omega_o+omega_i);
		Real clamped_dot = dot(normal, half_vector); clamped_dot = (clamped_dot>0.f)?clamped_dot:0.f;
		Real cosine_lobe = pow(clamped_dot, Ns);

		R = (Kd + Ks*(.5f*(Ns+2.f)*cosine_lobe))*(1.f/M_PI);
	}

	Real p(const Vector3 &omega_i, const Vector3 &omega_o, const Vector3 &normal, const Vector2 &uv) const
	{
		//Sample hemisphere cosine weighted-> Need to be redone...
		return (dot(omega_o,normal)/static_cast<Real>(M_PI));
	}

	void sample_direction(const Vector3 &omega_i, Vector3 &omega_o, const Vector3 &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		//Sample hemisphere cosine weighted-> Need to be redone...
		Sampling.cosine_weighted(omega_o, pdf);
		omega_o.transform_matrix_to(VectorN<3>(0, 1), normal);

		f(-omega_i, -omega_o, normal, uv, R);
	}

	void sample_outgoing_ray(const Intersection3D &it, Ray3D &new_ray, Spectrum &R, Real &pdf) const
	{
		Vector3 new_omega_o;
		sample_direction(it.get_ray().get_direction(), new_omega_o, it.get_normal(),it.get_uv(), R, pdf);

		new_ray = Ray3D(it.get_position(), new_omega_o, true, it.get_ray().get_level()+1, 
			it.get_ray().get_ior(), it.get_ray().get_medium());
	}

	Spectrum get_absorption( const Vector2 &uv ) const
	{
		return 1-Kd;
	}
}; //Phong

#endif //_PHONG_H_

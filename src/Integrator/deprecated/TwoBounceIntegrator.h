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

#ifndef _TWO_BOUNCE_INTEGRATOR_H_
#define _TWO_BOUNCE_INTEGRATOR_H_

#include "Integrator/Integrator.h"
#include "RayTracing/World.h"
#include "Sampling/SamplingTechniques.h"
#include "RayTracing/RayTraceDirection.h"
#include "Path.h"
#include "PathMIS.h"

//=============================================================
template<class Radiance, class RadianceAttenuation>
class TwoBounceIntegrator : public Integrator<3, Radiance, RadianceAttenuation>
{

public:
	TwoBounceIntegrator(World<3, Radiance> *w, unsigned int incoming_samples = 1, FILE *_f_log = NULL, FilmR *_m_film = NULL)
		:Integrator<3, Radiance, RadianceAttenuation>(w, incoming_samples, _f_log, _m_film)
		{	}
	~TwoBounceIntegrator()
	{}

	// ------------------------------------------------------------------------------------------
	void operator()(const VectorN<3> &p, RadSampleList &samples, const Real delta_time, const unsigned int nb_time_samples) const;

	void operator()(const Ray<3> &r, RadSampleList &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Radiance operator()(const Ray<3> &r) const;// {return Spectrum(0.);};

}; //TwoBounceIntegrator


//=============================================================================
// Computes the light transport from an initial camera-sampled ray r.
template<class Radiance, class RadianceAttenuation>
Radiance TwoBounceIntegrator<Radiance, RadianceAttenuation>::operator()(const Ray<3> &r) const
{
	return Radiance(0);
}


template<class Radiance, class RadianceAttenuation>
void TwoBounceIntegrator<Radiance, RadianceAttenuation>::operator()(const Ray<3> &r, RadSampleList &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	//Spectrum reflected_radiance(0.); // FIX
	Ray<3> curr_ray = r;

	for (int sample = 0; sample < this->m_incoming_samples; ++sample)
	{
		Intersection<3> it;
		this->m_world->first_intersection(curr_ray, it);

		if (!it.did_hit())
			return;

		Real albedo = 1.;
		Radiance f;
		Real p;
		Real t;

		
		bool angular = true;
		if (angular)
		{
			VectorN<3> new_omega_o;
			//Sampling::cosine_weighted(new_omega_o, p);
			Sampling::direction_uniformly(new_omega_o, p);

			//new_omega_o = new_omega_o[0] * it.get_tangent_x() + new_omega_o[1] * it.get_normal() + new_omega_o[2] * it.get_tangent_y();

			f = albedo * M_1_PI * dot_clamped(new_omega_o, it.get_normal());
			t = curr_ray.get_parameter();
			
			//printf("Intersection point 1: %f, %f, %f\n", it.get_position()[0], it.get_position()[1], it.get_position()[2]);

			curr_ray = Ray<3>(it.get_position(), new_omega_o);
			curr_ray.shift();
			

			it = Intersection<3>();
			m_world->first_intersection(curr_ray, it);

			if (!curr_ray.did_hit())
				continue;

			//printf("Intersection point2 : %f, %f, %f\n\n", it.get_position()[0], it.get_position()[1], it.get_position()[2]);

			LightSample<3, Radiance> light_sample;
			Real pl;

			if (!m_world->light(0)->sample(it.get_position(), light_sample, pl))
				continue;

			//printf("Intersection point1 : %f, %f, %f\n\n", it.get_position()[0], it.get_position()[1], it.get_position()[2]);
			//printf("Intersection point2 : %f, %f, %f\n\n", curr_ray.get_position()[0], curr_ray.get_position()[1], curr_ray.get_position()[2]);

			f *= light_sample.irradiance * albedo * M_1_PI * dot_clamped(-light_sample.dir, it.get_normal());
			p *= pl;
			
			t += light_sample.dist + curr_ray.get_parameter() + light_sample.instant;

		}
		else
		{
			// Sample point in the occluded geometry
			Real x = RNG::StdRandom::get_real()*.66 - .33;
			Real y = RNG::StdRandom::get_real()*.66 - .33;
			Real z = 2.2;

			p = 1 / (.66*.66);

			Vector3 pi(x, y, z);
			Vector3 ni(0, 0, -1);

			//printf("Intersection point: %f, %f, %f\n", x,y,z);

			// Contribution light to point
			LightSample<3, Radiance> light_sample;
			Real pl;

			// Check light visibility
			if (!m_world->light(0)->sample(pi, light_sample, pl))
				continue;

			// Check surface visibility
			// Assumed...

			f = light_sample.irradiance * albedo * M_1_PI * dot_clamped(-light_sample.dir, ni);
			
			Vector3 w2i = pi - it.get_position();
			Real d_w2i_2 = w2i.length2();
			w2i.normalize();

			f *= albedo * M_1_PI * dot_clamped(-w2i, ni) * dot_clamped(w2i, it.get_normal()) / (d_w2i_2);

			t = m_world->time_of_flight(light_sample.dist + sqrt(d_w2i_2) + curr_ray.get_parameter())*m_world->get_ior() + light_sample.instant;
			p *= pl;
		}

		samples.push_back(RadianceSampleR(f / p*m_inv_incoming_samples, t, 0));
	}
}


template<class Radiance, class RadianceAttenuation>
void TwoBounceIntegrator<Radiance, RadianceAttenuation>::operator()(const VectorN<3> &p, RadSampleList &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
}

#endif //_BIDIRECTIONAL_PATH_TRACING_INTEGRATOR_H_
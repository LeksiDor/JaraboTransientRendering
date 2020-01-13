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

#ifndef _PHOTON_ESTIMATION_H_
#define _PHOTON_ESTIMATION_H_

#include "bunnykiller.h"

#include <vector>

#include "Integrator/Integrator.h"
#include "Integrator/BidirectionalPathTracing.h"
#include "LinearAlgebra/VectorN.h"
#include "RayTracing/World.h"

template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
class PhotonEstimation : public BidirectionalPathTracing<D, Radiance, RadianceAttenuation>
{
protected:
	using BDPTR = BidirectionalPathTracing<D, Radiance, RadianceAttenuation>;
	using RadianceSampleR = typename BDPTR::RadianceSampleR;
	using RadSampleList = typename BDPTR::RadSampleList;
	using FilmR = typename BDPTR::FilmR;
	using LightSampleR = typename BDPTR::LightSampleR;
    using PathR = typename BDPTR::PathR;
    using VertexR = typename BDPTR::VertexR;
public:
	enum BeamType
	{
		LongBeam, ShortBeam
	};
protected:
    /* Beam record */
    struct Beam
	{
		/* Ray information */
		VectorN<D> m_origin;
		VectorN<D> m_direction;
		Real m_length;
		/* Photon */
		Real m_time;
		Radiance m_gathered;
		Medium<D>* m_medium;
	};

    std::vector<Beam> m_beams;
protected:
    Real m_radius;
    unsigned int m_nb_photon_shots;

    BeamType m_beam_type;
#ifdef _USE_EMBREE_
    RTCDevice g_device = nullptr;
    RTCScene g_scene = nullptr;
public:
	/* We pass custom info using the Embree ray */
	struct BeamRay
	{
		RTCRay rtc_ray;
		Ray<D> ray;
		Radiance gathered;
	};
#endif // _USE_EMBREE_
public:
    PhotonEstimation(World<D, Radiance>& w,  Real radius, unsigned int nb_photon_shots,
    		unsigned int incoming_samples = 1, FILE *_f_log = nullptr,
			FilmR *_m_film = nullptr) :
		BDPTR(w, incoming_samples, _f_log, _m_film),
			m_radius(radius),
			m_nb_photon_shots(nb_photon_shots),
			m_beam_type(BeamType::LongBeam)
	{
#ifdef _USE_EMBREE_
		/* Create new Embree device */
		g_device = rtcNewDevice(nullptr);
		EmbreeFunctions::rtc_error_handler(nullptr, rtcDeviceGetError(g_device));

		/* Set error handler (reuse World's one) */
		rtcDeviceSetErrorFunction2(g_device, EmbreeFunctions::rtc_error_handler, nullptr);

		/* Create scene */
		g_scene = rtcDeviceNewScene(g_device,
				RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT,
				RTC_INTERSECT1 | RTC_INTERPOLATE);
#endif // _USE_EMBREE_
	}

	virtual ~PhotonEstimation()
	{}
public:
	virtual void preprocess() override;

	virtual Radiance operator()(const Ray<D>& r) const override;
protected:
	void store_beams(const std::vector<Beam>& beams) const;

	void evaluate_beam(const BeamRay& r, const Beam& beam) const;
protected:
	Radiance volumetric_radiance(const Ray<D>& r) const;
};


#ifdef _USE_EMBREE_
namespace EmbreeFunctions {
/* Error reporting function for embree2 */
template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
void rtc_beam_filter_function(void*, RTCRay& rtc_ray)
{
	/* Perform custom processing */
	using BeamRay = typename PhotonEstimation<D, Radiance, RadianceAttenuation, VREstimator>::BeamRay;
	BeamRay& beam_ray = (BeamRay&)rtc_ray;

	beam_ray.gathered += Radiance(1., 1., 1.);
//	beam_ray.gathered += Radiance(rtc_ray.u, rtc_ray.v, 1.);
//	beam_ray.gathered += Radiance(rtc_ray.geomID, rtc_ray.primID, 1.);

	/* Reject all the beams intersections */
	rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
	rtc_ray.tfar = (float)beam_ray.ray.get_parameter();
}
} // namespace EmbreeFunctions
#endif // _USE_EMBREE_

template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
void PhotonEstimation<D, Radiance, RadianceAttenuation, VREstimator>::store_beams(
		const std::vector<Beam>& beams) const
{
#ifdef _USE_EMBREE_
	for (size_t i = 0; i < beams.size(); i++) {
		const Beam& beam = beams[i];

		unsigned int rtc_lines = rtcNewLineSegments2(g_scene, RTC_GEOMETRY_STATIC, 1, 2, 1);

		struct rtc_Vertex { float x, y, z, r; };
		rtc_Vertex* rtc_vertices = (rtc_Vertex*)rtcMapBuffer(g_scene, rtc_lines, RTC_VERTEX_BUFFER);

		VectorN<D> start = beam.m_origin;

		rtc_vertices[0].x = start[0];
		rtc_vertices[0].y = start[1];
		rtc_vertices[0].z = start[2];
		rtc_vertices[0].r = m_radius;

		VectorN<D> end = beam.m_origin + beam.m_length*beam.m_direction;

		rtc_vertices[1].x = end[0];
		rtc_vertices[1].y = end[1];
		rtc_vertices[1].z = end[2];
		rtc_vertices[1].r = m_radius;

		rtcUnmapBuffer(g_scene, rtc_lines, RTC_VERTEX_BUFFER);

		int* segments = (int*) rtcMapBuffer(g_scene, rtc_lines, RTC_INDEX_BUFFER);

		segments[0] = 0;

		rtcUnmapBuffer(g_scene, rtc_lines, RTC_INDEX_BUFFER);

		rtcSetIntersectionFilterFunction(g_scene, rtc_lines,
				(RTCFilterFunc)&EmbreeFunctions::rtc_beam_filter_function<D, Radiance, RadianceAttenuation, VREstimator>);
		rtcSetUserData(g_scene, rtc_lines, (void*)&beams[i]);
	}

	rtcCommit(g_scene);
#endif // _USE_EMBREE_
}

template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
void PhotonEstimation<D, Radiance, RadianceAttenuation, VREstimator>::evaluate_beam(
		const BeamRay& r, const Beam& beam) const
{

}

template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
Radiance PhotonEstimation<D, Radiance, RadianceAttenuation, VREstimator>::volumetric_radiance(
		const Ray<D>& r) const
{
#ifdef _USE_EMBREE_
	BeamRay beam_ray;

	/* Initialize Embree ray */
	RTCRay& rtc_ray = (RTCRay&)beam_ray;
	rtc_ray.org[0] = (float)r.get_origin()[0];
	rtc_ray.org[1] = (float)r.get_origin()[1];
	rtc_ray.org[2] = (float)r.get_origin()[2];
	rtc_ray.dir[0] = (float)r.get_direction()[0];
	rtc_ray.dir[1] = (float)r.get_direction()[1];
	rtc_ray.dir[2] = (float)r.get_direction()[2];
	rtc_ray.tnear = 0.0f;
	rtc_ray.tfar = (float)r.get_parameter();
	rtc_ray.time = 0.0f;
	rtc_ray.mask = 0XFFFFFFFF;
	rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
	rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;

	/* Initialize Bunny data */
	beam_ray.ray = r;
	beam_ray.gathered = Radiance(0.);

	rtcIntersect(g_scene, rtc_ray);

//	if (rtc_ray.geomID != RTC_INVALID_GEOMETRY_ID){
//		return Radiance(1., 0., 0.);
//	}

	return beam_ray.gathered;
#else
	return Radiance(0.);
#endif // _USE_EMBREE_
}
template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
void PhotonEstimation<D, Radiance, RadianceAttenuation, VREstimator>::preprocess()
{
#ifdef _USE_EMBREE_
	Real inv_nb_photon_shots = 1. / (Real)m_nb_photon_shots;
	const Real max_beam_length = (BDPTR::m_termination_sampling_technique == BDPTR::TerminationSampling::TotalTime) ?
		BDPTR::ITR::m_film->get_time_length() :
		BDPTR::ITR::m_world.get_bounding_volume().get_diagonal()*1.5;

	m_beams.clear();

	for (size_t i = 0; i < m_nb_photon_shots; i++) {
		/* Start in a light source */
		Real pl1, pl2;
		LightSample<D, Radiance> light_sample;
		BDPTR::m_world.sample_light(pl1)->sample(light_sample, pl2);

		const Radiance f0 = light_sample.irradiance * inv_nb_photon_shots;
		const Real p0 = pl1 * pl2;
		const Real t0 = light_sample.instant;

		// Prepare variables for the path-tracing loop
		Ray<D> curr_ray(light_sample.pos, light_sample.dir, BDPTR::ITR::m_world.get_ior(),
				BDPTR::ITR::m_world.get_medium());

		Intersection<D> curr_it;
		VectorN<D> position;

		Spectrum u_t, u_s;
		Radiance scat_albedo;
		Medium<D>* previous_medium = nullptr;
		Real avg_albedo;

		// Sampled mean free path in free-space
		Real sampled_distance = std::numeric_limits<Real>::infinity();

		// Throughput and probability of the path
		RadianceAttenuation f(1.);
		set_coordinate_system(f0, f);
		Real p = 1., t = t0;

		// Variables used through the loop
		RadianceAttenuation fa(1.);
		Real pa = 1., ta = 0.;

		RadianceAttenuation ft(1.);

		bool first_iteration = true;

		// Iterate path
		unsigned iter = 0;
		while (true) {
			iter++;

			u_t = Spectrum(0.);
			u_s = Spectrum(0.);

			// -----------------------------------------------------------------------------------------
			// With the sampled distance, sample the new ray
			fa = RadianceAttenuation(1.);
			pa = 1., ta = 0.;
			previous_medium = curr_ray.get_medium();

			if (first_iteration) {
				first_iteration = false;
			} else {
				BDPTR::sample_outgoing_ray(curr_ray, curr_it, false, fa, ta, pa);
			}

			// -----------------------------------------------------------------------------------------
			// Trace next ray
			curr_it = Intersection<D>();
			BDPTR::ITR::m_world.first_intersection(curr_ray, curr_it);

			// -----------------------------------------------------------------------------------------
			// If we are currently in a medium, get scattering parameters and calculate attenuation
			ft = RadianceAttenuation(1.);

			if (curr_ray.get_medium()) {
				// If within media...
				u_t = curr_ray.get_medium()->get_extinction(curr_ray.get_origin());
				u_s = curr_ray.get_medium()->get_scattering(curr_ray.get_origin());

				ft = curr_ray.get_medium()->get_transmittance(curr_ray);

				// -----------------------------------------------------------------------------------------
				// Generate a beam along the ray
				Beam beam;
				beam.m_origin = curr_ray.get_origin();
				beam.m_direction = curr_ray.get_direction();
				beam.m_length = std::min(curr_ray.get_parameter(), max_beam_length) - m_radius;
				beam.m_time = t;
				beam.m_gathered = f*f0;
				beam.m_medium = curr_ray.get_medium();

				m_beams.push_back(beam);
			}

			f = (f * fa * ft) / pa;
			p *= pa;

			// Test if there's a scattering event in media or surfaces. If not then the
			// path is finished...
			if (!curr_ray.get_medium() && !curr_it.did_hit())
				break;

			break;
		}
	}

	printf("BEAMS: %d\n", (int)m_beams.size());

	store_beams(m_beams);
#endif // _USE_EMBREE_
}

template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
Radiance PhotonEstimation<D, Radiance, RadianceAttenuation, VREstimator>::operator()(
		const Ray<D> &r) const
{
	Radiance reflected_radiance(0.);

	Ray<D> curr_ray(r.get_origin(), r.get_direction(), false, r.get_level(),
				r.get_refraction_index(), r.get_medium());
	Intersection<D> curr_it;
	VectorN<D> position;

	// Throughput and probability of the path
	RadianceAttenuation f(1.);
	Real p = 1., pi = 1., t = 0.;

	// Variables used through the loop
	RadianceAttenuation fa(1.);
	Real pa = 1., ta = 0.;

	RadianceAttenuation ft(1.);

	Radiance scat_albedo;
	Real avg_albedo;

	bool first_iteration = true;

	// Iterate path
	unsigned iter = 0;
	while (true) {
		iter++;

		// -----------------------------------------------------------------------------------------
		// With the sampled distance, sample the new ray
		fa = RadianceAttenuation(1.);
		pa = 1., ta = 0.;
		previous_medium = curr_ray.get_medium();

		if (first_iteration) {
			first_iteration = false;
		} else {
			BDPTR::sample_outgoing_ray(curr_ray, curr_it, false, fa, ta, pa);
		}

		// -----------------------------------------------------------------------------------------
		// Trace next ray
		curr_it = Intersection<D>();
		BDPTR::ITR::m_world.first_intersection(curr_ray, curr_it);

		// -----------------------------------------------------------------------------------------
		// If we are currently in a medium, get scattering parameters, calculate
		// attenuation and add volumetric radiance...
		ft = RadianceAttenuation(1.);

		if (curr_ray.get_medium()) {
			ft = curr_ray.get_medium()->get_transmittance(curr_ray);

//			reflected_radiance += f * volumetric_radiance(curr_ray);
		}

		// Test if there's a scattering event in surfaces. If not then the path is finished...
		if (!curr_it.did_hit())
			break;
		// -----------------------------------------------------------------------------------------
		// Update path throughput
		p = (pa * pi);
		f *= (fa * ft);
		t += ta + BDPTR::ITR::m_world.time_of_flight(curr_ray.get_parameter())*curr_ray.get_refraction_index();

		// -----------------------------------------------------------------------------------------
		// Connect with light and calculate albedo
		LightSample<D, Radiance> light_sample;
		Real pl1 = 1., pl2 = 1.;

		position = curr_it.get_position();

		if (BDPTR::ITR::m_world.sample_light(pl1)->sample(position, light_sample, pl2)) {
			RadianceAttenuation fe(0.);
			Real pe = 1., te = 0.;

			curr_it.material()->f(curr_ray.get_direction(), -light_sample.dir, curr_it.get_normal(), curr_it.get_uv(), fe, te, pe);

			/* Check surface normal */
			if (curr_it.material()->is_type(Reflectance::TWO_SIDED)) {
				fe *= dot_abs(curr_it.get_normal(), -light_sample.dir);
			} else {
				fe *= dot_clamped(curr_it.get_normal(), -light_sample.dir);
			}

			t = BDPTR::ITR::m_world.time_of_flight(light_sample.dist) * BDPTR::ITR::m_world.get_ior() + te + light_sample.instant;

			/* Compute attenuation if needed */
			if (curr_ray.get_medium()) {
				fe *= curr_ray.get_medium()->get_transmittance(
						Ray<D>(light_sample.dist, light_sample.pos, light_sample.dir));
			}
			set_attenuation_tracing_direction(fe, TraceDirection::FROM_EYE);

			reflected_radiance += f * (fe * light_sample.irradiance) / (p * pl1 * pl2 * pe);
		}

		scat_albedo = (1. - curr_it.material()->get_absorption(curr_it));

		// -----------------------------------------------------------------------------------------
		// Finally, evaluate Sample Termination!
		avg_albedo = scat_albedo.avg();

		if (BDPTR::sample_termination(avg_albedo, t, curr_ray.get_level(), pi))
			break;
	}

	return reflected_radiance;
}

#endif // _PHOTON_ESTIMATION_H_

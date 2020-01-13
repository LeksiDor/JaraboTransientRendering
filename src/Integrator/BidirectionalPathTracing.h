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

#ifndef _BIDIRECTIONAL_PATH_TRACING_INTEGRATOR_H_
#define _BIDIRECTIONAL_PATH_TRACING_INTEGRATOR_H_

#include "bunnykiller.h"

#include <algorithm>

#include "Integrator/Integrator.h"
#include "RayTracing/World.h"
#include "Sampling/SamplingTechniques.h"
#include "RayTracing/RayTraceDirection.h"
#include "Integrator/Path.h"
#include "Integrator/PathMIS.h"
#include "Color/PolarizedAttenuation.h"

template<unsigned D, class Radiance, class RadianceAttenuation>
class BidirectionalPathTracing : public Integrator<D, Radiance, RadianceAttenuation>
{
protected:
	using ITR = Integrator<D, Radiance, RadianceAttenuation>;
	using RadianceSampleR = typename ITR::RadianceSampleR;
	using RadSampleList = typename ITR::RadSampleList;
	using FilmR = typename ITR::FilmR;
	using LightSampleR = LightSample<D, Radiance>;
protected:
    using PathR = Path<D, Radiance, RadianceAttenuation>;
    using VertexR = typename PathR::Vertex;
public:
	enum SamplingMode
	{
		SteadyState, Transient
	};
	enum NextEventSampling
	{
		MeanFreePath, Time
	};
	enum AngularSampling
	{
		PhaseFunction, ShortestTime
	};
	enum TerminationSampling
	{
		RussianRoulette, NbBounces, TotalTime
	};
	enum VertexConnection
	{
		ShadowConnection, TimeSampling
	};
	enum WeightingTechnique
	{
		MIS_Balance, MIS_Power, Uniform
	};
protected:
	NextEventSampling m_next_event_sampling_technique;
	AngularSampling m_angular_sampling_technique;
	TerminationSampling m_termination_sampling_technique;
	VertexConnection m_vertex_connection_technique;
	WeightingTechnique m_weighting_technique;

	Real m_beta_MIS;

	unsigned int m_avg_nb_vertices_path;
	unsigned int m_max_path_size;

	bool m_path_trace_only;
	bool m_keep_track_sample_count;

	/* Path Sampling Functions */
	Real sample_next_event_distance(const Ray<D> &r, Real &p) const;

	void correct_next_event_distance(Ray<D> &r, Real &distance, bool &media_interaction, Real &p) const;

	void sample_outgoing_ray(Ray<D> &ray, Intersection<D> &it, bool media_interaction,
			Real sampled_distance, Real current_delay, const VectorN<D> &target,
			RadianceAttenuation &f, Real &time, Real &p) const;

	void sample_outgoing_ray(Ray<D> &ray, Intersection<D> &it, bool media_interaction,
			RadianceAttenuation &f, Real &time, Real &p) const;

	bool sample_termination(Real albedo, Real time, unsigned int nb_bounces, Real &p) const;

	void generate_path(const Ray<D> &r, const Radiance &f0, const Real p0, const Real t0,
			const VectorN<D> &target, const TraceDirection dir, PathR &path) const;

	/* Path Connecting Functions */
	bool connect_vertex_with_light(const VertexR& eye_vertex, Radiance &f, Real &p,
			Real *time = nullptr) const;

	bool connect_vertices(const VertexR& eye_vertex, const VertexR& light_vertex,
			RadianceAttenuation &f, Real *time = nullptr, Real *pdf_time = nullptr) const;

	Real compute_weights(const PathR &eye_path, const PathR &light_path, const unsigned int ie,
			const unsigned int il) const;

	Radiance connect_paths(const PathR &eye_path, const PathR &light_path) const;

	void connect_paths(const PathR &eye_path, const PathR &light_path, RadSampleList &samples,
			const Real delta_time, const unsigned int nb_time_samples) const;
public:
	BidirectionalPathTracing(World<D, Radiance>& w, unsigned int incoming_samples = 1,
			FILE *_f_log = nullptr, FilmR *_m_film = nullptr) :
		ITR(w, incoming_samples, _f_log, _m_film),
			m_next_event_sampling_technique(NextEventSampling::MeanFreePath),
			m_angular_sampling_technique(AngularSampling::PhaseFunction),
			m_termination_sampling_technique(TerminationSampling::RussianRoulette),
			m_vertex_connection_technique(VertexConnection::ShadowConnection),
			m_weighting_technique(WeightingTechnique::MIS_Power),
			m_beta_MIS(2.),
			m_avg_nb_vertices_path(4),
			m_max_path_size(30),
			m_path_trace_only(true),
			m_keep_track_sample_count(false)
	{
		set_mode(SamplingMode::SteadyState);
	}

	virtual ~BidirectionalPathTracing()
	{}

	// ------------------------------------------------------------------------------------------
	void preprocess()
	{
		if (ITR::f_log) {
			fprintf(ITR::f_log, "Bidirectional Path Tracing\n");
			fprintf(ITR::f_log, "- Samples: %d\n", ITR::m_incoming_samples);
			fprintf(ITR::f_log, "- Next Event: ");
			switch (m_next_event_sampling_technique) {
				case NextEventSampling::Time:
					fprintf(ITR::f_log, "Time Sampling - #Vertices %d\n", m_avg_nb_vertices_path);
					break;
				case NextEventSampling::MeanFreePath:
				default:
					fprintf(ITR::f_log, "Mean Free Path Sampling\n");
			}
			fprintf(ITR::f_log, "- Angular Sampling: ");
			switch (m_angular_sampling_technique) {
				case AngularSampling::PhaseFunction:
					fprintf(ITR::f_log, "Phase Function\n");
					break;
				case AngularSampling::ShortestTime:
					fprintf(ITR::f_log, "Shortest Time\n");
					break;
			}
			fprintf(ITR::f_log, "- Path Termination (Max Path Size %d): ", m_max_path_size);
			switch (m_termination_sampling_technique) {
				case TerminationSampling::TotalTime:
					fprintf(ITR::f_log, "Total Time\n");
					break;
				case TerminationSampling::RussianRoulette:
					fprintf(ITR::f_log, "Russian Roulette\n");
					break;
				default:
					break;
			}
			fprintf(ITR::f_log, "- Vertices Connection: ");
			switch (m_vertex_connection_technique) {
				case VertexConnection::TimeSampling:
					fprintf(ITR::f_log, "Point-to-Line\n");
					break;
				case VertexConnection::ShadowConnection:
					fprintf(ITR::f_log, "Shadow Connection\n");
					break;
			}
			fprintf(ITR::f_log, "\n");
		}
	}
	
	void set_mode(SamplingMode mode)
	{
		switch (mode) {
			case Transient:
				m_next_event_sampling_technique = NextEventSampling::Time;
				m_angular_sampling_technique = AngularSampling::ShortestTime;
				m_termination_sampling_technique = TerminationSampling::TotalTime;
				m_vertex_connection_technique = VertexConnection::TimeSampling;
				m_weighting_technique = WeightingTechnique::MIS_Balance;
				break;
			case SteadyState:
				m_next_event_sampling_technique = NextEventSampling::MeanFreePath;
				m_angular_sampling_technique = AngularSampling::PhaseFunction;
				m_termination_sampling_technique = TerminationSampling::RussianRoulette;
				m_vertex_connection_technique = VertexConnection::ShadowConnection;
				m_weighting_technique = WeightingTechnique::MIS_Balance;
		};
	}

	void set_max_path_size(int max_path_size)
	{
		m_max_path_size = max_path_size;
	}

	void set_path_tracing_only(bool path_trace_only)
	{
		m_path_trace_only = path_trace_only;
	}

	void keep_track_on_nb_samples(bool on = true)
	{
		m_keep_track_sample_count = on;
	}

	/* Transient radiance functions */
	void operator()(const VectorN<D> &p, RadSampleList &samples, const Real delta_time,
			const unsigned int nb_time_samples) const;

	virtual void operator()(const Ray<D> &r, RadSampleList &samples, const Real delta_time = 0.,
			const unsigned int nb_time_samples = 1) const override;

	/* Steady-state radiance functions */
	virtual Radiance operator()(const Ray<D> &r) const override;
}; // BidirectionalPathTracing

/*
 * Computes sampling distance. Note that here we assume that the ray has not been
 * traced, and therefore the distance and pdf needs to be corrected afterwards
 * (see function 'correct_next_event_distance').
 */
template<unsigned D, class Radiance, class RadianceAttenuation>
Real BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::sample_next_event_distance(
		const Ray<D> &r, Real &p) const
{
	Real lambda = 0.;
	switch (m_next_event_sampling_technique) {
		case NextEventSampling::Time: {
			Real total_time = ITR::m_film->get_time_resolution() * ITR::m_film->get_exposure_time();
			Real total_distance = total_time / ITR::m_world.time_of_flight(1)
					/ ITR::m_world.get_ior();

			lambda = Real(m_avg_nb_vertices_path) / total_distance;
			}
			break;
		case NextEventSampling::MeanFreePath:
			lambda = r.get_medium()->get_extinction(r.get_origin()).avg();
	}

	// Sample distance
	Real epsilon = Random::StdRNG.next_real();
	Real distance = -log(epsilon) / lambda;

	p = exp(-distance * lambda) * lambda;

	return distance;
}

/*
 * Corrects sampling probabilities when the sampled distance in the medium is
 * beyond the first surface hit.
 *
 * This function assumes that the ray 'r' has already been traced, and that it
 * is within participating media.
 */
template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::correct_next_event_distance(
		Ray<D> &r, Real &distance, bool &media_interaction, Real &p) const
{
	if (!r.get_medium()) {
		distance = r.get_parameter();
		media_interaction = false;
		return;
	}

	if (r.did_hit() && (r.get_parameter() < distance)) {
		distance = r.get_parameter();
		media_interaction = false;

		Real lambda = 0.;
		switch (m_next_event_sampling_technique) {
			case NextEventSampling::Time: {
				Real total_time = ITR::m_film->get_time_resolution()
						* ITR::m_film->get_exposure_time();
				Real total_distance = total_time / ITR::m_world.get_ior();

				lambda = Real(m_avg_nb_vertices_path) / total_distance;
			}
				break;
			case NextEventSampling::MeanFreePath:
				lambda = r.get_medium()->get_extinction(r.get_origin()).avg();
		}

		p = exp(-distance * lambda);
	} else {
		r.set_parameter(distance);
		media_interaction = true;
	}
}

/*
 * Sample next event's angular domain
 */
template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::sample_outgoing_ray(
		Ray<D> &ray, Intersection<D> &it, bool media_interaction,
		RadianceAttenuation &f, Real &time, Real &p) const
{
	if (media_interaction) {
		// Random walk's next step
		VectorN<D> position = ray.get_position();
		VectorN<D> new_omega_i;
		Medium<D>* m = ray.get_medium();
		Spectrum u_s = m->get_scattering(position);

		m->sample_direction(position, ray.get_direction(), new_omega_i, f, time, p);
		f *= u_s;

		ray = Ray<D>(position, new_omega_i, ray.get_level() + 1, ray.get_refraction_index(),
				ray.get_medium());
	} else {
		// SURFACE HIT

		// Random walk's next step
		it.material()->sample_outgoing_ray(it, ray, f, time, p);
		if (!it.material()->is_type(Reflectance::DELTA))
			f *= dot_abs(it.get_normal(), ray.get_direction());
	}
}

/*
 * Sample next event's angular domain
 */
template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::sample_outgoing_ray(
		Ray<D> &ray, Intersection<D> &it, bool media_interaction, Real sampled_distance,
		Real current_delay, const VectorN<D> &target, RadianceAttenuation &f,
		Real &time, Real &p) const
{
	if (media_interaction) {
		// Random walk's next step
		VectorN<D> position = ray.get_position();
		VectorN<D> new_omega_i;
		Medium<D>* m = ray.get_medium();
//		Spectrum u_t = m->get_extinction(position);
		Spectrum u_s = m->get_scattering(position);

		// Select sampling technique to use
		Real p_ts = 0.;
		switch (m_angular_sampling_technique) {
			case AngularSampling::PhaseFunction:
				p_ts = 0.;
				break;
			case AngularSampling::ShortestTime:
				p_ts = 1.;
		}
		
		// Create the angular time sampler
		// - Compute direction and distance to target
		VectorN<D> l = target - position;
		Real mod = l.length();

		// - Compute sampling boundaries
		Real ta = ITR::m_world.time_of_flight(mod) * ITR::m_world.get_ior() + current_delay;
		Real tb = ITR::m_film->get_time_resolution() * ITR::m_film->get_exposure_time()
				+ ITR::m_film->get_time_offset();

		// - And create the sampler
		Probability::TimeAngular<D> sa(ta, tb, sampled_distance, mod, ITR::m_world.get_ior(),
				1. / ITR::m_world.time_of_flight(1.));

		// Check if between boundaries
		if ((tb - (ta + ITR::m_world.time_of_flight(sampled_distance) * ITR::m_world.get_ior()))
				< ITR::m_film->get_exposure_time()) {
			p_ts = 0.;
		}

		Real xi = Random::StdRNG.next_real();

		// MIS between phase function sampling and angular time sampling
		if (xi < p_ts) {
			// Time sampling
			Real theta = sa.cdf_inv(Random::StdRNG.next_real());
			Real phi = 2. * M_PI * Random::StdRNG.next_real(); // Not sure about this

			VectorN<D> dir(cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta));
			new_omega_i = dir.transform_matrix_to(Vector3(0., 1., 0.), l);

			// Compute f and p, with p = p/w using MIS
			m->f(position, ray.get_direction(), new_omega_i, f, time, p);
			f *= u_s;

			p = m->p(position, ray.get_direction(), new_omega_i, time) * (1. - p_ts)
					+ sa.pdf(theta) / (2. * M_PI) * p_ts;
			p /= p_ts;

			if (p < 1.e-5)
				p_ts = 0.;
		}

		// Phase Function Sampling
		if (xi >= p_ts) {
			m->sample_direction(position, ray.get_direction(), new_omega_i, f, time, p);
			f *= u_s;
			if (p_ts > 0.) {
				Real cos_omega_i = Utils::clamp(dot(l, new_omega_i), -1., 1.);

				p = p * (1. - p_ts) + sa.pdf(acos(cos_omega_i)) / (2. * M_PI) * p_ts;
				p /= (1. - p_ts);
			}
		}

		ray = Ray<D>(position, new_omega_i, ray.get_level() + 1, ray.get_ior(),
				ray.get_medium());
	} else {
		// SURFACE HIT

		// Random walk's next step
		it.material()->sample_outgoing_ray(it, ray, f, time, p);
		if (!it.material()->is_type(Reflectance::DELTA))
			f *= dot_abs(it.get_normal(), ray.get_direction());
	}
}

/*
 * Decide the termination of the random walk.
 *	The function decides if the path finishes or not. The method to select the
 *	termination is specified by 'm_termination_sampling_technique' from the
 *	following:
 *		- TerminationSampling::TotalTime - Evaluates if the path total time is
 *		  longer than the total time of the m_film, and if it is, terminates the
 *		  random walk.
 *		- WeightingTechnique::NbBounces - If the number of bounces of the path is
 *		  higher than the max the path is terminated.
 *		- WeightingTechnique::RussianRoulette - Uses stochastic Russian roulette
 *		  to determine whether the path is terminated (absorbed) or not. This is
 *		  the only one unbiased.
 */
template<unsigned D, class Radiance, class RadianceAttenuation>
bool BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::sample_termination(
		Real albedo, Real time, unsigned int nb_bounces, Real &p) const
{
	p = 1.;

	if (nb_bounces > m_max_path_size)
		return true;

	switch (m_termination_sampling_technique) {
		case TerminationSampling::NbBounces:
			if (nb_bounces > m_avg_nb_vertices_path) {
				return true;
			}
			break;
		case TerminationSampling::TotalTime:
			if (ITR::m_film) {
				if (time > ITR::m_film->get_time_length()) {
					return true;
				}
				break;
			}
			// fall through
		case TerminationSampling::RussianRoulette:
			Real epsilon = Random::StdRNG.next_real();

			if (epsilon > albedo) {
				return true;
			}

			p = albedo;
	}

	return false;
}

//=============================================================================
// Compute path
template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::generate_path(const Ray<D> &r,
		const Radiance &f0, const Real p0, const Real t0, const VectorN<D> &target,
		const TraceDirection dir, PathR &path) const
{
	// Prepare variables for the loop
	Ray<D> curr_ray(r.get_origin(), r.get_direction(), false, r.get_level(),
			r.get_ior(), r.get_medium());
	Intersection<D> curr_it;
	VectorN<D> position;
	
	Spectrum u_t, u_s;
	Radiance scat_albedo;
	Medium<D>* previous_medium = nullptr;
	Real avg_albedo;

	// Sampled mean free path in free-space
	Real sampled_distance = std::numeric_limits<Real>::infinity();

	// Throughput and probability of the path
	path.set_value(f0);
	RadianceAttenuation f(1.);
	set_coordinate_system(f0, f);

	Real p = 1., pi = p0;
	Real time = t0;

	// Variables used through the loop
	RadianceAttenuation fa(1.);
	Real pd = 1., pa = 1., ta = 0.;

	bool media_interaction = false;
	bool first_iteration = true;

	// Iterate path
	unsigned iter = 0;
	while (true) {
		iter++;

		u_t = Spectrum(0.);
		u_s = Spectrum(0.);
		
		pd = 1.;
		
		// -----------------------------------------------------------------------------------------
		// First of all test if we are currently in a medium, and sample distance there...	
		if (curr_ray.get_medium()) {
			// If within media...
			u_t = curr_ray.get_medium()->get_extinction(curr_ray.get_origin());
			u_s = curr_ray.get_medium()->get_scattering(curr_ray.get_origin());
			
			// Sample Distance! 
			sampled_distance = sample_next_event_distance(curr_ray, pd);
		}

		// -----------------------------------------------------------------------------------------
		// With the sampled distance, sample the new ray
		fa = RadianceAttenuation(1.);
		pa = 1., ta = 0.;
		previous_medium = curr_ray.get_medium();
		
		if (first_iteration) {
			first_iteration = false;
		} else {
			sample_outgoing_ray(curr_ray, curr_it, media_interaction, sampled_distance, time,
					target, fa, ta, pa);
		}

		// -----------------------------------------------------------------------------------------
		// Trace next ray
		curr_it = Intersection<D>();
		this->m_world.first_intersection(curr_ray, curr_it);

		// Test if there's a scattering event in media or surfaces. If not then the
		// path is finished...
		if (!curr_ray.get_medium() && !curr_it.did_hit())
			return;

		// -----------------------------------------------------------------------------------------
		// In case of a surface interaction, it mights occur that the ray is now interacting with a
		// new media. Thus, we need to take into account the new media...
		if (previous_medium != curr_ray.get_medium()) {
			if (curr_ray.get_medium()) {
				// If within media...
				u_t = curr_ray.get_medium()->get_extinction(curr_ray.get_origin());
				u_s = curr_ray.get_medium()->get_scattering(curr_ray.get_origin());

				// Sample Distance! 
				sampled_distance = sample_next_event_distance(curr_ray, pd);
			} else {
				u_t = Spectrum(0.);
				u_s = Spectrum(0.);
				
				pd = 1.;
				// Sampled mean free path in free-space
				sampled_distance = std::numeric_limits<Real>::infinity();
			}
		}

		// -----------------------------------------------------------------------------------------
		// And correct distances 
		correct_next_event_distance(curr_ray, sampled_distance, media_interaction, pd);

		// -----------------------------------------------------------------------------------------
		// Update contribution, probability, accumulated probability, and time.
		pi *= (pd * pa);
		p *= pi;

		if (dir & PATH_TRACING) {
			set_attenuation_tracing_direction(dir, fa);
			try {
				f = (f * fa * exp(-sampled_distance * u_t)) / pi;
			} catch (std::runtime_error& e) {
				fprintf(stderr, "Failed to connect path!\n");
				f.print(stderr);
				fa.print(stderr);

				fprintf(stderr, "--------------------------------\n");
				fprintf(stderr, "PATH:\n");
				for (size_t i = 0; i < path.size(); ++i) {
					path[i].get_vertex_value().print(stderr);
				}

				throw e;
			}
		} else {
			f = (f * fa * exp(-sampled_distance * u_t)) / pi;
		}

		p = 1.;
		time += ta + ITR::m_world.time_of_flight(sampled_distance)*curr_ray.get_ior();

		// -----------------------------------------------------------------------------------------
		// Store vertex
		if (media_interaction) {
			position = curr_ray.get_origin() + curr_ray.get_direction() * sampled_distance;

			path.add_media_vertex(position, curr_ray.get_direction(), curr_ray.get_medium(), f, p, pi,
					time);
			
			scat_albedo = (u_s / u_t);
		} else {
			position = curr_it.get_position();

			if (!curr_it.material()->is_type(Reflectance::DELTA)) {
				path.add_surface_vertex(position, curr_ray.get_direction(), curr_it.get_normal(),
						curr_it.get_uv(), curr_it.material(), curr_ray.get_medium(), f, p, pi, time);
			}

			scat_albedo = (1. - curr_it.material()->get_absorption(curr_it));
		}

		// -----------------------------------------------------------------------------------------
		// Finally, evaluate Sample Termination!
		avg_albedo = scat_albedo.avg();

		if (sample_termination(avg_albedo, time, curr_ray.get_level(), pi))
			return;
	}
}

//=============================================================================
// Connect an eye vertex with a randomly sampled light.
// TODO: Make use of "connect_vertices" to use time sampling methods...
template<unsigned D, class Radiance, class RadianceAttenuation>
bool BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::connect_vertex_with_light(
		const VertexR& eye_vertex, Radiance &f, Real &p, Real *time) const
{
	LightSample<D, Radiance> light_sample;
	Real pl1 = 1., pl2 = 1.;

	const VectorN<D>& vertex_position = eye_vertex.get_vertex_position();

	if (ITR::m_world.sample_light(pl1)->sample(vertex_position, light_sample, pl2)) {
		RadianceAttenuation fe;

		Real pe = 1., te = 0.;
		if (time != nullptr) {
			eye_vertex.compute_scattering(-light_sample.dir, fe, te, pe);
			// Careful, need to change the "get_ior" shit...
			(*time) = ITR::m_world.time_of_flight(light_sample.dist) * ITR::m_world.get_ior() + te
					+ light_sample.instant;
		} else {
			eye_vertex.compute_scattering(-light_sample.dir, fe);
		}

		// Compute attenuation if needed
		Spectrum att = eye_vertex.compute_attenuation(
				Ray<D>(light_sample.dist, light_sample.pos, light_sample.dir));

		set_attenuation_tracing_direction(TraceDirection::FROM_EYE, fe);

		f = fe * (light_sample.irradiance * att);
		p = pl1 * pl2 * pe;

		return true;
	} else {
		f = Radiance(0.);
		p = 1.;

		if (time != nullptr) {
			(*time) = ITR::m_world.time_of_flight(light_sample.dist) * ITR::m_world.get_ior()
					+ light_sample.instant;
		}

		return false;
	}
}

//=============================================================================
// Connect two vertices
template<unsigned D, class Radiance, class RadianceAttenuation>
bool BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::connect_vertices(
		const VertexR& eye_vertex, const VertexR& light_vertex, RadianceAttenuation &f,
		Real *time, Real *pdf_time) const
{
	VectorN<D> l2e = eye_vertex.get_vertex_position() - light_vertex.get_vertex_position();
	Real d_l2e = l2e.length();
	l2e /= d_l2e;
	
	Ray<D> ray;

	if (time
			&& (d_l2e + eye_vertex.get_subpath_delay() + light_vertex.get_subpath_delay())
					> ITR::m_film->get_time_resolution() * ITR::m_film->get_exposure_time()
							+ ITR::m_film->get_time_offset()) {
		return false;
	}

	switch (m_vertex_connection_technique) {
		//-----------------------------------------------------------------------------
		case VertexConnection::TimeSampling:
			if ((eye_vertex.get_type() & PathR::Vertex::Type::MEDIUM) && time && ITR::m_film) {
				// Only makes sense if we are in a medium, and interested in time.
				// Otherwise, we use regular MFP sampling.

				TraceDirection connection_direction = TraceDirection::UNKNOWN;

				RadianceAttenuation fo, ft, fv;
				Real po, pt, pv, psp, pd;
				Real to, tt, tv;

				const typename PathR::Vertex* o;
				const typename PathR::Vertex* t;

				// We first obtain the intermediate vertex
				// - Select if the origin is the light or the eye vertex.
				if (Random::StdRNG.next_real() < .5) {
					o = &eye_vertex;
					t = &light_vertex;
					l2e = -l2e;
					connection_direction = TraceDirection::FROM_EYE;
				} else {
					o = &light_vertex;
					t = &eye_vertex;
					connection_direction = TraceDirection::FROM_LIGHT;
				}
				psp = .5;

				// - Sample a new direction
				VectorN<D> dir;
				o->sample_direction(dir, fo, to, po);

				// - Sample distance using Time Sampling!
				Real dist = -1.;

				Real ta = ITR::m_world.time_of_flight(d_l2e) * ITR::m_world.get_ior()
						+ o->get_subpath_delay() + t->get_subpath_delay();

				Real tb = ITR::m_film->get_time_resolution() * ITR::m_film->get_exposure_time()
						+ ITR::m_film->get_time_offset();
				if (ta >= tb) {
					return false;
				}

				Probability::TimeLineToPoint<D> st(ta, tb, dir, l2e, ITR::m_world.get_ior(),
						1. / ITR::m_world.time_of_flight(1));
				while (dist < 0.) {
					dist = st.cdf_inv(Random::StdRNG.next_real());
				}
				pd = st.pdf(dist);

				// - Check if hitting a surface
				Intersection<D> it;
				ray = Ray<D>(o->get_vertex_position(), dir, true);

				bool media_interaction = true;
				ITR::m_world.first_intersection(ray, it);

				if (ray.did_hit() && (ray.get_parameter() < dist)) {
					dist = ray.get_parameter();
					pd = st.cdf(dist);
					media_interaction = false;
				}

				typename PathR::Vertex v;

				// - Create the new vertex
				if (media_interaction) {
					v = VertexR(o->get_vertex_position() + dir * dist, dir, ray.get_medium(),
							RadianceAttenuation(1.), 1., 1., 0.);
				} else {
					if (!it.material()->is_type(Reflectance::DELTA)) {
						v = VertexR(it.get_position(), dir, it.get_normal(), it.get_uv(),
								it.material(), ray.get_medium(), RadianceAttenuation(1.),
								1., 1., 0.);
					} else {
						return false;
					}
				}

				// Perform shadow connection
				VectorN<D> v2r = t->get_vertex_position() - v.get_vertex_position();
				Real d_v2r = v2r.length();
				v2r /= d_v2r;

				ray = Ray<D>(t->get_vertex_position(), v2r, false);

				if (ITR::m_world.intersects(ray, d_v2r)) {
					v.compute_scattering(v2r, fv, tv, pv);
					t->compute_scattering(-v2r, ft, tt, pt);
					
					(*time) = ITR::m_world.time_of_flight(dist) * ITR::m_world.get_ior()
							+ ITR::m_world.time_of_flight(d_v2r) * ITR::m_world.get_ior() + to + tv
							+ tt;

					(*pdf_time) = psp * po * pd * pv * pt;

					Spectrum att = eye_vertex.compute_attenuation(eye_vertex.get_vertex_position(),
							light_vertex.get_vertex_position());
					att /= ((D == 2) ? d_v2r : (d_v2r * d_v2r));

					if (connection_direction && TraceDirection::FROM_EYE) {
						set_attenuation_tracing_direction(TraceDirection::FROM_EYE, fo);
						set_attenuation_tracing_direction(TraceDirection::FROM_EYE, fv);
						f = (fo * fv * ft) * att;
					} else {
						set_attenuation_tracing_direction(TraceDirection::FROM_EYE, ft);
						f = (ft * fv * fo) * att;
					}

					return true;
				}

				// else...
				f = RadianceAttenuation();
				(*time) = ITR::m_world.time_of_flight(d_l2e) * ITR::m_world.get_ior();
				(*pdf_time) = 1.;

				return false;
			}
			// fall through
		//-----------------------------------------------------------------------------
		// Regular shadow connection
		case VertexConnection::ShadowConnection:
		default:
			RadianceAttenuation fe, fl;
			Real pe, pl;
			Real te, tl;

			ray = Ray<D>(light_vertex.get_vertex_position(), l2e, false);

			if (ITR::m_world.intersects(ray, d_l2e)) {
				if (time) {
					eye_vertex.compute_scattering(-l2e, fe, te, pe);
					light_vertex.compute_scattering(l2e, fl, tl, pl);

					(*time) = ITR::m_world.time_of_flight(d_l2e) * ITR::m_world.get_ior() + te + tl;

					(*pdf_time) = pe * pl;
				} else {
					eye_vertex.compute_scattering(-l2e, fe);
					light_vertex.compute_scattering(l2e, fl);
				}

				// Compute attenuation if needed
				Spectrum att = eye_vertex.compute_attenuation(eye_vertex.get_vertex_position(),
						light_vertex.get_vertex_position());

				set_attenuation_tracing_direction(TraceDirection::FROM_EYE, fe);
				set_attenuation_tracing_direction(TraceDirection::FROM_LIGHT, fl);

				f = fe * fl * (att / (D == 2 ? d_l2e : (d_l2e * d_l2e)));

				return true;
			} else {
				f = RadianceAttenuation();
				if (time) {
					(*time) = ITR::m_world.time_of_flight(d_l2e) * ITR::m_world.get_ior();
					(*pdf_time) = 1.;
				}

				return false;
			}
	}
}

//=============================================================================
// Compute the weight of the path.
//	The function weights the contribution of the path resulting on connecting
//	the subpaths eye_path[ie] and light_path[il], using the technique defined 
//	in 'm_weighting_technique', which can be:
//		- WeightingTechnique::Uniform - It applies a uniform weight to all
//			paths (w=1).
//		- WeightingTechnique::Balance - Computes the weight using MIS, using
//			the balance heuristic [1], as described in Veach's thesis [2].
//		- WeightingTechnique::Power - Computes the weight using MIS with the
//			the power heuristic [1,2], using the parameter 'm_beta_MIS'.
//
//   [1] E. Veach and L. Guibas. 1995. Optimally Combining Sampling Techniques 
//		for Monte Carlo Rendering. In SIGGRAPH '95.
//      http://www-graphics.stanford.edu/papers/combine/
//   [2] E. Veach. 1997. Robust Monte Carlo Methods for Light Transport
//      Simulation. PhD Thesis (Stanford University).
//      http://graphics.stanford.edu/papers/veach_thesis/
//
template<unsigned D, class Radiance, class RadianceAttenuation>
Real BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::compute_weights(
		const PathR &eye_path, const PathR &light_path, const unsigned int ie,
		const unsigned int il) const
{
	Real beta = m_beta_MIS;

	switch (m_weighting_technique) {
		//------------------------------------------------------------------------
		// MIS
		case WeightingTechnique::MIS_Balance:
			beta = 1.;
			// fall through
		case WeightingTechnique::MIS_Power:
			return PathMIS::compute_weight(eye_path, light_path, ie, il, beta);
			break;
		//------------------------------------------------------------------------
		// Uniform
		case WeightingTechnique::Uniform:
		default:
			return 1.;
			break;
	}
}

// Connect paths in steady state
template<unsigned D, class Radiance, class RadianceAttenuation>
Radiance BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::connect_paths(
		const PathR &eye_path, const PathR &light_path) const
{
	Radiance radiance(0.);
	Radiance f(0.);

	Real p = 1., w = 1.;

	// Iterate over all eye vertices
	for (size_t ie = 0; ie < eye_path.size(); ie++) {
		if ((ie + 1) > m_max_path_size)
			break;

		RadianceAttenuation fe = eye_path[ie].get_vertex_value();

		// First vertex of light path is the light source.
		// Following [Veach 95] we avoid correlation by randomly sampling the direct light
		// for each vertex of the eye path.
		if (((ITR::m_components & Scattering::Level::SINGLE) || ie > 0)
				&& connect_vertex_with_light(eye_path[ie], f, p)) {
			f = fe * f;
			p *= eye_path[ie].get_subpath_pdf();
			
			radiance += f / p;
		}

		if (!(ITR::m_components & Scattering::Level::MULTIPLE))
			break;

		// Iterate over all light vertices
		for (size_t il = 0; il < light_path.size(); il++) {
			RadianceAttenuation fv(0.);
			if (connect_vertices(eye_path[ie], light_path[il], fv)) {
				// Note that for polarization or fluorescense the order of the operations matter.
				f = fe * (fv * (light_path[il].get_vertex_value() * light_path.get_value()));
				p = (eye_path[ie].get_subpath_pdf()) * (light_path[il].get_subpath_pdf());
				w = compute_weights(eye_path, light_path, ie, il);

				radiance += f * (w / p);
			}
		}
	}

	return radiance;
}

// Connect paths in transient state
template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::connect_paths(
		const PathR &eye_path, const PathR &light_path, RadSampleList &samples, const Real,
		const unsigned int) const
{
	Radiance f;
	RadianceAttenuation fv(0.);
	Real p = 1.0;
	Real w = 1.0;
	Real t = 0.0;

	Real total_time = ITR::m_film->get_time_length();
	
	// Iterate over all eye vertices
	for (size_t ie = 0; ie < eye_path.size(); ie++) {
		// If the current vertex delay is longer than the m_film, all its events
		// and these in the subsequent vertices of the path will happen outside
		// the sensor window.
		if (eye_path[ie].get_subpath_delay() > total_time)
			return;

		if ((ie + 1) > m_max_path_size)
			return;

		RadianceAttenuation fe = eye_path[ie].get_vertex_value();

		// First vertex of light path is the light source.
		// Following [Veach 95] we avoid correlation by randomly sampling the direct light
		// for each vertex of the eye path.
		bool connected = connect_vertex_with_light(eye_path[ie], f, p, &t);

		if (connected && ((ITR::m_components & Scattering::Level::SINGLE) || ie > 0)) {
			f = fe * f;
			p *= eye_path[ie].get_subpath_pdf();
			t += eye_path[ie].get_subpath_delay();

			samples.push_back(RadianceSampleR(f / p * ITR::m_inv_incoming_samples, t, ie));
		} else if (m_keep_track_sample_count) {
			t += eye_path[ie].get_subpath_delay();

			samples.push_back(RadianceSampleR(Radiance(0.), t, ie));
		}
		
		if (!(ITR::m_components & Scattering::Level::MULTIPLE))
			break;

		// Iterate over all light vertices
		for (size_t il = 0; il < light_path.size(); il++) {
			if ((ie + il + 2) > m_max_path_size)
				break;
			
			// Check if within the animation...	
			if ((eye_path[ie].get_subpath_delay() + light_path[il].get_subpath_delay())
					> total_time)
				break;

			if (connect_vertices(eye_path[ie], light_path[il], fv, &t, &p)) {
				// Note that for polarization or fluorescence the order of the operations matter.
				f = fe * (fv * (light_path[il].get_vertex_value() * light_path.get_value()));
				p = eye_path[ie].get_subpath_pdf() * light_path[il].get_subpath_pdf();
				t += eye_path[ie].get_subpath_delay() + light_path[il].get_subpath_delay();
				w = compute_weights(eye_path, light_path, ie, il);
				samples.push_back(
						RadianceSampleR(f / (p * w) * ITR::m_inv_incoming_samples, t, ie + il + 1));
			} else if (m_keep_track_sample_count) {
				t += eye_path[ie].get_subpath_delay() + light_path[il].get_subpath_delay();
				samples.push_back(RadianceSampleR(Radiance(0.), t, ie + il + 1));
			}
		}
	}
}

//=============================================================================
// Computes the light transport from an initial camera-sampled ray r.
template<unsigned D, class Radiance, class RadianceAttenuation>
Radiance BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::operator()(
		const Ray<D> &r) const
{
	Radiance reflected_radiance(0.);

	// Preallocate a conservative amount of vertices
    unsigned n_eye = std::min<unsigned>(m_max_path_size, (unsigned) (MAX_BPT_VERTICES));
	unsigned n_light = (!m_path_trace_only) ? n_eye : 0;

	PathR eye_path(n_eye), light_path(n_light);

	for (size_t sample = 0; sample < this->m_incoming_samples; sample++) {
		// Sample light source. We do this first since it can be used
		// to guide the eye path when:
		// m_angular_sampling_technique = AngularSampling::ShortestTime.
		LightSample<D, Radiance> light_sample;
		Real pl1 = 0., pl2 = 0.;
		if (m_angular_sampling_technique == AngularSampling::ShortestTime || !m_path_trace_only) {
			ITR::m_world.sample_light(pl1)->sample(light_sample, pl2);
		}
		
		// Generate eye path
		generate_path(r, Radiance(1.), 1., 0., light_sample.pos, TraceDirection::FROM_EYE,
				eye_path);

		if (eye_path.size() == 0)
			continue;

		// Generate light path
		if (!m_path_trace_only) {
			generate_path(Ray<D>(light_sample.pos, light_sample.dir, false),
					light_sample.irradiance, pl1 * pl2, light_sample.instant, r.get_origin(),
					TraceDirection::FROM_LIGHT, light_path);
		}

		// Connect vertices
		Radiance r1 = connect_paths(eye_path, light_path);
		
		//Average...
		reflected_radiance += r1;

		// Clear paths
		eye_path.clear();
		light_path.clear();
	}

	return reflected_radiance * ITR::m_inv_incoming_samples;
}

template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::operator()(const Ray<D> &r,
		RadSampleList &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	//Spectrum reflected_radiance(0.); // FIX
	
	// Preallocate a conservative amount of vertices
    unsigned n_eye = std::min<unsigned>(m_max_path_size, (unsigned) (MAX_BPT_VERTICES));
	unsigned n_light = (!m_path_trace_only) ? n_eye : 0;

	PathR eye_path(n_eye), light_path(n_light);

	for (size_t sample = 0; sample < ITR::m_incoming_samples; sample++) {
		// Sample light source. We do this first since it can be used
		// to guide the eye path when:
		// m_angular_sampling_technique = AngularSampling::ShortestTime.
		LightSample<D, Radiance> light_sample;
		Real pl1 = 0.0, pl2 = 0.0;
		if (m_path_trace_only || (m_angular_sampling_technique == AngularSampling::ShortestTime))
			ITR::m_world.sample_light(pl1)->sample(light_sample, pl2);
		
		// Generate eye path
		generate_path(r, Radiance(1), 1, 0, light_sample.pos, TraceDirection::FROM_EYE, eye_path);
		if (eye_path.size() == 0)
			continue;

		// Generate light path
		if (!m_path_trace_only) {
			// Need to fix the IOR and Medium of this path!
			generate_path(Ray<D>(light_sample.pos, light_sample.dir, false),
					light_sample.irradiance, pl1 * pl2, light_sample.instant, r.get_origin(),
					TraceDirection::FROM_LIGHT, light_path);
		}

		// Connect vertices
		connect_paths(eye_path, light_path, samples, delta_time, nb_time_samples);

		// Clear paths
		eye_path.clear();
		light_path.clear();
	}
}

template<unsigned D, class Radiance, class RadianceAttenuation>
void BidirectionalPathTracing<D, Radiance, RadianceAttenuation>::operator()(const VectorN<D> &p,
		RadSampleList &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	// Preallocate a conservative amount of vertices
    unsigned n_path = std::min(m_max_path_size, (unsigned) (MAX_BPT_VERTICES));
	PathR eye_path(n_path), light_path(n_path);

	for (size_t sample = 0; sample < ITR::m_incoming_samples; sample++) {
		// Sample light source. We do this first since it can be used
		// to guide the eye path when:
		// m_angular_sampling_technique = AngularSampling::ShortestTime.
		LightSample<D, Radiance> light_sample;
		Real pl1, pl2;
		ITR::m_world.sample_light(pl1)->sample(light_sample, pl2);

		// Generate eye path
		// - sample eye ray
		VectorN<D> dir_ray;
		Real pdf_ray;
		Sampling.direction_uniformly(dir_ray, pdf_ray);
		// - ...and generate the path
		generate_path(Ray<D>(p, dir_ray, false), Spectrum(ITR::m_inv_incoming_samples), pdf_ray,
				light_sample.pos, eye_path);

		if (eye_path.size() == 0)
			continue;

		// Generate light path
		generate_path(Ray<D>(light_sample.pos, light_sample.dir, false),
				light_sample.irradiance * ITR::m_inv_incoming_samples, pl1 * pl2, p, light_path);
		
		// Connect vertices
		connect_paths(eye_path, light_path, samples, delta_time, nb_time_samples);

		// Clear paths
		eye_path.clear();
		light_path.clear();
	}
}

#endif // _BIDIRECTIONAL_PATH_TRACING_INTEGRATOR_H_

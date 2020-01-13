/*
 * Copyright (C) 2017, Julio Marco (http://webdiis.unizar.es/~juliom/)
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

#ifndef _PARTICLE_TRACING_H_
#define _PARTICLE_TRACING_H_

#include "bunnykiller.h"

#include <list>
#include <cstdlib>

#include "Integrator/Integrator.h"
#include "RayTracing/World.h"
#include "RayTracing/Intersection.h"
#include "Integrator/BidirectionalPathTracing.h"
#include "Integrator/Particle.h"

/**	Virtual class to implement all techniques based on particle tracing,
	including Photon Mapping [Jensen 2001] and Instant Radiosity 
	[Keller 1997], and the algorithms that inherit from them.	*/
template<unsigned D, class Radiance, class RadianceAttenuation>
class ParticleTracing: public BidirectionalPathTracing<D, Radiance, RadianceAttenuation>
{
protected:
	using BDPTR = BidirectionalPathTracing<D, Radiance, RadianceAttenuation>;
	using RadianceSampleR = typename BDPTR::RadianceSampleR;
	using RadSampleList = typename BDPTR::RadSampleList;
	using FilmR = typename BDPTR::FilmR;
	using PathR = typename BDPTR::PathR;
	using VertexR = typename BDPTR::VertexR;
protected:
	using ParticleR = Particle<D, Radiance, RadianceAttenuation>;
protected:
	Real MAX_LONGBEAM_DIST;
	int m_max_nb_shots;
	int m_max_global_spoints, m_max_caustic_spoints;
	int m_max_global_vpoints, m_max_caustic_vpoints;
	int m_max_global_longbeams, m_max_caustic_longbeams;
	int m_max_global_shortbeams, m_max_caustic_shortbeams;
	Real m_spoint_radius, m_vpoint_radius, m_shortbeam_radius, m_longbeam_radius;
	Real m_max_time;
	bool m_use_russian_roulette;
	typename ParticleR::ParticleType m_particleTypes;
	typename Scattering::Level m_scatt_level;

	/**	Trace particles and store their intersections	*/
	void trace_photons(std::vector<ParticleR> &global_particles, std::vector<ParticleR> &caustic_particles) const;
public:
	struct TracingSetup{
		int max_shots;
		int max_spoints_global;
		int max_vpoints_global;
		int max_shortbeams_global;
		int max_longbeams_global;
		float spoint_radius;
		float vpoint_radius;
		float shortbeam_radius;
		float longbeam_radius;
		typename Scattering::Level scatt_level;
		int max_nbounces;
		typename ParticleR::ParticleType p_types;
	};

	ParticleTracing(World<D, Radiance> *w, TracingSetup ts, FilmR *film) :
		BDPTR(w, 1, NULL, film),
		m_max_nb_shots(ts.max_shots),
		m_max_global_spoints(ts.max_spoints_global),
		m_max_global_vpoints(ts.max_vpoints_global),
		m_max_global_shortbeams(ts.max_shortbeams_global),
		m_max_global_longbeams(ts.max_longbeams_global),
		m_spoint_radius(ts.spoint_radius),
		m_vpoint_radius(ts.vpoint_radius),
		m_shortbeam_radius(ts.shortbeam_radius),
		m_longbeam_radius(ts.longbeam_radius),
		m_particleTypes(ts.p_types),
		m_scatt_level(ts.scatt_level),
		MAX_LONGBEAM_DIST(1000),
		m_max_time(1e20),
		m_use_russian_roulette(true)
	{
		switch (m_scatt_level) {
			case Scattering::SINGLE:
				this->set_mode(2);
				this->set_max_path_size(ts.max_nbounces);
				break;
			case Scattering::ALL:
				this->set_max_path_size(30);
				break;
		}
	}

	virtual ~ParticleTracing() {};

	/** Virtual functions to store or load the particles. Each method
		store the particles in a different way, depending if they
		are Photons, VPLs, Rays or even full hierarchies.	*/
	virtual void store(const std::string &namefile)const = 0;
	virtual void load(const std::string &namefile)const = 0;

	// Redefinition of generate path to store vertices at delta surfaces
	void generate_path(const Ray<D> &r, const Radiance &f0, const Real p0, const VectorN<D> &target, PathR &path ) const;

	bool do_trace (typename ParticleR::ParticleType _t) const
	{
		return (m_particleTypes & _t) == _t;
	}

	void set_max_time(Real max_time)
	{
		m_max_time = max_time;
		m_use_russian_roulette = false;
	}
}; //ParticleTracing

//template<unsigned D, class Radiance, class RadianceAttenuation>
//void ParticleTracing<D, Radiance, RadianceAttenuation>::generate_path(const Ray<D> &r, const Radiance &f0, const Real p0, const VectorN<D> &target, PathR &path )const
//{
//	//Prepare variables for the loop
//	Ray<D> curr_ray(r.get_origin(),r.get_direction(), false, r.get_level(),
//			        r.get_refraction_index(), r.get_medium());
//	Intersection<D> curr_it;
//	VectorN<D> position;
//
//	Spectrum u_t, u_s, scat_albedo;
//	Real avg_albedo;
//	Real sampled_distance = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
//
//	// Throughput and probability of the path
//	RadianceAttenuation f = f0;
//	Real p = 1., pi = p0;
//
//	Real time = 0;
//
//	// Variables used through the loop
//	RadianceAttenuation fa(1.);
//	Real pd(1.), pa(1.);
//	Real ta = 0.;
//
//
//	bool media_interaction = false;
//	bool first_iteration = true;
//
//	// Iterate path
//	while(1)
//	{
//		u_t = Spectrum(0.);
//		u_s = Spectrum(0.);
//
//		pd = 1.;
//		// ------------------------------------------------------------------------------------------
//		// First of all, evaluate surface reflection, in case we need to switch to a different media.
//		// Note that this is suboptimal for rough interfaces, since we don't apply time-sampling...
//		fa = Spectrum(1.); pa=1.; ta=0.;
//		if( !( first_iteration || media_interaction) )
//		{
//			fa = curr_it.material()->sample_outgoing_ray(curr_it, curr_ray, ta, pa );
//			fa *= dot_abs(curr_it.get_normal(), curr_ray.get_direction());
//		}
//
//
//		// ------------------------------------------------------------------------------------------
//		// First of all test if we are currently in a medium, and sample distance there...
//		if(curr_ray.get_medium())
//		{
//			// If within media...
//			u_t = curr_ray.get_medium()->get_extinction(curr_ray.get_origin());
//			u_s = curr_ray.get_medium()->get_scattering(curr_ray.get_origin());
//
//			// Sample Distance!
//
//			sampled_distance = sample_next_event_distance(curr_ray, pd);
//			//media_interaction = true;
//		}
//		// ------------------------------------------------------------------------------------------
//		// With the sampled distance, sample the new ray
//		fa = Spectrum(1.); pa=1.; ta=0.;
//		if( first_iteration )
//			first_iteration = false;
//		else
//			sample_outgoing_ray( curr_ray, curr_it, media_interaction, sampled_distance, time, target, fa, ta, pa );
//
//
//		// ------------------------------------------------------------------------------------------
//		// Trace next ray
//		curr_it = Intersection<D>();
//		this->m_world->first_intersection(curr_ray, curr_it);
//
//
//		// Test if there's no scattering event in media or surfaces. If don't then the
//		// path is finished...
//		if( !curr_ray.get_medium() && !curr_it.did_hit())
//			return;
//
//		// ------------------------------------------------------------------------------------------
//		// And correct distances
//		correct_next_event_distance( curr_ray, sampled_distance, media_interaction, pd );
//
//		// ------------------------------------------------------------------------------------------
//		// Evaluate if the sampled distance is in the infinity
//		if( sampled_distance == std::numeric_limits<Real>::infinity() )
//		return;
//
//		// ------------------------------------------------------------------------------------------
//		// Update contribution, probability, accomulated probability, and time.
//		f *= (fa*exp(-sampled_distance*u_t));
//		//f *= (media_interaction) ? 1 : dot_clamped(curr_it.get_normal(), -curr_ray.get_direction());
//		pi *= (pd*pa);
//		p *= pi;
//		time += ta + BDPTR::m_world->time_of_flight(sampled_distance)*curr_ray.get_refraction_index();;
//
//		// ------------------------------------------------------------------------------------------
//		// Store vertex
//		if (media_interaction) {
//			position = curr_ray.get_origin() + curr_ray.get_direction()*sampled_distance;
//			path.add_vertex(position, curr_ray.get_direction(), curr_ray.get_medium(),
//							f, p, pi, time);
//
//			scat_albedo = (u_s/u_t);
//			if (m_scatt_level == Scattering::SINGLE)
//				return;
//		} else {
//			position = curr_it.get_position();
//			//if( !curr_it.material()->is_type(Reflectance::DELTA) )
//				path.add_vertex(position, curr_ray.get_direction(), curr_it.get_normal(), curr_it.get_uv(), curr_it.material(),
//							f, p, pi, time);
//
//			scat_albedo = (1.-curr_it.material()->get_absorption(curr_it.get_uv()));
//		}
//
//		// ------------------------------------------------------------------------------------------
//		// Finally, evaluate Sample Termination!
//		pi = 1.f;
//		avg_albedo = scat_albedo.avg();
//
//		if (sample_termination(avg_albedo, time, curr_ray.get_level(), pi))
//			return;
//	}
//}

#if 1
template<unsigned D, class Radiance, class RadianceAttenuation>
void ParticleTracing<D, Radiance, RadianceAttenuation>::trace_photons(std::vector<ParticleR> &global_particles, std::vector<ParticleR> &caustic_particles) const
{
	unsigned int nb_shots = 0;
	unsigned int nb_shots_caustics = 0;
	unsigned int nb_shots_global = 0;
	bool shooting = true;
	int j = 0;
	unsigned int next_id = 0;
	while (shooting) {
		if (++nb_shots > m_max_nb_shots) {
			nb_shots_caustics = nb_shots_caustics ?
					nb_shots_caustics :
					(nb_shots - 1);
			nb_shots_global = nb_shots_global ?
					nb_shots_global :
					(nb_shots - 1);
			break;
		}

		// Generate light path
		PathR path;
		LightSample<D, Radiance> light_sample;
		Real pl1, pl2;
		BDPTR::m_world->sample_light(pl1)->sample( light_sample, pl2);

		generate_path(Ray<D>(light_sample.pos, light_sample.dir,false,0, BDPTR::m_world->get_ior(), BDPTR::m_world->get_medium()),
							 light_sample.irradiance, pl1*pl2, VectorN<D>(0.), path);
		
		for (typename PathR::Vertex& v : path) {
			if (do_trace(ParticleR::ShortBeam)) { // BUILD SHORT BEAM
				// Scatter through outgoing direction, update irradiance and probability
				RadianceAttenuation f0;
				Real t0, p0;
				v0.compute_scattering(v.m_direction, f0, t0, p0);
				v0.m_f *= f0;
				v0.m_p *= p0;
				v0.m_pi *= p0;

				// Change outgoing direction
				v0.m_direction = path[i]->m_direction; // in first vertex, it equal to light_sample.dir

				//Compute track length for short beam
				VectorN<D> vc = path[i]->get_vertex_position() - v0->get_vertex_position();
				Real tracklength = vc.length();

				global_particles.push_back(ParticleR(v, i, ParticleR::ShortBeam, m_shortbeam_radius, tracklength, 0.0f, next_id++));
			}
		}

		for (int i=0; i < path.size(); i++){
			/*if (get_surface_vertex(path[i])) // Vertex is on surface
			{
				//Surface photon point
				const Path<D>::SurfaceVertex *v = get_surface_vertex(path[i]);
				if (do_trace(Particle<D>::SurfPoint) && i > 0)
				{
					global_particles.push_back(Particle<D>(v,i,Particle<D>::SurfPoint));
				}
			}
			else if (get_volume_vertex(path[i])) // Vertex is on medium
			{
				//Medium photon point
				const Path<D>::VolumeVertex *v = get_volume_vertex(path[i]);
				if ( do_trace(Particle<D>::VolPoint) )
				{
					global_particles.push_back(Particle<D>(v,i,Particle<D>::VolPoint, m_vpoint_radius));
					//if( global_particles.size() == m_max_global_spoints )
					//	{	nb_shots_global = nb_shots; break;	}
				}
			}*/

			//No matter vertex is on medium or surface, build photon beam crossing medium
			if (do_trace(ParticleR::ShortBeam)) { // BUILD SHORT BEAM
				typename PathR::Vertex *v0;
				if (i==0)
					v0 = new typename PathR::Vertex(light_sample.pos, light_sample.dir, m_world->get_medium(), light_sample.irradiance, pl1*pl2, pl1*pl2);
				else if (get_volume_vertex(path[i-1])) 
					v0 = new PathR::VolumeVertex(*get_volume_vertex(path[i-1]));
				else if (get_surface_vertex(path[i-1]))
				{
					const PathR::SurfaceVertex *v = get_surface_vertex(path[i-1]);
					// TO-DO: Fix for transmissive materials, check if subpath is inside object and change medium
					v0 = new PathR::VolumeVertex(v->m_position, v->m_direction, m_world->get_medium(), v->get_vertex_value(), v->get_subpath_pdf(), v->get_vertex_pdf());
				}
				
				// Scatter through outgoing direction, update irradiance and probability
				Radiance f0;
				Real t0, p0; 
				v0->compute_scattering(path[i]->m_direction, f0, t0, p0);
				v0->m_f *= f0;
				v0->m_p *= p0;
				v0->m_pi *= p0;

				// Change outgoing direction
				v0->m_direction = path[i]->m_direction; // in first vertex, it equal to light_sample.dir

				//Compute track length for short beam
				VectorN<D> vc = path[i]->get_vertex_position() - v0->get_vertex_position();
				Real tracklength = vc.length();

				global_particles.push_back(ParticleR(v0, i, ParticleR::ShortBeam, m_shortbeam_radius, tracklength, 0.0f, next_id++));
			}

			if ( do_trace(ParticleR::LongBeam) )
			{	
				PathR::VolumeVertex *v0 = NULL;

				Radiance f0;
				Real t0, p0; 
				if (i==0)
				{	
					v0 = new Path<D>::VolumeVertex(light_sample.pos/*+light_sample.dir*this->m_longbeam_radius*/, light_sample.dir, m_world->get_medium(), light_sample.irradiance, pl1*pl2, pl1*pl2);
				}
				/*else if (get_volume_vertex(path[i-1])) 
				{
					if (fabs(get_volume_vertex(path[i-1])->get_vertex_value().avg()/get_volume_vertex(path[i-1])->get_subpath_pdf()) == std::numeric_limits<float>::infinity()
						||get_volume_vertex(path[i-1])->get_vertex_value().avg()/get_volume_vertex(path[i-1])->get_subpath_pdf() != get_volume_vertex(path[i-1])->get_vertex_value().avg()/get_volume_vertex(path[i-1])->get_subpath_pdf())
						break;
					
					v0 = new Path<D>::VolumeVertex(*get_volume_vertex(path[i-1])); 
					v0->compute_scattering(path[i]->m_direction, f0, t0, p0);
					v0->m_f *= f0;
					v0->m_p *= p0;
					v0->m_pi *= p0;
				}//*/
				else if (get_surface_vertex(path[i-1]))
				{
					Path<D>::SurfaceVertex v(*get_surface_vertex(path[i-1])); 
					
					if (fabs(v.get_vertex_value().avg()/v.get_subpath_pdf()) == std::numeric_limits<float>::infinity() 
						|| v.get_vertex_value().avg()/v.get_subpath_pdf() != v.get_vertex_value().avg()/v.get_subpath_pdf())
						break;
					v.compute_scattering(path[i]->m_direction, f0, t0, p0);
					v.m_f *= f0;// * dot_clamped(-path[i-1]->m_direction, v.m_normal);
					v.m_p *= p0;
					v.m_pi *= p0;
					// TO-DO: Fix for transmissive materials, check if subpath is inside object and change medium
					v0 = new Path<D>::VolumeVertex(v.m_position, v.m_direction, m_world->get_medium(), v.get_vertex_value(), v.get_subpath_pdf(), v.get_vertex_pdf());
				}//*/
				
				if (v0) {
					// Change outgoing direction
					v0->m_direction = path[i]->m_direction; // in first vertex, it equal to light_sample.dir
					//Compute track length for short beam
					Real tracklength = MAX_LONGBEAM_DIST;

					Ray<D> r0(v0->m_position, v0->m_direction, false);
					Intersection<D> it0;
					m_world->first_intersection(r0, it0);

					v0->m_position = v0->m_position + v0->m_direction*m_longbeam_radius;

					Real safe_cap = it0.did_hit() ? (2 * this->m_longbeam_radius / dot_abs(-r0.get_direction(), it0.get_normal())) : 0;

					tracklength = min(r0.get_parameter(), tracklength) - safe_cap;
					global_particles.push_back(ParticleR(v0, i, ParticleR::LongBeam, m_longbeam_radius, tracklength, 0.0f, next_id++));
				}
			}
		}
		if( nb_shots_global && nb_shots_caustics )
			shooting = false;
	}

	Real scale = 1./static_cast<Real>(nb_shots_caustics);
	scale = 1./static_cast<Real>(nb_shots-1);
	
	typename std::vector<ParticleR>::iterator it;
	
	//Scale caustics' photons
	for( it = caustic_particles.begin(); it != caustic_particles.end(); it++ )
		it->m_power = (_finite(it->m_power.avg()*scale))?(it->m_power*scale):Radiance(0.);

	//scale = 1./static_cast<Real>(nb_shots_global);
	//Scale global' photons
	//global_photons.clear();
	for( it = global_particles.begin(); it != global_particles.end(); it++ )
		it->m_power = (_finite(it->m_power.avg()*scale))?(it->m_power*scale):Radiance(0.);
}
#else

template<int D>
void ParticleTracing<D>::trace_photons(std::vector<Particle<D> > &global_photons, std::vector<Particle<D> > &caustic_photons) const
{
	
	unsigned int nb_shots = 0;
	unsigned int nb_shots_caustics = 0;
	unsigned int nb_shots_global = 0;
	bool shooting = true;

	//Shoot photons!
	//Check if finished...
	while(shooting)
	{
		//Check if max number of shots done...
		if( ++nb_shots > m_max_nb_shots)
		{
			nb_shots_caustics = nb_shots_caustics?nb_shots_caustics:(nb_shots-1);
			nb_shots_global = nb_shots_global?nb_shots_global:(nb_shots-1);
			break;
		}
		
		// If not, start a new photon path
		// Sample light source, and from it sample the photon origin
		LightSample<D> photon_origin;
		Real pdf_light1, pdf_light2;
		this->m_world->sample_light(pdf_light1)->sample(photon_origin,pdf_light2);
		
		// Compute irradiance photon's energy
		Spectrum energy(photon_origin.irradiance/(pdf_light1*pdf_light2));
		
		Real time = 0; //picoseconds(10^-9)

		Ray<D> photon_ray(photon_origin.pos, photon_origin.dir);
		photon_ray.shift();

		bool is_caustic_particle = false;

		//Iterate the path
		while(1)
		{
			// Throw ray and update current_it
			Intersection<D> it;
			this->m_world->first_intersection(photon_ray, it);

			if( !it.did_hit() )
				break;

			//Update time
			time += this->m_world->time_of_flight(it.get_ray().get_parameter())*it.get_ray().get_refraction_index();

			//Check if has hit a delta material...
			if( it.material()->is_type(Reflectance::DELTA) )
			{
				// If delta material, then is caustic...
				// Don't store the photon!
				is_caustic_particle = true;
			}
			else if (photon_ray.get_level() > 0)
			{
				//If non-delta material, store the photon!
				if( is_caustic_particle )	//If caustic particle, store in caustics
				{				
					if( !nb_shots_caustics )
					{
						if( caustic_photons.size() == m_max_caustic_longbeams)
						{	nb_shots_caustics = nb_shots; break;	}

						caustic_photons.push_back( Particle<D>(it, energy, time, it.get_ray().get_level()) );
					}
				}
				else						//If non-caustic particle, store in global
				{
					if( !nb_shots_global ) 
					{
						if( global_photons.size() == m_nb_global_photons )
						{	nb_shots_global = nb_shots; break;	}

						global_photons.push_back( Particle<D>(it, energy, time, it.get_ray().get_level()) );
					}
				}
				is_caustic_particle = false;
			}	
			
			Real pdf;

			Spectrum surf_albedo = Spectrum(1.)-it.material()->get_absorption(it.get_uv());
			Real avg_surf_albedo = surf_albedo.avg();

			if( !m_use_russian_roulette )
				avg_surf_albedo = 1.;


			Real epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
			while (epsilon2 < 0.)
				epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
			
			if (epsilon2 > avg_surf_albedo || time > m_max_time ) 
				break;
				
			// Random walk's next step
			// Get sampled direction plus pdf, and update attenuation
			energy *= it.material()->sample_outgoing_ray(it, photon_ray, pdf );
			if( it.material()->is_type(Reflectance::DELTA) )
				energy *= dot_abs(it.get_normal(), photon_ray.get_direction()); 			
	
			energy /= (pdf*avg_surf_albedo);
		}
		
		if( nb_shots_global && nb_shots_caustics )
			shooting = false;
	}

	// Finalize photon shooting
	Real scale = 1./static_cast<Real>(nb_shots_caustics);
	scale = 1./static_cast<Real>(nb_shots);
	
	typename std::list<Particle<D> >::iterator it;
	
	//Scale caustics' photons
	for( it = caustic_photons.begin(); it != caustic_photons.end(); it++ )
		it->m_power *= scale;

	//scale = 1./static_cast<Real>(nb_shots_global);
	//Scale global' photons
	//global_photons.clear();
	for( it = global_photons.begin(); it != global_photons.end(); it++ )
		it->m_power *= scale;
}//*/
#endif
#endif //_PARTICLE_TRACING_H_

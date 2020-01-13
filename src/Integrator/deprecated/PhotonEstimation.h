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


#ifndef _PHOTON_ESTIMATION_H_
#define _PHOTON_ESTIMATION_H_

#include "bunnykiller.h"

#include "RayTracing/Intersection.h"
#include "RayTracing/Ray.h"
#include "Integrator/Integrator.h"
#include "Integrator/Particle.h"
#include "Integrator/ParticleTracing.h"
#include "Color/Spectrum.h"
#include "DataStructures/KDTree.h"
#include "DensityEstimation/VolumetricRadianceEstimators.h"
#include "DensityEstimation/DensityEstimationKernel.h"

#include <cmath>
#include <list>
#include <cstdio>

#ifndef _TEMPORAL_KERNEL_SAMPLES_ 
#define _TEMPORAL_KERNEL_SAMPLES_ 60
#endif

typedef DensityEsimationKernels::Box<1> Kernel1D;
typedef DensityEsimationKernels::Box<2> Kernel2D;

template<unsigned D, class Radiance, class RadianceAttenuation, class VREstimator>
class PhotonEstimation: public ParticleTracing<D, Radiance, RadianceAttenuation>
{
protected:
	using PTR = ParticleTracing<D, Radiance, RadianceAttenuation>;
	using RadianceSampleR = typename PTR::RadianceSampleR;
	using RadSampleList = typename PTR::RadSampleList;
	using FilmR = typename PTR::FilmR;
	using PathR = typename PTR::PathR;
	using VertexR = typename PTR::VertexR;
	using ParticleR = typename PTR::ParticleR;
protected:
	struct RadianceEstimationSample 
	{
		Spectrum radiance;
		Intersection<D> it;
		Real time;
		unsigned int bounce;
		RadianceEstimationSample(Spectrum _r, Intersection<D> _it, Real _time, unsigned int _b)
			:radiance(_r), it(_it), time(_time), bounce(_b){};
	} ;
	DensityEsimationKernels::Box<D-1> K;
	Kernel1D Kr;
	Kernel2D Kt;
	VREstimator m_vre;
	
	KDTree<ParticleR, D> m_global_spoints_tree, m_caustic_spoints_tree;
	KDTree<ParticleR, D> m_global_vpoints_tree, m_caustic_vpoints_tree;
	
	std::vector<ParticleR> m_global_longbeams, m_caustic_longbeams;
	std::vector<ParticleR> m_global_shortbeams, m_caustic_shortbeams;

	Real m_spoint_kernel_norm; 
	Real m_vpoint_kernel_norm; 
	Real m_time_bandwidth, m_temporal_kernel_norm;
	bool m_spatial_only_kernel;

	unsigned int m_nb_spoints, m_nb_samples_ray;
	bool m_final_gathering;
	DirectIntegrator<D> *direct;
	
	Scattering::Level m_scattering_level;
	World<D, Radiance> *m_w;
	typedef typename KDTree<ParticleR, D>::Node KDTreeNode;
public:
	PhotonEstimation(World<D> *w, 
					ParticleTracing<D>::TracingSetup tracing_setup, 
					int search_knn_spoints,
					int samples_ray = 10,
					bool final_gathering = false, 
					Film *film = NULL,
					VREstimator vre = VREstimator())
					:ParticleTracing(w, tracing_setup, film), m_w(w), direct(new DirectIntegrator<D>(w)), m_nb_spoints(search_knn_spoints),  
					m_final_gathering(final_gathering), 
					m_nb_samples_ray(samples_ray),
					m_spatial_only_kernel(true), 
					m_time_bandwidth(0.),
					m_scattering_level(tracing_setup.scatt_level),
					m_vre(vre)
	{}

	~PhotonEstimation(){}

	// Functions to set irradiance estimation parameters
	void set_max_dist(const Real max_dist){ m_spoint_radius = max_dist; m_spoint_kernel_norm = 1./(M_PI*m_spoint_radius*m_spoint_radius);}
	void set_nb_photons(const int nb_photons){ m_nb_spoints = nb_photons;}
	
	void set_temporal_kernel(bool on){m_spatial_only_kernel=!on;}
	void set_scattering_level(Scattering::Level level){ m_scattering_level = level; }

	virtual void preprocess();
	virtual void store(const std::string &namefile)const;
	virtual void load(const std::string &namefile)const;
	
	virtual Spectrum operator()(const Ray<D> &r);
	virtual Spectrum operator()(const Intersection<D> &it) const;
	
	virtual void operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
	virtual void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;

}; //PhotonEstimation

// Preprocess functions
template<unsigned D, class VREstimator>
void PhotonEstimation<D, VREstimator>::preprocess()
{
	unsigned int nb_spoints, nb_vpoints, nb_longbeams, nb_shortbeams, nb_volparticles;
	nb_spoints = nb_vpoints = nb_longbeams = nb_shortbeams = nb_volparticles = 0;

	//Clear trees
	m_global_spoints_tree.clear();
	m_caustic_spoints_tree.clear();
	m_global_vpoints_tree.clear();
	m_caustic_vpoints_tree.clear();

	//Trace photons
	std::vector<Particle<D> > global_particles, caustic_particles;
	
	this->trace_photons(global_particles, caustic_particles);

	//Store photons into KDTrees and balance KDTrees
	typename std::vector<Particle<D> >::iterator it;
	//global_particles.clear();

	printf("Total Particles: %7d global particles, %7d caustic particles\n", global_particles.size(), caustic_particles.size());

	
	FILE *f = fopen("global_particles.txt", "w+");
	
	if( global_particles.size() )
	{
		std::vector<Particle<D> > vol_particles;
		for( it = global_particles.begin(); it != global_particles.end(); it ++)
		{	
			switch (it->m_type){
			case Particle<D>::SurfPoint:
				nb_spoints++;
				m_global_spoints_tree.store(std::vector<Real>(it->m_v->get_vertex_position().m_data, 
										it->m_v->get_vertex_position().m_data+D), *it);
				break;
			default:
				//fprintf(f, "%f %f %f %f %f %f %f\n", it->m_v->m_position[0], it->m_v->m_position[1], it->m_v->m_position[2], it->m_v->m_direction[0], it->m_v->m_direction[1], it->m_v->m_direction[2], it->m_v->m_f.avg());
				nb_volparticles++;
				vol_particles.push_back(*it);
				break;
			}
		}
		m_vre = VREstimator(vol_particles, m_nb_samples_ray);
		m_vre.m_time_exp = film->get_exposure_time();
		m_vre.m_world = m_world;
		printf("Total surface points %7d\n", nb_spoints);
		printf("Total media particles: %7d\n", nb_volparticles);

		if (nb_spoints) 
			m_global_spoints_tree.balance();
		
		global_particles.clear();
	}
	
		fclose(f);
	/*if( caustic_particles.size() )
	{
		for( it = caustic_particles.begin(); it != caustic_particles.end(); it ++)
		{	
			switch (it->m_type){
			case Particle<D>::SurfPoint:
				m_caustic_spoints_tree.store(std::vector<Real>(it->m_v->get_vertex_position().m_data, 
										it->m_v->get_vertex_position().m_data+D), *it);
				break;
			case Particle<D>::VolPoint:
				break;
			case Particle<D>::ShortBeam:
				break;
			case Particle<D>::LongBeam:
				m_global_longbeams.push_back(*it);
				break;
			}
		}
		caustic_particles.clear();
		m_caustic_spoints_tree.balance();
	}*/
}

template<unsigned D, class VREstimator>
void PhotonEstimation<D, VREstimator>::store(const std::string &namefile)const
{
	
	if( m_global_spoints_tree.is_empty() && m_caustic_spoints_tree.is_empty() )
		return;

	FILE* f = fopen(namefile.c_str(), "wb");
	int nb_globals = m_global_spoints_tree.nb_elements(), 
		nb_caustics = m_caustic_spoints_tree.nb_elements();

	// Store dimensions
	int dimensions[3] = {D, sizeof(Spectrum)/sizeof(Real), sizeof(Real)};
	fwrite(&dimensions, sizeof(int),3,f);
	// Store global photons
	fwrite(&nb_globals, sizeof(int),1,f);
	for( unsigned int i=0; i<nb_globals; ++i )
	{
		typename KDTree<Particle<D>, D>::Node n = m_global_spoints_tree[i];
		
		fwrite((void*)(&n.data().m_v->m_position),sizeof(VectorN<D>),1,f);
		fwrite((void*)(&n.data().m_power),sizeof(Spectrum),1,f);
	}

	// Store caustics photons
	fwrite(&nb_caustics, sizeof(int),1,f);
	for( unsigned int i=0; i<nb_caustics; ++i )
	{
		typename KDTree<Particle<D>, D>::Node n = m_caustic_spoints_tree[i];
		
		fwrite((void*)(&n.data().m_v->m_position),sizeof(VectorN<D>),1,f);
		fwrite((void*)(&n.data().m_power),sizeof(Spectrum),1,f);
	}

	fclose(f);
}

template<int D, class VREstimator>
void PhotonEstimation<D, VREstimator>::load(const std::string &namefile)const
{
	// TODO
}

//*****************************************************************************
//=============================================================================
// Integration functions

// Steady
template<int D, class VREstimator>
Spectrum PhotonEstimation<D, VREstimator>::operator()(const Ray<D> &r)
{
	Ray<D> r0 (r);
	Intersection<D> it;
	Spectrum f_s(0.), f_v(0.), f(0.), Tr_s(1.);

	int its = 0;
	bool run = true;

	while( run ) 
	{
		f_s = 0;
		if (m_w->first_intersection(r0,it) && its < 5)
		{
			/*if( !it.material()->is_type(Reflectance::DELTA) )
				f_s = (*this)(it);*/
		}
		else
			run = false;
		
		++its;

		Medium<D> *m = r0.get_medium();
		if (m)
		{
			//VolPoints
			f_v = m_vre(r0)*Tr_s;

			//Surface
			if (it.did_hit()) 
				Tr_s *= m->get_transmittance(r0);
		}
		f += f_s*Tr_s + f_v;


		if( it.did_hit() && it.material()->is_type(Reflectance::DELTA) )
		{
			Real pdf = 1.;
			Tr_s *= it.material()->sample_outgoing_ray( it, r0, pdf );
			Tr_s /= pdf;
		}				
	}
	return f;
}

// Transient
template<unsigned D, class VREstimator>
void PhotonEstimation<D, VREstimator>::operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const
{
	Ray<D> r0 (r);
	Intersection<D> it;
	Spectrum Tr_s(1.);
	std::list<RadianceSample> surface_samples, volume_samples;
	std::list<RadianceSample>::iterator rsit;
	
	Real time(0);
	int its = 0;
	bool run = true;

	while( run ) 
	{
		if( m_w->first_intersection(r0,it) && its < 5 )
		{
			if( !it.material()->is_type(Reflectance::DELTA) )
				(*this)(it, surface_samples);
		}
		else
			run = false;
		
		++its;

		Medium<D> *m = r0.get_medium();
		if (m)
		{
			//VolPoints
			m_vre(r0, volume_samples);
			
			for (rsit=volume_samples.begin(); rsit != volume_samples.end(); rsit++)
			{	
				rsit->radiance*=Tr_s; 
				rsit->time+=time;
				samples.push_back(*rsit);
			}
			volume_samples.clear();

			//Surface
			if (it.did_hit()) 
				Tr_s *= m->get_transmittance(r0);
		}

		for (rsit=surface_samples.begin(); rsit != surface_samples.end(); rsit++)
		{	
			rsit->radiance*=Tr_s; 
			rsit->time+=time;
			samples.push_back(*rsit);
		}

		surface_samples.clear();
		
		
		if( it.did_hit() && it.material()->is_type(Reflectance::DELTA) )
		{
			time += this->m_world->time_of_flight(r0.get_parameter())*r0.get_refraction_index();
			
			Real pdf = 1.;
			Tr_s *= it.material()->sample_outgoing_ray( it, r0, pdf );
			Tr_s /= pdf;
		}	
	}
}


template<unsigned D, class VREstimator>
Spectrum PhotonEstimation<D, VREstimator>::operator()(const Intersection<D> &it) const
{
	Spectrum reflected_radiance(0.), attenuation(1.);
	//printf("Mielda1...\n");

	//Handle delta materials (e.g. perfect reflection/refraction)
	Intersection<D> shaded_it(it);
	Ray<D> new_ray;
	int iter=0;

	while( shaded_it.material()->is_type(Reflectance::DELTA) )
	{
		Real pdf2(1.);

		//Sample outgoing ray
		attenuation *= shaded_it.material()->sample_outgoing_ray(shaded_it, new_ray, pdf2 );
		//attenuation *= dot_abs(shaded_it.get_normal(), new_ray.get_direction()); 			
		attenuation /= (pdf2);

		//Get new intersection
		this->m_world->first_intersection(new_ray, shaded_it);

		if( !shaded_it.did_hit() || iter > 30 )
			return Spectrum(0.);

		++iter;
	}

	//	Get all samples where the radiance estimation is needed.
	//	If we are performing final gathering, compute the reflected
	//	samples. Else, just add the current intersection.
	list<pair<Spectrum, Intersection<D> > > samples;
	samples.push_back(pair<Spectrum,Intersection<D> >(Spectrum(1.),shaded_it));

	
	//Compute radiance estimation! 
	if( !(m_scattering_level & Scattering::MULTIPLE) )
		return (*direct)(shaded_it)*attenuation;

	typename std::list<pair< Spectrum,Intersection<D> > >::iterator cit;
	for( cit = samples.begin(); cit != samples.end(); cit++ )
	{
		Spectrum sample_radiance(0.);
		VectorN<D> p = cit->second.get_position();
		VectorN<D> d = -cit->second.get_ray().get_direction();

		typename KDTree<Particle<D>,D>::Node;
	
#ifdef NN
		std::vector<const KDTreeNode*> photons;
		typename std::vector<const KDTreeNode*>::iterator ph;
		Real max_distance;

		// Get global contribution
		//	Find closest global photons
		m_global_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_spoints,  photons, max_distance);
		//	Compute their contribution
		for( ph = photons.begin(); ph != photons.end(); ph++ ) 
			sample_radiance += cit->second.material()->f( (*ph)->data().m_direction, d, cit->second.get_normal(), cit->second.get_uv())
											*((*ph)->data().m_power)
											* dot_clamped(-(*ph)->data().m_direction,cit->second.get_normal());

		//	Finally, add it to the reflected radiance, 
		//	multiplied by the density estimator
		reflected_radiance += sample_radiance*cit->first*(1./(M_PI*max_distance*max_distance));

		sample_radiance = Spectrum();
		// Get caustics contribution
		//	Find closest caustic photons
		photons.clear();
		m_caustic_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_spoints,  photons, max_distance);
		//	Compute their contribution
		for( ph = photons.begin(); ph != photons.end(); ph++ )
			sample_radiance += cit->second.material()->f( (*ph)->data().m_direction, d, cit->second.get_normal(), cit->second.get_uv())
											*((*ph)->data().m_power)
											* dot_clamped(-(*ph)->data().m_direction,cit->second.get_normal());


		//	Finally, add it to the reflected radiance, 
		//	multiplied by the density estimator
		reflected_radiance += sample_radiance*cit->first*(1./(M_PI*max_distance*max_distance));
#else
		std::vector<const KDTreeNode*> photons;
		typename std::vector<const KDTreeNode*>::iterator ph;
		int photons_found; 
		//-----------------------------------------------------------------------------
		// Get global contribution
		if( !m_global_spoints_tree.is_empty() )
		{
			//	Find closest global photons
			Real max_dist = m_spoint_radius;
			m_global_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_spoints, photons, max_dist);
			Spectrum f;
			//	Compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ ) 
			{
				(*ph)->data().m_v->compute_scattering(d, f);
				sample_radiance += f* K(((*ph)->data().m_v->get_vertex_position() - p).length(), max_dist)
									* (*ph)->data().m_power;
			}
			//	Finally, add it to the reflected radiance, 
			//	multiplied by the density estimator
			reflected_radiance += sample_radiance*cit->first;
		}
		
		//-----------------------------------------------------------------------------
		// Get caustics contribution
		if( !m_caustic_spoints_tree.is_empty() )
		{
			sample_radiance = Spectrum();
			
			//	Find closest caustic photons
			photons.clear();
			//photons_found = m_caustic_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_spoint_radius, &photons);
			Real max_dist = m_spoint_radius;
			m_global_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_spoints, photons, max_dist);
			Spectrum f;
			//	Compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ )
			{	
				(*ph)->data().m_v->compute_scattering(d, f);
				sample_radiance += f* K(((*ph)->data().m_v->get_vertex_position() - p).length(), max_dist)
									*(*ph)->data().m_power;
				/*sample_radiance += cit->second.material()->f( (*ph)->data().m_direction, d, cit->second.get_normal(), cit->second.get_uv())
												*((*ph)->data().m_power)
												* dot_clamped(-(*ph)->data().m_direction,cit->second.get_normal())
												* Kr(((*ph)->data().m_position - p).length(), m_spoint_radius);
				*/
			}
			//	Finally, add it to the reflected radiance, 
			//	multiplied by the density estimator
			reflected_radiance += sample_radiance*cit->first;
		}

#endif

	}

	if( !(m_scattering_level & Scattering::SINGLE) )
		return reflected_radiance;

	return (reflected_radiance + (*direct)(shaded_it))*attenuation;
}


//=============================================================================
template<unsigned D, class VREstimator>
void PhotonEstimation<D, VREstimator>::operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples =1) const
{
	// Workaround for RadianceSensor in steady state!! This must be changed asap by a proper transient version
	//samples.push_back(RadianceSample(m_vre(p), 0., 0., 0.));
}

//=============================================================================
template<unsigned D, class VREstimator>
void PhotonEstimation<D, VREstimator>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
	//printf("Entering PM integrator\n");
	
	//-----------------------------------------------------------------------------
	//Handle delta materials (e.g. perfect reflection/refraction)
	Spectrum attenuation(1.);
	Intersection<D> shaded_it(it);
	Ray<D> new_ray;
	int iter=0;

	Real time = m_world->time_of_flight(it.get_ray().get_parameter());

	//printf("Starting random walk\n");
	while( shaded_it.material()->is_type(Reflectance::DELTA) )
	{
		Real pdf2(1.);
		Real delta_time;

		//Sample outgoing ray
		attenuation *= shaded_it.material()->sample_outgoing_ray(shaded_it, new_ray, delta_time, pdf2 );
		//attenuation *= dot_abs(shaded_it.get_normal(), new_ray.get_direction()); 			
		attenuation /= (pdf2);

		//Get new intersection
		this->m_world->first_intersection(new_ray, shaded_it);

		if( !shaded_it.did_hit() || iter > 30 )
			return;

		time += delta_time+this->m_world->time_of_flight(new_ray.get_parameter())*new_ray.get_refraction_index();
		++iter;
	}

	//printf("Finished RW\n");

	std::list<RadianceEstimationSample> estimation_samples;
	estimation_samples.push_back(RadianceEstimationSample(Spectrum(1.),shaded_it, time, shaded_it.get_ray().get_level()));

	//-----------------------------------------------------------------------------
	//Compute radiance estimation! 
	typename std::list<RadianceEstimationSample>::iterator cit;
	for( cit = estimation_samples.begin(); cit != estimation_samples.end(); cit++ )
	{
		
		if( (m_scattering_level & Scattering::MULTIPLE) )
		{
		//printf("Starting Radiance Estimate\n");

		VectorN<D> p = cit->it.get_position();
		VectorN<D> d = -cit->it.get_ray().get_direction();

		// Initialize density estimation variables...
		std::list<const KDTreeNode*> photons;
		typename std::list<const KDTreeNode*>::iterator ph;
		
		Spectrum fc_global = attenuation; // Here should go the new kernel...

		
		// 
		int nb_iterations = _TEMPORAL_KERNEL_SAMPLES_;//4* //static_cast<int>(m_time_bandwidth/ + 1.);
		Real kernel_delta_time = 2*m_time_bandwidth/static_cast<Real>(nb_iterations-1);
		Real time_marching_norm = m_spatial_only_kernel?1.:(2*m_time_bandwidth/static_cast<Real>(nb_iterations));
		int photons_found;
		
		//-----------------------------------------------------------------------------
		// Get global contribution
		if( !m_global_spoints_tree.is_empty() )
		{
			// - Find closest global photons... (fixed radius)
			photons_found = m_global_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_spoint_radius, &photons);

			// - ... and compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ )
			{
				// Spatial kernel
				Real delta_time, pdf_time;
				Spectrum f;
				(*ph)->data().m_v->compute_scattering(d, f, delta_time, pdf_time);
				Spectrum f_p = f* K(((*ph)->data().m_v->get_vertex_position() - p).length(), m_spoint_radius)
									*(*ph)->data().m_power 
									* attenuation
									* time_marching_norm / pdf_time;
				
				Real total_time = delta_time + cit->time + (*ph)->data().m_v->m_time;
				//int total_level = cit->bounce + (*ph)->data().m_level;

				if( m_spatial_only_kernel )
					samples.push_back( RadianceSample( f_p, total_time, 0, cit->time));
				else
					for( Real current_time = -m_time_bandwidth; current_time < m_time_bandwidth; current_time += kernel_delta_time )
						samples.push_back( RadianceSample( f_p*Kt(fabs(current_time),m_time_bandwidth), 
	 														total_time + current_time));
			}
		}

		//printf("Look at caustics photons\n");
		//-----------------------------------------------------------------------------
		// Get caustics contribution
		
		/*if( !m_caustic_spoints_tree.is_empty() )
		{
			// - Find closest caustic photons...
			photons.clear();
			photons_found = m_caustic_spoints_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_spoint_radius, &photons);
			//printf("-%d photons found!\n", photons_found );
			// - ... and compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ )
			{
				// Spatial kernel
				Real delta_time, pdf_time;
				Spectrum f_p = cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv(),
													delta_time, pdf_time)
							 * ((*ph)->data().m_power) * fc_global 
							 * ( dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal()) 
							 * Kr(((*ph)->data().m_position - p).length(), m_spoint_radius)
							 * time_marching_norm / pdf_time); 
				
				Real total_time = delta_time + cit->time + (*ph)->data().m_time;
				int total_level = cit->bounce + (*ph)->data().m_level;

				//printf("--spatial kernel done\n");
				// Temporal kernel
				if( m_spatial_only_kernel )
					samples.push_back( RadianceSample( f_p, total_time, total_level)); // Here should go the new kernel...
				else
					for( Real current_time = -m_time_bandwidth; current_time < m_time_bandwidth; current_time += kernel_delta_time )
						samples.push_back( RadianceSample( f_p*Kt(fabs(current_time),m_time_bandwidth), 
															total_time + current_time, total_level)); // Here should go the new kernel...
			}
		}*/
		}
		if( (m_scattering_level & Scattering::SINGLE) )
		{
			// Add the (attenuated) direct samples
			std::list<RadianceSample> direct_samples;
			std::list<RadianceSample>::iterator ds_it;
			(*direct)(cit->it,direct_samples);
			for(ds_it = direct_samples.begin(); ds_it != direct_samples.end(); ds_it++)
				samples.push_back( RadianceSample( ds_it->radiance*attenuation, ds_it->time + cit->time, ds_it->bounce + cit->bounce ));
		}
	}
}

#endif // _PHOTON_ESTIMATION_H_

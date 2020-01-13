#ifndef _PHOTON_MAPPING_H_
#define _PHOTON_MAPPING_H_

#include "Photon.h"
#include "ParticleTracing.h"
#include "RayTracing/Intersection.h"
#include "RayTracing/Ray.h"
#include "Color/Spectrum.h"
#include "DataStructures/KDTree.h"
#include "DensityEstimationKernel.h"
#include "Integrator/VolumeIntegrator.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <list>
#include <stdio.h>

#ifndef _TEMPORAL_KERNEL_SAMPLES_ 
#define _TEMPORAL_KERNEL_SAMPLES_ 60
#endif 

typedef DensityEsimationKernels::Perlin<1> Kernel1D;
typedef DensityEsimationKernels::Perlin<2> Kernel2D;

template<class TermCriteria, int D>
class PhotonMapping: public ParticleTracing<TermCriteria,D>
{
	Kernel1D Kt;
	Kernel2D Kr;

	KDTree<Photon<D>, D> m_global_tree, m_caustics_tree;
	Real m_radius, m_spatial_kernel_norm; 
	Real m_time_bandwidth, m_temporal_kernel_norm;
	bool m_spatial_only_kernel;

	int m_nb_photons;
	bool m_final_gathering;
	DirectIntegrator<D> *direct;
	
	Scattering::Level m_scattering_level;
	
	typedef typename KDTree<Photon<D>, D>::Node KDTreeNode;
public:
	PhotonMapping(World<D> *w, Real radius, int nb_global_photons, int nb_caustics_photons,
				  int max_nb_global_photons = 5000000, unsigned int incoming_samples = 1 )
				  :ParticleTracing<TermCriteria,D>(w,incoming_samples,nb_global_photons,nb_caustics_photons,max_nb_global_photons),
					m_radius(radius), m_spatial_kernel_norm((D==2)?(1./m_radius):(1./(M_PI*m_radius*m_radius))), 
					m_nb_photons(400), m_final_gathering(false),m_scattering_level(Scattering::ALL),
					direct(new DirectIntegrator<D>(w)), m_time_bandwidth(0.), m_temporal_kernel_norm(0), m_spatial_only_kernel(true)
	{}
	PhotonMapping(World<D> *w, Real radius, Real time_bandwidth, int nb_global_photons, int nb_caustics_photons,
				  int max_nb_global_photons = 5000000, unsigned int incoming_samples = 1 )
				  :ParticleTracing<TermCriteria,D>(w,incoming_samples,nb_global_photons,nb_caustics_photons,max_nb_global_photons),
					m_radius(radius), m_spatial_kernel_norm((D==2)?(1./m_radius):(1./(M_PI*m_radius*m_radius))), 
					m_nb_photons(400), m_final_gathering(false),
					direct(new DirectIntegrator<D>(w)), m_time_bandwidth(time_bandwidth/2), m_temporal_kernel_norm(1./time_bandwidth),
					m_spatial_only_kernel(false),m_scattering_level(Scattering::ALL)
	{}

	PhotonMapping(World<D> *w, Real radius, int nb_photons, bool final_gathering,
				  int nb_global_photons = std::numeric_limits<int>::max(),
				  int nb_caustics_photons = std::numeric_limits<int>::max(),
				  int max_nb_global_photons = 5000000,
				  unsigned int incoming_samples = 1)
				  :ParticleTracing<TermCriteria,D>(w,incoming_samples,nb_global_photons,nb_caustics_photons,max_nb_global_photons),
				  m_radius(radius), m_spatial_kernel_norm((D==2)?(1./m_radius):(1./(M_PI*m_radius*m_radius))), m_nb_photons(nb_photons), m_final_gathering(final_gathering),
					direct(new DirectIntegrator<D>(w)),m_scattering_level(Scattering::ALL)
	{}
	PhotonMapping(World<D> *w, Real radius, int nb_photons,
				  int nb_global_photons = std::numeric_limits<int>::max(),
				  int nb_caustics_photons = std::numeric_limits<int>::max(),
				  int max_nb_global_photons = 5000000,
				  unsigned int incoming_samples = 1)
				  :ParticleTracing<TermCriteria,D>(w,incoming_samples,nb_global_photons,nb_caustics_photons,max_nb_global_photons),
				  m_radius(radius),m_spatial_kernel_norm((D==2)?(1./m_radius):(1./(M_PI*m_radius*m_radius))), m_nb_photons(nb_photons), m_final_gathering(false),
				  direct(new DirectIntegrator<D>(w)),m_scattering_level(Scattering::ALL)
	{}
	~PhotonMapping(){}

	// Functions to set irradiance estimation parameters
	void set_max_dist(const Real max_dist){ m_radius = max_dist; m_spatial_kernel_norm = 1./(M_PI*m_radius*m_radius);}
	void set_nb_photons(const int nb_photons){ m_nb_photons = nb_photons;}
	
	void set_temporal_kernel(bool on){m_spatial_only_kernel=!on;}
	void set_scattering_level(Scattering::Level level){ m_scattering_level = level; }

	virtual void preprocess();
	virtual void store(const std::string &namefile)const;
	virtual void load(const std::string &namefile)const;

	virtual Spectrum operator()(const Intersection<D> &it) const;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
}; //PhotonMapping


// Preprocess functions
template<class TermCriteria, int D>
void PhotonMapping<TermCriteria,D>::preprocess()
{
	//Clear trees
	m_global_tree.clear();
	m_caustics_tree.clear();

	//Trace photons
	std::list<Particle<D> > global_photons, caustic_photons;
	this->trace_photons(global_photons, caustic_photons);

	//Store photons into KDTrees and balance KDTrees
	typename std::list<Particle<D> >::iterator it;
	//global_photons.clear();

	printf("Total Photons: %7d global photons, %7d caustic photons\n", global_photons.size(), caustic_photons.size());

	if( global_photons.size() )
	{
		for( it = global_photons.begin(); it != global_photons.end(); it ++)
			m_global_tree.store(std::vector<Real>(it->m_hit_point.get_position().m_data, 
												  it->m_hit_point.get_position().m_data+D), Photon<D>(*it));
		global_photons.clear();
		m_global_tree.balance();
	}
	if( caustic_photons.size() )
	{
		for( it = caustic_photons.begin(); it != caustic_photons.end(); it ++)
			m_caustics_tree.store(std::vector<Real>(it->m_hit_point.get_position().m_data, 
													it->m_hit_point.get_position().m_data+D), Photon<D>(*it));
		
		caustic_photons.clear();
		m_caustics_tree.balance();
	}

}

template<class TermCriteria, int D>
void PhotonMapping<TermCriteria,D>::store(const std::string &namefile)const
{
	if( m_global_tree.is_empty() && m_caustics_tree.is_empty() )
		return;

	FILE* f = fopen(namefile.c_str(), "wb");
	int nb_globals = m_global_tree.nb_elements(), 
		nb_caustics = m_caustics_tree.nb_elements();

	// Store dimensions
	int dimensions[3] = {D, sizeof(Spectrum)/sizeof(Real), sizeof(Real)};
	fwrite(&dimensions, sizeof(int),3,f);
	// Store global photons
	fwrite(&nb_globals, sizeof(int),1,f);
	for( unsigned int i=0; i<nb_globals; ++i )
	{
		typename KDTree<Photon<D>, D>::Node n = m_global_tree[i];
		
		fwrite((void*)(&n.data().m_position),sizeof(VectorN<D>),1,f);
		fwrite((void*)(&n.data().m_power),sizeof(Spectrum),1,f);
	}

	// Store caustics photons
	fwrite(&nb_caustics, sizeof(int),1,f);
	for( unsigned int i=0; i<nb_caustics; ++i )
	{
		typename KDTree<Photon<D>, D>::Node n = m_caustics_tree[i];
		
		fwrite((void*)(&n.data().m_position),sizeof(VectorN<D>),1,f);
		fwrite((void*)(&n.data().m_power),sizeof(Spectrum),1,f);
	}

	fclose(f);
}

template<class TermCriteria, int D>
void PhotonMapping<TermCriteria,D>::load(const std::string &namefile)const
{
	// TODO
}

//*****************************************************************************
//=============================================================================
// Integration functions
template<class TermCriteria, int D>
Spectrum PhotonMapping<TermCriteria,D>::operator()(const Intersection<D> &it) const
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

		typename KDTree<Photon<D>,D>::Node;
		
#ifdef NN
		std::vector<const KDTreeNode*> photons;
		typename std::vector<const KDTreeNode*>::iterator ph;
		Real max_distance;

		// Get global contribution
		//	Find closest global photons
		m_global_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_photons,  photons, max_distance);
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
		m_caustics_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_photons,  photons, max_distance);
		//	Compute their contribution
		for( ph = photons.begin(); ph != photons.end(); ph++ )
			sample_radiance += cit->second.material()->f( (*ph)->data().m_direction, d, cit->second.get_normal(), cit->second.get_uv())
											*((*ph)->data().m_power)
											* dot_clamped(-(*ph)->data().m_direction,cit->second.get_normal());

		//	Finally, add it to the reflected radiance, 
		//	multiplied by the density estimator
		reflected_radiance += sample_radiance*cit->first*(1./(M_PI*max_distance*max_distance));
#else
		std::list<const KDTreeNode*> photons;
		typename std::list<const KDTreeNode*>::iterator ph;
		int photons_found; 
		//-----------------------------------------------------------------------------
		// Get global contribution
		if( !m_global_tree.is_empty() )
		{
			//	Find closest global photons
			photons_found = m_global_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_radius, &photons);
		
			//	Compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ ) 
				sample_radiance += cit->second.material()->f( (*ph)->data().m_direction, d, cit->second.get_normal(), cit->second.get_uv())
												*((*ph)->data().m_power)
												* dot_clamped(-(*ph)->data().m_direction,cit->second.get_normal())
												* Kr(((*ph)->data().m_position - p).length(), m_radius);

			//	Finally, add it to the reflected radiance, 
			//	multiplied by the density estimator
			reflected_radiance += sample_radiance*cit->first;
		}
		
		//-----------------------------------------------------------------------------
		// Get caustics contribution
		if( !m_caustics_tree.is_empty() )
		{
			sample_radiance = Spectrum();
			
			//	Find closest caustic photons
			photons.clear();
			photons_found = m_caustics_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_radius, &photons);
			
			//	Compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ )
				sample_radiance += cit->second.material()->f( (*ph)->data().m_direction, d, cit->second.get_normal(), cit->second.get_uv())
												*((*ph)->data().m_power)
												* dot_clamped(-(*ph)->data().m_direction,cit->second.get_normal())
												* Kr(((*ph)->data().m_position - p).length(), m_radius);


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
template<class TermCriteria, int D>
void PhotonMapping<TermCriteria,D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
	//printf("Entering PM integrator\n");
	
	//-----------------------------------------------------------------------------
	//Handle delta materials (e.g. perfect reflection/refraction)
	Spectrum attenuation(1.);
	Intersection<D> shaded_it(it);
	Ray<D> new_ray;
	int iter=0;

	Real time = 0;

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

	std::list<RadianceEstimationSample<D>> estimation_samples;
	estimation_samples.push_back(RadianceEstimationSample<D>(Spectrum(1.),shaded_it, time, shaded_it.get_ray().get_level()));

	//-----------------------------------------------------------------------------
	//Compute radiance estimation! 
	typename std::list<RadianceEstimationSample<D>>::iterator cit;
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
		
		int nb_iterations = _TEMPORAL_KERNEL_SAMPLES_;//4*//static_cast<int>(m_time_bandwidth/ + 1.);
		Real kernel_delta_time = 2*m_time_bandwidth/static_cast<Real>(nb_iterations-1);
		Real time_marching_norm = m_spatial_only_kernel?1.:(2*m_time_bandwidth/static_cast<Real>(nb_iterations));
		int photons_found;
		
		//-----------------------------------------------------------------------------
		// Get global contribution
		if( !m_global_tree.is_empty() )
		{
			// - Find closest global photons... (fixed radius)
			photons_found = m_global_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_radius, &photons);

			// - ... and compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ )
			{
				// Spatial kernel
				Real delta_time, pdf_time;
				Spectrum f_p = cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv(),
													delta_time, pdf_time)
							 * ((*ph)->data().m_power) * fc_global 
							 * ( dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal())
							 * Kr(((*ph)->data().m_position - p).length(), m_radius)
							 * time_marching_norm / pdf_time);
				
				Real total_time = delta_time + cit->time + (*ph)->data().m_time;
				int total_level = cit->bounce + (*ph)->data().m_level;

				if( m_spatial_only_kernel )
					samples.push_back( RadianceSample( f_p, total_time, total_level));
				else
					for( Real current_time = -m_time_bandwidth; current_time < m_time_bandwidth; current_time += kernel_delta_time )
						samples.push_back( RadianceSample( f_p*Kt(fabs(current_time),m_time_bandwidth), 
	 														total_time + current_time, total_level));
			}
		}

		//printf("Look at caustics photons\n");
		//-----------------------------------------------------------------------------
		// Get caustics contribution
		if( !m_caustics_tree.is_empty() )
		{
			// - Find closest caustic photons...
			photons.clear();
			photons_found = m_caustics_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_radius, &photons);
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
							 * Kr(((*ph)->data().m_position - p).length(), m_radius)
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
		}
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

#endif //_PHOTON_MAPPING_H_

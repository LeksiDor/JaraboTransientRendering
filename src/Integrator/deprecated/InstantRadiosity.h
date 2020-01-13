#ifndef _INSTANT_RADIOSITY_H_
#define _INSTANT_RADIOSITY_H_

#include "ParticleTracing.h"
#include "DirectIntegrator.h"

#include "CosineLightSource/LightSource.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <list>


template<class TermCriteria, int D>
class InstantRadiosity: public ParticleTracing<TermCriteria,D>
{
	Real m_clamp_value;
	std::vector< CosineLightSource<D> > m_vpls;
	DirectIntegrator<D> m_direct_light_integrator;

public:
	InstantRadiosity(World<D> *w, 
				  int nb_vpls = std::numeric_limits<int>::max(),
				  Real clamp_value = .1,
				  int max_nb_shots = 5000000,
				  unsigned int incoming_samples = 1)
				  :ParticleTracing(w,incoming_samples, nb_vpls, 0, max_nb_shots), m_clamp_value(clamp_value),
				  m_direct_light_integrator(w)
	{}
	~InstantRadiosity(){}

	// Functions to set clamping parameters
	void set_clamp_value( Real clamp_value ){m_clamp_value = clamp_value; }

	virtual void preprocess();
	virtual Spectrum operator()(const Intersection<D> &it) const;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
}; //InstantRadiosity

template<class TermCriteria, int D>
void InstantRadiosity<TermCriteria,D>::preprocess()
{
	//Clear vpls
	m_vpls.clear();

	//Trace photons
	std::list< Particle<D> > global_photons, caustic_photons;
	trace_photons(global_photons, caustic_photons);

	//Store photons into KDTrees and balance KDTrees
	typename std::list<Particle<D> >::iterator it;
	//global_photons.clear();

	printf("Total VPLs: %7d vpls", global_photons.size());

	if( global_photons.size() )
	{
		m_vpls.reserve( global_photons.size() );

		for( it = global_photons.begin(); it != global_photons.end(); it ++)
			m_vpls.push_back(CosineLightSource(m_world, it->m_hit_point.get_position(), 
				it->m_hit_point.get_normal(), 
				it->m_power*it->m_hit_point.f(it->m_hit_point.get_ray().get_direction(),it->m_hit_point.get_normal()), 
				/*	We assume that the surface is lambertian, so its contribution will be the light outgoing on the normal.
					This should be changed in the future. Two different avenues: (1)getting the 'rho' from the material,
					and using it as the color of the surface, or (2) create a new type of VPL, that stores the normal,
					the material and the incoming direction. The latter will make Lightcuts more complex (a new bounding value 
					for the light will be required), but it'd allow much more flexibility on the materials supported 
					by VPL based methods. */
				it->m_time));
		
		global_photons.clear();
	}
}
template<class TermCriteria, int D>
Spectrum InstantRadiosity<TermCriteria,D>::operator()(const Intersection<D> &it) const
{
	Spectrum reflected_radiance(0.), attenuation(1.);

	//Handle delta materials (e.g. perfect reflection/refraction)
	Intersection<3> shaded_it(it);
	while( shaded_it.material()->is_type(Reflectance::DELTA) )
	{
		Real pdf;
		//Check termination condition...
		if ( m_termCriteria(shaded_it, pdf))
				return Spectrum(0.);	

		//Get output direction, plus the attenuation
		VectorN<D> new_omega_i;
		attenuation *= shaded_it.sample_direction(shaded_it.get_ray().get_direction(),
			new_omega_i, pdf );
		attenuation *= dot_abs(shaded_it.get_normal(), new_omega_i)/pdf;

		//Create outgoing ray
		Ray<D> new_ray(shaded_it.get_position(), new_omega_i, shaded_it.get_ray().get_level()+1,  shaded_it.get_ray().get_refraction_index());
		new_ray.shift();

		Intersection<D> new_it;
		//Get new intersection
		m_world->first_intersection(new_ray, new_it);

		if( !new_it.did_hit() )
				return Spectrum(0.);

		shaded_it = new_it;
	}
	// Compute direct light contribution
	reflected_radiance += m_direct_light_integrator(shaded_it);
	
	//Get VPL contribution
	LightSample<D> light_sample;	
	Real pdf_light; 

	//For each VPL 
	for( int light = 0; light < m_vpls.size(); ++light)
	{
		//Compute incoming light
		if( m_vpls[light].sample( shaded_it, light_sample, pdf_light ) )
		{
			Spectrum radiance = light_sample.irradiance*
				shaded_it.material()->f(light_sample.dir,-shaded_it.get_ray().get_direction())
				*dot_clamped(-light_sample.dir,shaded_it.get_normal())
				/pdf_light;

			reflected_radiance += radiance;
			//reflected_radiance += Vector3(0.,0.,1.);
		}
	}
	return reflected_radiance*attenuation;
}


template<class TermCriteria, int D>
void InstantRadiosity<TermCriteria,D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
	//Handle delta materials (e.g. perfect reflection/refraction)
	Spectrum attenuation(1.);
	Intersection<D> shaded_it(it);
	Real time = 0;

	while( shaded_it.material()->is_type(Reflectance::DELTA)  )
	{
		Real pdf;
		//Check termination condition...
		if ( m_termCriteria(shaded_it, pdf))
			return;	
		
		//Get output direction, plus the attenuation
		VectorN<D> new_omega_i;
		attenuation = shaded_it.sample_direction(shaded_it.get_ray().get_direction(),
			new_omega_i, pdf );
		attenuation *= dot_abs(shaded_it.get_normal(), new_omega_i)/pdf;

		//Create outgoing ray
		Ray<D> new_ray(shaded_it.get_position(), new_omega_i, shaded_it.get_ray().get_level()+1,  shaded_it.get_ray().get_refraction_index());
		new_ray.shift();

		Intersection<D> new_it;
		//Get new intersection
		m_world->first_intersection(new_ray, new_it);

		if( !new_it.did_hit() )
				return;

		shaded_it = new_it;
		time += m_world->time_of_flight(new_it.get_ray().get_parameter());
	}

	//Get VPL contribution
	LightSample<D> light_sample;	
	Real pdf_light; 
	//For each VPL 
	for( int light = 0; light < m_vpls.size(); ++light)
	{
		//Compute incoming light
		if( m_vpls[light].sample( shaded_it, light_sample, pdf_light ) )
		{
			Spectrum radiance = light_sample.irradiance*
				shaded_it.f(light_sample.dir,-shaded_it.get_ray().get_direction())
				*dot_clamped(-light_sample.dir,shaded_it.get_normal())
				/pdf_light;

			samples.push_back( RadianceSample(radiance, 
				m_world->time_of_flight(light_sample.dist)/it.get_ray().get_refraction_index(),  
				shaded_it.get_ray().get_level())); //Take care with the bounce of the VPL... fix it!
		}
	}
}

#endif //_INSTANT_RADIOSITY_H_
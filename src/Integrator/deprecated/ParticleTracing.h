#ifndef _PARTICLE_TRACING_H_
#define _PARTICLE_TRACING_H_

#include "Integrator/SurfaceIntegrator.h"
#include "RayTracing/World.h"
#include "RayTracing/Intersection.h"

#include "Particle.h"

#include <list>

/**	Virtual class to implement all techniques based on particle tracing,
	including Photon Mapping [Jensen 2001] and Instant Radiosity 
	[Keller 1997], and the algorithms that inherit from them.	*/
template<class TermCriteria, int D>
class ParticleTracing: public SurfaceIntegrator<D>
{
protected:
	TermCriteria m_termCriteria;
	int m_max_nb_photons_shot, m_nb_global_photons, m_nb_caustics_photons;
	
	Real m_max_time;
	bool m_use_russian_roulette;


	/**	Trace particles and store their intersections	*/
	void trace_photons(std::list<Particle<D> > &photons, std::list<Particle<D> > &caustic_photons) const;

public:
	ParticleTracing(World<D> *w,  unsigned int incoming_samples = 1,  
		int nb_global_photons = std::numeric_limits<int>::max(),
				  int nb_caustics_photons = std::numeric_limits<int>::max(),
				  int max_nb_shots = 500000)
		: SurfaceIntegrator<D>(w,incoming_samples), m_max_nb_photons_shot(max_nb_shots),
		m_nb_global_photons(nb_global_photons), m_nb_caustics_photons(nb_caustics_photons),
		m_max_time(1e20),m_use_russian_roulette(true)
	{}

	~ParticleTracing() {};

	/** Virtual functions to store or load the particles. Each method
		store the particles in a different way, depending if they
		are Photons, VPLs, Rays or even full hierarchies.	*/
	virtual void store(const std::string &namefile)const = 0;
	virtual void load(const std::string &namefile)const = 0;


	void set_max_time(Real max_time){m_max_time=max_time; m_use_russian_roulette=false;}
}; //ParticleTracing

template<class TermCriteria, int D>
void ParticleTracing<TermCriteria,D>::trace_photons(std::list<Particle<D> > &global_photons, std::list<Particle<D> > &caustic_photons) const
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
		if( ++nb_shots > m_max_nb_photons_shot )
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
						if( caustic_photons.size() == m_nb_caustics_photons )
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
}

#endif //_PARTICLE_TRACING_H_
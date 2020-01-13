/*
 *  BeamTracing.h
 *  transient
 *
 *  Created by Julio on 27/06/2013.
 *
 */

#ifndef _BEAM_TRACING_H_
#define _BEAM_TRACING_H_

#include "Integrator/VolumeIntegrator.h"
#include "RayTracing/World.h"
#include "RayTracing/Intersection.h"
#include "Photon.h"
#include "Beam.h"
#include "Utils/RandomNumbers.h"

#include <list>
#define isNaN(num) (num != num) //IEEE standard

template<class TermCriteria, int D, BeamBlurRegion BR>
class BeamTracing: public VolumeIntegrator<D>
{
protected:
	TermCriteria m_termCriteria;
	int m_max_nb_shots, m_nb_global_beams, m_nb_caustics_beams, m_nb_global_photons, m_nb_caustics_photons;
	Real global_beam_radius, caustic_beam_radius, beam_max_length;
	bool rm_single_scattering;
	int b_subdivide;
	Scattering::Level scat_level;
	Real TMAX;
    Film *film;
	/**	Trace beams	*/
public:
	BeamTracing(World<D> *w, Real _global_beam_radius=0.3, Real _caustic_beam_radius=0.03, unsigned int incoming_samples = 1,  
				int nb_global_photons = std::numeric_limits<int>::max(),
				int nb_caustics_photons = std::numeric_limits<int>::max(),
				int nb_global_beams = std::numeric_limits<int>::max(),
				int nb_caustics_beams = std::numeric_limits<int>::max(),
				int max_nb_shots = 500000, int _b_subdivide = 0, bool _rm_single_scattering = false, 
				Scattering::Level _scat_level = Scattering::ALL, Real _TMAX = std::numeric_limits<Real>::infinity(), Film *_film = NULL)
				: VolumeIntegrator<D>(w), m_max_nb_shots(max_nb_shots), 
	m_nb_global_photons(nb_global_photons), m_nb_caustics_photons(nb_caustics_photons), 
	m_nb_global_beams(nb_global_beams), m_nb_caustics_beams(nb_caustics_beams), 
	global_beam_radius(_global_beam_radius), caustic_beam_radius(_caustic_beam_radius), beam_max_length(300), 
	b_subdivide(_b_subdivide), rm_single_scattering(_rm_single_scattering), scat_level(_scat_level), TMAX(_TMAX), film(_film)
	{}
	
	~BeamTracing() {};
	
	bool is_level (Scattering::Level level) const { return (level & scat_level); }
	//virtual void preprocess();
	virtual void store(const std::vector<Beam<D,BR> > beams, const std::string &namefile)const;
	//virtual void load(const std::string &namefile)const = 0;
	
	void trace_beams(std::vector<Beam<D,BR> > &global_beams, 
					 std::vector<Beam<D,BR> > &caustic_beams,
					 std::list<Photon<D> > &global_photons, 
					 std::list<Photon<D> > &caustic_photons) const;
	void operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const{}// = 0;
	Spectrum operator()(const Ray<D> &r) const {return Spectrum(0.);};
	
	void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const {}
	Spectrum operator()(const Intersection<D> &it) const {return Spectrum(0.);};
	
}; //BeamTracing

template<class TermCriteria, int D, BeamBlurRegion BR>
void BeamTracing<TermCriteria,D,BR>::store(const std::vector<Beam<D,BR> > beams, const std::string &namefile) const
{
	FILE *f = fopen(namefile.c_str(), "w");
	typename std::vector<Beam<D,BR> >::const_iterator it;
	for (it = beams.begin(); it != beams.end(); it++)
	{
		//[ORIGIN] [DIRECTION] [PARAMETER RADIUS]
		VectorN<D> o = it->get_origin(), w = it->get_direction();
		Real r = it->get_radius(), t = it->get_parameter();
		if (D == 2)
			fprintf(f, "[%f %f] [%f %f] [%f %f] power: %f level: %d time: %f\n", o[0], o[1], w[0], w[1], 
					r, t, it->power_at(0.)[0], it->get_level(), it->get_time_of_flight());
		else
			fprintf(f, "[%f %f %f] [%f %f %f] [%f %f] power: %f level: %d time: %f\n", o[0], o[1], o[2], w[0], w[1], w[2],
					r, t, it->power_at(0.)[0], it->get_level(), it->get_time_of_flight());
		
	}
	fclose(f);
}

template<class TermCriteria, int D, BeamBlurRegion BR>
void BeamTracing<TermCriteria,D,BR>::trace_beams(std::vector<Beam<D,BR> > &global_beams,
												 std::vector<Beam<D,BR> > &caustic_beams,
												 std::list<Photon<D> > &global_photons, 
												 std::list<Photon<D> > &caustic_photons) const
{
	unsigned int nb_shots = 0;
	//unsigned int nb_shots_beams = 0, nb_shots_photons = 0;
	unsigned int nb_shots_caustics = 0, nb_shots_caustics_photons = 0;
	unsigned int nb_shots_global = 0, nb_shots_global_photons = 0;
	bool shooting = true;
	int i = 0;
	bool uniform_sampling = TMAX != std::numeric_limits<Real>::infinity();

	//Shoot beams!
	//Check if finished...
	while(shooting)
	{
		i++;
		//printf("Shooting photon number %d...\n", i++);
		//Check if max number of shots done...
		nb_shots++;
		//if(!nb_shots_caustics && !nb_shots_global) ++nb_shots_beams;
		//if(!nb_shots_caustics_photons && !nb_shots_global_photons) ++nb_shots_photons;

		if( nb_shots > m_max_nb_shots)
		{
			nb_shots_caustics = nb_shots_caustics?nb_shots_caustics:(nb_shots-1);
			nb_shots_global = nb_shots_global?nb_shots_global:(nb_shots-1);
			nb_shots_caustics_photons = nb_shots_caustics_photons?nb_shots_caustics_photons:(nb_shots-1);
			nb_shots_global_photons = nb_shots_global_photons?nb_shots_global_photons:(nb_shots-1);
			break;
		}
		
		// If not, start a new photon-beam path
		// Sample light source, and from it sample the photon origin
		LightSample<D> photon_origin;
		Real pdf_light1, pdf_light2;
		this->m_world->sample_light(pdf_light1)->sample(photon_origin,pdf_light2);
		
		// Compute irradiance photon's energy
		Spectrum energy(photon_origin.irradiance/(pdf_light1*pdf_light2));
		
		Real time = photon_origin.instant; //picoseconds(10^-12)
		
		Ray<D> photon_ray(photon_origin.pos, photon_origin.dir, DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
		photon_ray.shift();
		
		bool is_caustic_particle = false, hits_surface = false;
		//Iterate the path
		while(1)
		{
			// Throw ray and update current_it
			Intersection<D> it;
			
			Real epsilon1;
			Spectrum u_t, u_s, scat_prob;
			Real avg_scat_prob;
			Real sampled_mfp = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
			
			Real dist_TMAX;
			Real c = 1./this->m_world->time_of_flight(1);
			hits_surface = false;

			if (photon_ray.get_medium()) //SOME SCATTERING MEDIA -> sample mean free path
			{
				u_t = photon_ray.get_medium()->get_extinction(photon_ray.get_origin());
				//dist_TMAX = 1./u_t.avg();
				dist_TMAX = (TMAX - time)/this->m_world->time_of_flight(1);
				epsilon1=RNG::StdRandom::get_real();
				if (uniform_sampling) 
					sampled_mfp = epsilon1*dist_TMAX;
				else
					sampled_mfp = -logf(epsilon1)/u_t.avg();
			}
			
			//intersect with some surface
			this->m_world->first_intersection(photon_ray, it);
			
			if (photon_ray.get_medium() 
				&& ((is_level(Scattering::SINGLE) && !rm_single_scattering && photon_ray.get_level() == 0)
					|| (is_level(Scattering::MULTIPLE) && photon_ray.get_level() > 0)))
			{
				//Real max_time_length = min((TMAX - time)/this->m_world->time_of_flight(1.), min(beam_max_length, photon_ray.get_parameter()));
				Real max_time_length = min(beam_max_length, photon_ray.get_parameter());

				if (is_caustic_particle){
					if( !nb_shots_caustics )
					{
						//if( i == m_nb_caustics_beams )
						if( caustic_beams.size() == m_nb_caustics_beams )
						{	nb_shots_caustics = nb_shots; break;	}
						Photon<D> *p = new Photon<D>(photon_ray.get_origin() ,photon_ray.get_direction(), 
													 energy, time, it.get_ray().get_level());
						Beam<D,BR> b (p, caustic_beam_radius, max_time_length, true, photon_ray.get_level(), 
													   DEFAULT_REFRACTION_INDEX, photon_ray.get_medium());
						caustic_beams.push_back(b);
					}
				}
				else{
					if( !nb_shots_global ) 
					{
						//if( i == m_nb_global_beams )
						if( global_beams.size() == m_nb_global_beams )
						{	nb_shots_global = nb_shots; break;	}
						
						Photon<D> *p = new Photon<D>(photon_ray.get_origin(), photon_ray.get_direction(), 
													 energy, time, it.get_ray().get_level());
						
						Beam<D,BR> b (p, global_beam_radius, max_time_length, true, photon_ray.get_level(),
													   DEFAULT_REFRACTION_INDEX, photon_ray.get_medium());
						global_beams.push_back(b);
					}
				}
			}
			
			if (it.did_hit() && sampled_mfp > photon_ray.get_parameter())
			{ sampled_mfp = photon_ray.get_parameter(); hits_surface = true; }
			
			if (!hits_surface && photon_ray.get_medium()) //SAMPLE NEXT SCATTERING EVENT
			{
				u_s = photon_ray.get_medium()->get_scattering(photon_ray.get_origin());
				scat_prob = u_s/u_t;
				avg_scat_prob = scat_prob.avg();
				
				/* This is just a test!! */
				Real epsilon2 = RNG::StdRandom::get_real();

				if (uniform_sampling && dist_TMAX  > 1e-3) 
				//if (uniform_sampling && time + this->m_world->time_of_flight(sampled_mfp)  < TMAX) 
				{
					Real pdf_dist_TMAX = 1./(dist_TMAX);
					time += this->m_world->time_of_flight(sampled_mfp)*photon_ray.get_refraction_index();
					VectorN<D> new_omega_i, new_position = photon_ray.get_origin()+photon_ray.get_direction()*sampled_mfp;
					Real pdf;
					energy *= exp(-u_t*sampled_mfp)*scat_prob*photon_ray.get_medium()->sample_direction(new_position, photon_ray.get_direction(), new_omega_i, pdf);
					energy /= pdf*pdf_dist_TMAX; 

					photon_ray = Ray<D>(new_position, new_omega_i, photon_ray.get_level()+1, photon_ray.get_refraction_index(), photon_ray.get_medium());
					is_caustic_particle = false;
				}
				else if (!uniform_sampling)
				{	
					/*if(epsilon2 <= avg_scat_prob) //RUSSIAN ROULETTE: SCATTERED => SAMPLE NEW DIRECTION
					{
						time += this->m_world->time_of_flight(sampled_mfp)*photon_ray.get_refraction_index();
						VectorN<D> new_omega_i, new_position = photon_ray.get_origin()+photon_ray.get_direction()*sampled_mfp;
						Real pdf;
						energy *= scat_prob*photon_ray.get_medium()->sample_direction(new_position, photon_ray.get_direction(), new_omega_i, pdf);
						energy /= (pdf*avg_scat_prob);
						photon_ray = Ray<D>(new_position, new_omega_i, photon_ray.get_level()+1, photon_ray.get_refraction_index(), photon_ray.get_medium());
						is_caustic_particle = false;
					}*/
					//if (photon_ray.get_level() < 4)
					if (film != NULL) 
					{
						if (time + this->m_world->time_of_flight(sampled_mfp) < film->get_time_length()) //STOP WHEN OUT OF FILM
						{
							time += this->m_world->time_of_flight(sampled_mfp)*photon_ray.get_refraction_index();
							VectorN<D> new_omega_i, new_position = photon_ray.get_origin()+photon_ray.get_direction()*sampled_mfp;
							Real pdf;
							energy *= scat_prob*photon_ray.get_medium()->sample_direction(new_position, photon_ray.get_direction(), new_omega_i, pdf);
							energy /= pdf;
							photon_ray = Ray<D>(new_position, new_omega_i, photon_ray.get_level()+1, photon_ray.get_refraction_index(), photon_ray.get_medium());
							is_caustic_particle = false;
						}
						else break;
					} 
					else 
					{
						if(epsilon2 <= avg_scat_prob) //RUSSIAN ROULETTE: SCATTERED => SAMPLE NEW DIRECTION
						{
							time += this->m_world->time_of_flight(sampled_mfp)*photon_ray.get_refraction_index();
							VectorN<D> new_omega_i, new_position = photon_ray.get_origin()+photon_ray.get_direction()*sampled_mfp;
							Real pdf;
							energy *= scat_prob*photon_ray.get_medium()->sample_direction(new_position, photon_ray.get_direction(), new_omega_i, pdf);
							energy /= (pdf*avg_scat_prob);
							photon_ray = Ray<D>(new_position, new_omega_i, photon_ray.get_level()+1, photon_ray.get_refraction_index(), photon_ray.get_medium());
							is_caustic_particle = false;
						}
						else break;
					}
				}
				else //ABSORBED => STOP 
					break;
			}
			else if (hits_surface) //STORE PHOTON AND SAMPLE DIRECTION FROM SURFACE
			{ 
				//Update time
				time += this->m_world->time_of_flight(photon_ray.get_parameter())*photon_ray.get_refraction_index();
				
				//Check if has hit a delta material...
				if( it.material()->is_type(Reflectance::DELTA) )
				{
					// If delta material, then is caustic...
					// Don't store the photon!
					is_caustic_particle = true;
				}
				else if ((is_level(Scattering::SINGLE) && !rm_single_scattering && photon_ray.get_level() == 0)
						|| (is_level(Scattering::MULTIPLE) && photon_ray.get_level() > 0))
						
				{
					//If non-delta material, store the photon!
					if( is_caustic_particle )	//If caustic particle, store in caustics
					{				
						if( !nb_shots_caustics_photons )
						{
							if( caustic_photons.size() == m_nb_caustics_photons )
							{	nb_shots_caustics_photons = nb_shots; break;	}
							
							caustic_photons.push_back( Particle<D>(it, energy, time, it.get_ray().get_level()) );
						}
					}
					else						//If non-caustic particle, store in global
					{
						if( !nb_shots_global_photons ) 
						{
							if( global_photons.size() == m_nb_global_photons )
							{	nb_shots_global_photons = nb_shots; break;	}
							
							global_photons.push_back( Particle<D>(it, energy, time, it.get_ray().get_level()) );
						}
					}
					is_caustic_particle = false;
				}
				
				//Test termination condition
				Real pdf;
				/*if ( m_termCriteria(it, pdf))
				 break;
				 energy /= pdf;
				 */
				
				//uncomment this for Lambertian
				Spectrum surf_albedo = Spectrum(1.)-it.material()->get_absorption(it.get_uv());
				
				Real avg_surf_albedo = surf_albedo.avg();
				
				Real epsilon2 = RNG::StdRandom::get_real();
				/*Real epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				while (epsilon2 < 0.)
					epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				*/
				if (epsilon2 > avg_surf_albedo) 
					break;
				
				// Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				if(it.material()->is_type(Reflectance::DELTA))
				{
					energy*=surf_albedo;
					//it.material()->sample_direction(it.get_ray().get_direction(), new_omega_i, 
					//								it.get_normal(), it.get_uv(),pdf );	
					it.material()->sample_outgoing_ray(it, photon_ray, pdf);
					if (dot(it.get_normal(), photon_ray.get_direction()) > 0.)
						photon_ray.set_medium(this->m_world->get_medium());
					else
						photon_ray.set_medium(NULL);

					energy /= avg_surf_albedo;
				}	
				else{
					energy *= it.material()->sample_outgoing_ray(it, photon_ray, pdf);
					energy *= dot_abs(it.get_normal(), photon_ray.get_direction());
					//energy *= it.material()->sample_direction(it.get_ray().get_direction(), new_omega_i, 
					//								it.get_normal(), it.get_uv(),pdf );
					//energy *= dot_abs(it.get_normal(), new_omega_i);
					
					energy /= avg_surf_albedo*pdf;
				}
				photon_ray.shift();
			}
			else break; //don't hit surface and no medium => the photon is lost
		}
		
		if( nb_shots_global && nb_shots_caustics && nb_shots_global_photons &&  nb_shots_caustics_photons)
			shooting = false;
	}
	
	// Finalize photon shooting
	Real scale_beams_global, scale_beams_caustics, scale_photons_global, scale_photons_caustics;
	
	printf("nb_shots_caustics: %d\n", nb_shots_caustics);
	printf("nb_shots_global: %d\n", nb_shots_global);
	printf("nb_shots_caustics_photons: %d\n", nb_shots_caustics_photons);
	printf("nb_shots_global_photons: %d\n", nb_shots_global_photons);
	
	typename std::list<Photon<D> >::iterator it;
	typename std::vector<Beam<D,BR> >::iterator b_it;
	
	scale_beams_global = 1./static_cast<Real>(nb_shots_global);
	scale_beams_caustics = 1./static_cast<Real>(nb_shots_caustics);
	scale_photons_global = 1./static_cast<Real>(nb_shots_global_photons);
	scale_photons_caustics = 1./static_cast<Real>(nb_shots_caustics_photons);
	//Scale caustics' photons
	for( it = caustic_photons.begin(); it != caustic_photons.end(); it++ )
		it->m_power *= scale_photons_caustics;
	
	for( b_it = caustic_beams.begin(); b_it != caustic_beams.end(); b_it++ )
		b_it->scale_power(scale_beams_caustics);
	
	//Scale global' photons
	for( it = global_photons.begin(); it != global_photons.end(); it++ )
		it->m_power *= scale_photons_global;
	
	for( b_it = global_beams.begin(); b_it != global_beams.end(); b_it++ )
		b_it->scale_power(scale_beams_global);
	
	//SUBDIVIDING BEAMS
	if (b_subdivide > 0)
	{
		std::vector<Beam<D,BR> > subbeams;
		typename std::vector<Beam<D,BR> >::iterator sub_it;
		typename std::vector<Beam<D,BR> >::iterator beam_it;
		
		std::vector<Beam<D,BR> > old_caustic_beams = caustic_beams, old_global_beams = global_beams;
		
		caustic_beams.clear();
		global_beams.clear();
		
		for(beam_it = old_global_beams.begin(); beam_it != old_global_beams.end(); beam_it++)
		{	
			beam_it->subdivide(subbeams, b_subdivide); 
			for (sub_it = subbeams.begin(); sub_it != subbeams.end(); sub_it++)
				global_beams.push_back(*sub_it);
			subbeams.clear();
		
		}
		
		for(beam_it = old_caustic_beams.begin(); beam_it != old_caustic_beams.end(); beam_it++)
		{	
			beam_it->subdivide(subbeams, b_subdivide); 
			for (sub_it = subbeams.begin(); sub_it != subbeams.end(); sub_it++)
				caustic_beams.push_back(*sub_it);
			subbeams.clear();
		}
	}
}

#endif //_BEAM_TRACING_H_

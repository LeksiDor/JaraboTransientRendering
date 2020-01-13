#ifndef _VOLUME_PATH_TRACING_INTEGRATOR_H_
#define _VOLUME_PATH_TRACING_INTEGRATOR_H_

#include "Integrator/VolumeIntegrator.h"
#include "RayTracing/World.h"
#include "deprecated/TimeSampler.h"

/* For the moment just a homogeneous medium path tracer */
//Always keep one of this flags to 1
#define MIS	0		
#define STRAT 0
#define TIME 1

template<class TermCriteria, int D>
class VolumePathTracing: public VolumeIntegrator<D>
{
    Film *film;
	TermCriteria m_termCriteria;
	unsigned int m_incoming_samples;
	Real m_inv_incoming_samples;
	FILE *f_log;
public:
	VolumePathTracing(World<D> *w, const TermCriteria &termination,  unsigned int incoming_samples = 1, FILE *_f_log = NULL, Film *_film = NULL)
	:VolumeIntegrator<D>(w), m_termCriteria(termination), m_incoming_samples(incoming_samples),
	m_inv_incoming_samples(1./static_cast<Real> (m_incoming_samples)), f_log(_f_log), film(_film)
	{}
	~VolumePathTracing()
	{}
	void preprocess() 
	{
		if (f_log)
		{
			fprintf(f_log, "Volumetric Path Tracing\n");
			fprintf(f_log, "=======================\n");
			fprintf(f_log, "# Samples: %d\n", m_incoming_samples);
			fprintf(f_log, "\n");
		}
	}
	void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const;

	void operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Spectrum operator()(const Ray<D> &r) const;// {return Spectrum(0.);};
	
	void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Spectrum operator()(const Intersection<D> &it) const;
};

template<class TermCriteria, int D>
Spectrum VolumePathTracing<TermCriteria,D>::operator()(const Ray<D> &r) const
{
	Spectrum reflected_radiance(0.);
	
	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		//NEW RAY
		Ray<D> new_ray(r.get_origin(), r.get_direction(), DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
		new_ray.set_parameter(r.get_parameter());
		
		//SAMPLE MEAN FREE PATH
		Real epsilon1, epsilon2;
		Spectrum u_t, u_s, scat_albedo;
		Real avg_scat_albedo;
		Real sampled_mfp = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
		
		Spectrum attenuation(1.);
		Intersection<D> current_it;

		int iteration = 0;
		while(1)
		{
			if(new_ray.get_medium())
			{	
				u_t = new_ray.get_medium()->get_extinction(new_ray.get_origin());
				u_s = new_ray.get_medium()->get_scattering(new_ray.get_origin());
				epsilon1 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				while (epsilon1 < 1e-5)
					epsilon1 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				sampled_mfp = -logf(epsilon1)/u_t.avg();
				scat_albedo = u_s/u_t;
				avg_scat_albedo = scat_albedo.avg();
			}
			else
			{
				u_t = Spectrum(0.);
				u_s = Spectrum(0.);
				sampled_mfp = std::numeric_limits<Real>::infinity();
			}
			Real pdf_light1, pdf_light2;
			if (!current_it.did_hit() || sampled_mfp < new_ray.get_parameter()) //SCATTERING EVENT
			{
				VectorN<D> p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
				
				//Get light source 
				
				const LightSource<D> *l = this->m_world->sample_light(pdf_light1);			
				LightSample<D> ls;
				Real pdf_light2;
				if(l->sample(p, ls, pdf_light2))
					reflected_radiance +=  (u_s)*attenuation*exp(-u_t*ls.dist)*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())/(pdf_light1*pdf_light2);
				
				Real pdf;
				/*if ( m_termCriteria(current_it, pdf))
				 break;
				 attenuation /= pdf;
				 */
				
				epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				while (epsilon1 < 1e-5)
					epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				
				if (epsilon2 > avg_scat_albedo) 
					break;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				Medium<D> *m = new_ray.get_medium();
				
				attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
				attenuation /= pdf*avg_scat_albedo;
				
				new_ray = Ray<D>(p, new_omega_i, DEFAULT_REFRACTION_INDEX, m);
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				
				current_it = new_it;
			}
			else //SURFACE HIT
			{
				const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);
				
				//Get direct contribution
				//if( !current_it.delta_material() )
				LightSample<D> light_sample;
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
					Spectrum radiance = light_sample.irradiance*current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())
										*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation*exp(-u_t*light_sample.dist)/(pdf_light1*pdf_light2);
					
					reflected_radiance += radiance;
				}
				
				//Test termination condition
				Real pdf;
				if ( m_termCriteria(current_it, pdf))
					break;
				attenuation /= pdf;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				attenuation *= current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
				attenuation *= dot_abs(current_it.get_normal(), new_omega_i)/pdf; 
				
				// Throw ray and update current_it
				if (dot(current_it.get_normal(), new_omega_i) >= 0.) //new ray outside the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
				else //new ray into the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
				
				new_ray.shift();
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				if (!new_it.did_hit() && !new_ray.get_medium()) 
				{	
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
					new_ray.shift();
				}
				current_it = new_it;
			}
		}
	}
	return reflected_radiance*m_inv_incoming_samples;
}


template<class TermCriteria, int D>
Spectrum VolumePathTracing<TermCriteria,D>::operator()(const Intersection<D> &it) const
{
	Spectrum reflected_radiance(0.);
	
	//ANALYTICAL SOLUTION
#if 0
{	Real pdf_light1, pdf_light2;
	const LightSource<D> *l = this->m_world->sample_light(pdf_light1);
	Spectrum u_t, u_s;
	LightSample<D> ls;
	
	VectorN<D> p = it.get_ray().get_position(); //point of scattering event
	
	l->sample(p, ls, pdf_light2);
	u_t = this->m_world->get_medium()->get_extinction(VectorN<D>(0.));
	u_s = this->m_world->get_medium()->get_scattering(VectorN<D>(0.));
	
	reflected_radiance += ls.irradiance*it.material()->f(ls.dir,-it.get_ray().get_direction(),it.get_normal(),it.get_uv())*dot_clamped(-ls.dir,it.get_normal())*exp(-u_t*it.get_ray().get_parameter())/(pdf_light1*pdf_light2);
	reflected_radiance += l->get_intensities()*this->m_world->get_medium()->f(p,ls.dir,-it.get_ray().get_direction())*u_s*(1.-exp(-u_t*it.get_ray().get_parameter()))/u_t;
	
	return reflected_radiance;
}
#endif
	//PATH TRACING
	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		//NEW RAY
		Intersection<D> current_it(it);
		Ray<D> new_ray(current_it.get_ray().get_origin(), current_it.get_ray().get_direction(), DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
		new_ray.set_parameter(current_it.get_ray().get_parameter());
		
		//SAMPLE MEAN FREE PATH
		Real epsilon1, epsilon2;
		Spectrum u_t, u_s, scat_albedo;
		Real avg_scat_albedo, pdf_light1, pdf_light2;
		Real sampled_mfp = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
		
		Spectrum attenuation(1.);
		while(1)
		{
			if(new_ray.get_medium())
			{	
				u_t = new_ray.get_medium()->get_extinction(new_ray.get_origin());
				u_s = new_ray.get_medium()->get_scattering(new_ray.get_origin());
				epsilon1 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				while (epsilon1 < 1e-5)
					epsilon1 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				sampled_mfp = -logf(epsilon1)/u_t.avg();
				scat_albedo = u_s/u_t;
				avg_scat_albedo = scat_albedo.avg();
			}
			else
			{
				u_t = Spectrum(0.);
				u_s = Spectrum(0.);
				sampled_mfp = std::numeric_limits<Real>::infinity();
			}
			
			
			if (!current_it.did_hit() || sampled_mfp < new_ray.get_parameter()) //SCATTERING EVENT
			{
				VectorN<D> p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
				
				//Get light source 
				const LightSource<D> *l = this->m_world->sample_light(pdf_light1);			
				LightSample<D> ls;
				Real pdf_light2;
				if(l->sample(p, ls, pdf_light2))
					reflected_radiance += (u_s)*attenuation*exp(-u_t*ls.dist)*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())/(pdf_light1*pdf_light2);
				
				
				Real pdf;
				/*if ( m_termCriteria(current_it, pdf))
				 break;
				 attenuation /= pdf;
				 */
				
				epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				while (epsilon1 < 1e-5)
					epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
				
				if (epsilon2 > avg_scat_albedo) 
					break;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				Medium<D> *m = new_ray.get_medium();
				
				attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
				attenuation /= pdf*avg_scat_albedo;
				
				new_ray = Ray<D>(p, new_omega_i, DEFAULT_REFRACTION_INDEX, m);
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				
				current_it = new_it;
			}
			else //SURFACE HIT
			{
				const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);

				//Get direct contribution
				//if( !current_it.delta_material() )
				LightSample<D> light_sample;
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
					Spectrum radiance = light_sample.irradiance*current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())
					*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation*exp(-u_t*light_sample.dist)/(pdf_light1*pdf_light2);
					
					reflected_radiance += radiance;
				}
				
				//Test termination condition
				Real pdf;
				if ( m_termCriteria(current_it, pdf))
					break;
				attenuation /= pdf;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				
				if(it.material()->is_type(Reflectance::DELTA_TRANSMISSION))
				{
					current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
					attenuation *= dot_abs(current_it.get_normal(), new_omega_i);
				}
				else{
					attenuation *= current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
					attenuation *= dot_abs(current_it.get_normal(), new_omega_i)/pdf; 
				}
				
				// Throw ray and update current_it
				if (dot(current_it.get_normal(), new_omega_i) >= 0.) //new ray outside the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
				else //new ray into the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
				
				new_ray.shift();
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				
				if (!new_it.did_hit() && !new_ray.get_medium()) //singularities when new out direction and normal are perpendicular
				{	
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
					new_ray.shift();
				}
				
				current_it = new_it;
			}
		}
	}
	return reflected_radiance*m_inv_incoming_samples;
}

template<class TermCriteria, int D>
void VolumePathTracing<TermCriteria,D>::operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
    Real offset_time =  film?film->get_time_offset():0.,
                        film_time_length = film?film->get_time_length():0., 
                        film_exp = film->get_exposure_time();
	Real max_time = (nb_time_samples-1) * delta_time;
	Real inv_delta_time = 1./delta_time, ray_dist(-1.);
	bool multiple_time_samples = (nb_time_samples > 1);
	
	std::vector<Spectrum> time_samples;
	if( multiple_time_samples )
		time_samples.resize(nb_time_samples);

    if (!film) printf("WARNING: Film should be != NULL\n");

	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		//NEW RAY
		Ray<D> new_ray(r.get_origin(), r.get_direction(), DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
		new_ray.set_parameter(r.get_parameter());
		Real time_mfp;
		//SAMPLE MEAN FREE PATH
		Real epsilon1, epsilon2;
		Spectrum u_t, u_s, scat_albedo;
		Real avg_scat_albedo;
		Real sampled_mfp = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
		Spectrum attenuation(1.);
		Intersection<D> current_it;
		
		Real time = 0., pf_time = 0.;
		ray_dist = -1.;

		//TIME SAMPLING VARIABLES
		Real t_path = 0.0f, T_MAX = film_time_length + offset_time;
		while(1)
		{
			if(new_ray.get_medium())
			{	
				u_t = new_ray.get_medium()->get_extinction(new_ray.get_origin());
				u_s = new_ray.get_medium()->get_scattering(new_ray.get_origin());
				
				epsilon1 = RNG::StdRandom::get_real();
				while (epsilon1 < 1e-5) 
					epsilon1 = RNG::StdRandom::get_real();
				
				sampled_mfp = -logf(epsilon1)/u_t.avg();
				scat_albedo = u_s/u_t;
				avg_scat_albedo = scat_albedo.avg();
				time_mfp = this->m_world->time_of_flight(1.)/u_t.avg();

			}
			else
			{
				u_t = Spectrum(0.);
				u_s = Spectrum(0.);
				sampled_mfp = std::numeric_limits<Real>::infinity();
			}
			
			Real pdf_light1, pdf_light2;
			//only in the beginning of the path starting in the camera ray
			if (!current_it.did_hit() || sampled_mfp < new_ray.get_parameter()) //SCATTERING EVENT
			{
				if( nb_time_samples == 0) //RADIANCE SAMPLING
				{
					if (ray_dist < 0.) ray_dist = std::min(sampled_mfp, r.get_parameter());
					time += this->m_world->time_of_flight(sampled_mfp);
                    //transient PF 
                    time += pf_time;

                    if (time > film_time_length + offset_time)
						break;
                    
					VectorN<D> p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
				
					//Get light source 
				
					const LightSource<D> *l = this->m_world->sample_light(pdf_light1);			
					LightSample<D> ls;
					Real pdf_light2, pf_pdf;
				
					if(l->sample(p, ls, pdf_light2))
					{	
                        /* //Old steady-state PF sampling.
                        Spectrum radiance = scat_albedo*attenuation*exp(-u_t*ls.dist)*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;
                        radiance /= (pf_pdf*pdf_light1*pdf_light2);
                         Real sample_time = time + this->m_world->time_of_flight(ls.dist)/new_ray.get_refraction_index() + ls.instant;
                         */
                        
                        //Transient PF sampling. To use the PF in steady state, just set the mean_time to 0 in the PF constructor.
                        Spectrum radiance = scat_albedo*attenuation*exp(-u_t*ls.dist)*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction(), pf_time, pf_pdf)*m_inv_incoming_samples;

                        radiance /= (pf_pdf*pdf_light1*pdf_light2);
                        //printf("Radiance = %f\n", radiance.avg());
                        
                        //sample_time takes into account the length of the path followed by several bounces and the time the light takes to travel until the scattering event
						Real sample_time = time + pf_time + this->m_world->time_of_flight(ls.dist)/new_ray.get_refraction_index() + ls.instant;
					
						if ( multiple_time_samples )
						{
							if( sample_time < max_time )
							{	
								Real time_sample = sample_time * inv_delta_time;
								int min_sample = static_cast<int>(time_sample);
							
								time_samples[min_sample] += radiance*(1.-time_sample+min_sample);
								time_samples[min_sample+1] += radiance*(time_sample-min_sample);
							}
						}
						else
							samples.push_back( RadianceSample(radiance, 
															  sample_time,  
															  new_ray.get_level(),
															  this->m_world->time_of_flight(ray_dist)));
					}
				
					Real pdf;
					/*if ( m_termCriteria(current_it, pdf))
					 break;
					 attenuation /= pdf;
					 */
					epsilon2 = RNG::StdRandom::get_real();
				
					/*epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
					while (epsilon2 < 1e-5)
						epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
					*/
					if (epsilon2 > avg_scat_albedo) 
						break;
				
					//Random walk's next step
					// Get sampled direction plus pdf, and update attenuation
					VectorN<D> new_omega_i;
					Medium<D> *m = new_ray.get_medium();
				
					//attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
                    
                    //transient PF 
                    attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pf_time, pdf);
					attenuation /= pdf*avg_scat_albedo;
				
					new_ray = Ray<D>(p, new_omega_i, new_ray.get_level()+1, DEFAULT_REFRACTION_INDEX, m);
				
					Intersection<D> new_it;
					this->m_world->first_intersection(new_ray, new_it);
				
					current_it = new_it;
				}
				else //TIME SAMPLING
				{
					TimeSampler<D> time_sampler (this->m_world, new_ray.get_refraction_index());
					Real ta, tb;
					const LightSource<D> *light = this->m_world->sample_light(pdf_light1);
					
					VectorN<D> l, omega;
					LightSample<D> ls;
					Real pdf_light2;
					light->sample(new_ray.get_origin(), ls, pdf_light2);
					
					l = ls.pos - new_ray.get_origin();
					omega = new_ray.get_direction();
                    ta = max(offset_time - t_path, this->m_world->time_of_flight(l.length())/new_ray.get_refraction_index());
                    //ta = this->m_world->time_of_flight(l.length())/new_ray.get_refraction_index();
					//tb = T_MAX;
					
					tb = (new_ray.did_hit()?min(this->m_world->time_of_flight(new_ray.get_parameter()), T_MAX):T_MAX) - t_path;
					Real sampled_dist, pdf_time, pf_pdf;
					VectorN<D> p;
#if MIS
					if (tb < ta || time > T_MAX) break;
					//MIS
					Real prob_mfp, prob_ts, sum_prob_pdf, w;
					prob_mfp = 1- pow(time/T_MAX,1);
					prob_ts = 1. - prob_mfp;
					if (RNG::StdRandom::get_real() < prob_mfp)
					{
						if (ray_dist < 0.) 
							ray_dist = std::min(sampled_mfp, r.get_parameter());
						pdf_time = time_sampler.pdf(sampled_mfp,ta,tb,l,omega);
						sum_prob_pdf = prob_mfp*exp(-u_t.avg()*sampled_mfp) + prob_ts*pdf_time;
						w = prob_mfp/sum_prob_pdf;
						p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
						if( light->sample(p, ls, pdf_light2) )
						{
							/* //Old steady-state PF sampling.
                             Spectrum radiance = w*attenuation*exp(-u_t*(ls.dist+sampled_mfp))*scat_albedo*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
							radiance /= pdf_light1*pdf_light2*prob_mfp;
							Real sample_time = time + this->m_world->time_of_flight(sampled_mfp+ls.dist);// + ls.instant;
                             */
                            
                            //Transient PF
                            Spectrum radiance = w*attenuation*exp(-u_t*(ls.dist+sampled_mfp))*scat_albedo*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction(), pf_time, pf_pdf)*m_inv_incoming_samples;	
							radiance /= pdf_light1*pdf_light2*prob_mfp*pf_pdf;
							Real sample_time = time + pf_time + this->m_world->time_of_flight(sampled_mfp+ls.dist);// + ls.instant;

							//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
							samples.push_back( RadianceSample(radiance, 
																sample_time,  
																new_ray.get_level(),
																this->m_world->time_of_flight(ray_dist)));
							

						}
					}
					else
					{
						time_sampler(l,omega,ta,tb,sampled_dist,pdf_time);
						sum_prob_pdf = prob_mfp*exp(-u_t.avg()*sampled_dist) + prob_ts*pdf_time;
						w = prob_ts/sum_prob_pdf;
						//printf("w = %f\n", w);
						p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
						if( light->sample(p, ls, pdf_light2) )
						{
							 //Old steady-state PF sampling.
                             /*Spectrum radiance = w*attenuation*exp(-u_t*(ls.dist+sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;
							radiance /= pdf_light1*pdf_light2*prob_ts;
							Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
                            */
							Spectrum radiance = w*attenuation*exp(-u_t*(ls.dist+sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction(), pf_time, pf_pdf)*m_inv_incoming_samples;
							radiance /= pdf_light1*pdf_light2*prob_ts*pf_pdf;
							Real sample_time = time + pf_time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
                            
                            //Transient PF
                            samples.push_back( RadianceSample(radiance, 
																sample_time,  
																new_ray.get_level(),
																this->m_world->time_of_flight(ray_dist < 0?sampled_dist:ray_dist)));
						}
						
						if (ray_dist < 0.) 
							ray_dist = std::min(sampled_mfp, r.get_parameter());
					}
#elif STRAT
					if (ta > tb) break;

					Real t_exp = this->m_world->time_of_flight(1/u_t.avg());
					//Real t_exp = delta_time/static_cast<Real>(nb_time_samples);
					Real ta_i, tb_i;
					ta_i = ta;
					tb_i = ta + t_exp;
					//printf("\n");
					while (tb_i < tb + t_exp){
						//printf("tb = %.3f, tb_i = %.3f\n",tb,tb_i);
						time_sampler(l,omega,ta_i, tb_i, sampled_dist, pdf_time);
						
						p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
						if( light->sample(p, ls, pdf_light2))
						{
							//printf("sampled_dist = %f, ta = %f, tb = %f\n", sampled_dist + ls.dist,ta,tb);
							Spectrum radiance = attenuation*exp(-u_t*(ls.dist + sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
							radiance /= pdf_light1*pdf_light2*pdf_time;
							//if (radiance.avg() < 1e-5) printf("attenuation = %f, u_s = %f, pdfs: %f %f %f\n", attenuation.avg(), u_s.avg(), pdf_light1, pdf_light2, pdf_time);

							Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
							//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

							//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
							samples.push_back( RadianceSample(radiance, 
																sample_time,  
																new_ray.get_level(),
																this->m_world->time_of_flight(ray_dist)));
						}
						ta_i += t_exp;
						tb_i += t_exp;
					}
#elif TIME
					if (ta >= tb) break;
					//sample per frame

					time_sampler(l,omega,ta,tb,sampled_dist,pdf_time);
					p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
					if( light->sample(p, ls, pdf_light2))
					{
						 //Old steady-state PF sampling
                        /* Spectrum radiance = attenuation*exp(-u_t*(ls.dist + sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
                        radiance /= pdf_light1*pdf_light2*pdf_time;
						Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
                         */
                        
                        // Transient PF
                        Spectrum radiance = attenuation*exp(-u_t*(ls.dist + sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction(), pf_time, pf_pdf)*m_inv_incoming_samples;	
                        radiance /= pdf_light1*pdf_light2*pdf_time*pf_pdf;
						Real sample_time = time + pf_time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
                        
						samples.push_back( RadianceSample(radiance, 
															sample_time,  
															new_ray.get_level(),
															this->m_world->time_of_flight(ray_dist)));
					}
#else				
                    //LINE TO LINE
					if (ta > tb) break;
					//printf("%f\n",ta);
					// Line-to-line sampling, j subscript is for light path, i subscript is for camera path
					Real pdf_tj, rj, rj_max, ri, ri_max, ra, rb;
					VectorN<D> wj, xj;
					// 1. Sample light direction
					Real pdf_wj;
					LightSample<D> ls_j;
					light->sample(ls_j, pdf_wj);
					wj = ls_j.dir;
					xj = ls_j.pos;
					rj_max =(T_MAX - t_path)/this->m_world->time_of_flight(1.);

					//2. compute interval [ra,rb] of first sampled path (we choose camera path)
					ri_max = (T_MAX -t_path)/this->m_world->time_of_flight(1.);
					ra = time_sampler.rx(ta,rj_max, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
					rb = time_sampler.rx(tb,0.0f, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
					
					//Real taux = time_sampler.t(ra, rj_max, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
					//if (fabs(taux - ta) > 100) break;// printf("rj = %f ra = %f, taux, ta = %f %f\n",rj, ra, taux, ta);
					//taux = time_sampler.t(rb, 0, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
					
					ra = max(ra,0.0f);
					rb = min(rb, ri_max);
					//3. Uniformly sample a position inside [ra,rb]
					ri = RNG::StdRandom::get_real()*(rb-ra);
					//printf("%f\n",ri);
					Real pdf_rba = 1./(rb-ra);
					VectorN<D> x_ri = new_ray.get_origin() + new_ray.get_direction()*ri;
					
					time_sampler(x_ri - xj, wj, ta, tb, rj, pdf_tj);
					VectorN<D> x_rj = xj + wj*rj;
					Real rij = (x_ri - x_rj).length();
					Real att_dist = 1./((D==3)?(rij*rij):(rij));
					//Assuming visibility, this should be changed
					Spectrum radiance = attenuation*att_dist*exp(-u_t*(ri+rj+rij))
										*u_s*u_s
										*ls_j.irradiance
										*this->m_world->get_medium()->f(x_rj, ls_j.dir, x_ri - x_rj)
										*this->m_world->get_medium()->f(x_ri, x_rj - x_ri, -new_ray.get_direction())
										*m_inv_incoming_samples;	
					radiance /= pdf_light1*pdf_light2*pdf_tj*pdf_rba*pdf_wj;
					//if (radiance.avg() < 1e-5) printf("attenuation = %f, u_s = %f, pdfs: %f %f %f\n", attenuation.avg(), u_s.avg(), pdf_light1, pdf_light2, pdf_time);

					Real sample_time = time + this->m_world->time_of_flight(ri+rj+rij);// + ls.instant;
					//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

					//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
					samples.push_back( RadianceSample(radiance, 
														sample_time,  
														new_ray.get_level(),
														this->m_world->time_of_flight(ray_dist)));
#endif

					p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp;
                    time +=  this->m_world->time_of_flight(sampled_mfp)/new_ray.get_refraction_index();
                    time += pf_time;
                    
                    t_path = time;
					Real pdf;
					//Random walk's next step
					// Get sampled direction plus pdf, and update attenuation
					VectorN<D> new_omega_i;
					Medium<D> *m = new_ray.get_medium();
                    
                    //Steady-state PF
                    attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
					
                    //Trasient PF
                    //attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pf_time, pdf);
                    attenuation /= pdf;

					new_ray = Ray<D>(p, new_omega_i, new_ray.get_level()+1, DEFAULT_REFRACTION_INDEX, m);
				
					Intersection<D> new_it;
					this->m_world->first_intersection(new_ray, new_it);
				
					current_it = new_it;
				}
				
			}
			else //SURFACE HIT
			{
				time += this->m_world->time_of_flight(current_it.get_ray().get_parameter());
				const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);
				
                if (time > film_time_length + offset_time)
                    break;
				//Get direct contribution
				//if( !current_it.delta_material() )
				LightSample<D> light_sample;
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
					Spectrum radiance = light_sample.irradiance*current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())
					*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation*exp(-u_t*light_sample.dist)/(pdf_light1*pdf_light2)*m_inv_incoming_samples;
					
					Real sample_time = time + this->m_world->time_of_flight(light_sample.dist)/current_it.get_ray().get_refraction_index() + light_sample.instant;
					
					/*if ( multiple_time_samples )
					{
						if( sample_time < max_time )
						{	
							Real time_sample = sample_time * inv_delta_time;
							int min_sample = static_cast<int>(time_sample);
							
							time_samples[min_sample] += radiance*(1.-time_sample+min_sample);
							time_samples[min_sample+1] += radiance*(time_sample-min_sample);
						}
					}
					else*/
						samples.push_back( RadianceSample(radiance, 
														  sample_time,  
														  current_it.get_ray().get_level(),
														  this->m_world->time_of_flight(ray_dist)));
				}
				
				//Test termination condition
				Real pdf;
				if ( m_termCriteria(current_it, pdf))
					break;
				attenuation /= pdf;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				attenuation *= current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
				attenuation *= dot_abs(current_it.get_normal(), new_omega_i)/pdf; 
				
				// Throw ray and update current_it
				if (dot(current_it.get_normal(), new_omega_i) >= 0.) //new ray outside the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
				else //new ray into the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
				
				new_ray.shift();
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				if (!new_it.did_hit() && !new_ray.get_medium()) 
				{	
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
					new_ray.shift();
				}
				current_it = new_it;
			}
		}
	}
	// Finally, if we're getting multiple time samples, we push their reconstructed 
	// radiance into the output smples.
	/*if( multiple_time_samples )
		for( int i=0; i < nb_time_samples; ++i )
			if( !time_samples[i].isZero() )
				samples.push_back( RadianceSample( time_samples[i], 
												  i*delta_time,  
												  0, ray_dist));*/
}




template<class TermCriteria, int D>
void VolumePathTracing<TermCriteria,D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
    Real offset_time =  film?film->get_time_offset():0.,
                        film_time_length = film?film->get_time_length():0., 
                        film_exp = film->get_exposure_time();

	Real max_time = (nb_time_samples-1) * delta_time;
	Real inv_delta_time = 1./delta_time, ray_dist(-1.);
	bool multiple_time_samples = (nb_time_samples > 1);
	bool sample_surface = false;
	
	std::vector<Spectrum> time_samples;
	if( multiple_time_samples )
		time_samples.resize(nb_time_samples);
	
    if (!film) printf("WARNING: Film should be != NULL\n");
    
	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		//NEW RAY
		Intersection<D> current_it(it);
		Ray<D> new_ray(current_it.get_ray().get_origin(), current_it.get_ray().get_direction(), current_it.get_ray().get_level(), DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
		new_ray.set_parameter(current_it.get_ray().get_parameter());
		
		//SAMPLE MEAN FREE PATH
		Real epsilon1, epsilon2;
		Spectrum u_t, u_s, scat_albedo;
		Real avg_scat_albedo, pdf_light1, pdf_light2;
		Real sampled_mfp = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
		
		Spectrum attenuation(1.);
		Real time = 0.;
		ray_dist = -1.;
		
		//TIME SAMPLING VARIABLES
		Real t_path = 0.0f, T_FILM = film_time_length + offset_time, T_MAX, mis_surface = 1.;
		while(1)
		{
			sample_surface = false;
			if(new_ray.get_medium())
			{	
				u_t = new_ray.get_medium()->get_extinction(new_ray.get_origin());
				u_s = new_ray.get_medium()->get_scattering(new_ray.get_origin());
				scat_albedo = u_s/u_t;
				avg_scat_albedo = scat_albedo.avg();
			}
			else
			{
				u_t = Spectrum(0.);
				u_s = Spectrum(0.);
				sampled_mfp = std::numeric_limits<Real>::infinity();
			}
			
			//only in the beginning of the path starting in the camera ray			
			if(nb_time_samples == 0) //RADIANCE SAMPLING
			{
				epsilon1 = RNG::StdRandom::get_real();
				while (epsilon1 < 1e-5) 
					epsilon1 = RNG::StdRandom::get_real();

				sampled_mfp = -logf(epsilon1)/u_t.avg();
				
				if (ray_dist < 0.) ray_dist = std::min(sampled_mfp, new_ray.get_parameter());	
				if (!current_it.did_hit() || sampled_mfp < new_ray.get_parameter()) //SCATTERING EVENT
				{
					time += this->m_world->time_of_flight(sampled_mfp);
                    
                    if (time > T_FILM)
						break;
					VectorN<D> p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
				
					//Get light source 
					const LightSource<D> *l = this->m_world->sample_light(pdf_light1);			
					LightSample<D> ls;
					Real pdf_light2;
					if(l->sample(p, ls, pdf_light2))
					{	
						Spectrum radiance = scat_albedo*attenuation*exp(-u_t*ls.dist)*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())/(pdf_light1*pdf_light2)*m_inv_incoming_samples;
					
						//sample_time takes into account the length of the path followed by several bounces and the time the light takes to travel until the scattering event
						Real sample_time = time + this->m_world->time_of_flight(ls.dist)/new_ray.get_refraction_index() + ls.instant;
						samples.push_back( RadianceSample(radiance, 
																sample_time,  
																new_ray.get_level(),
																this->m_world->time_of_flight(ray_dist)));
					}
				
					Real pdf;
					/*if ( m_termCriteria(current_it, pdf))
						break;
						attenuation /= pdf;
						*/
					epsilon2 = RNG::StdRandom::get_real();

					/*epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
					while (epsilon2 < 1e-5)
						epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
					*/
					if (epsilon2 > avg_scat_albedo) 
						break;
				
					//Random walk's next step
					// Get sampled direction plus pdf, and update attenuation
					VectorN<D> new_omega_i;
					Medium<D> *m = new_ray.get_medium();
				
					attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
					attenuation /= pdf*avg_scat_albedo;
				
					new_ray = Ray<D>(p, new_omega_i, new_ray.get_level()+1, DEFAULT_REFRACTION_INDEX, m);
				
					Intersection<D> new_it;
					this->m_world->first_intersection(new_ray, new_it);
				
					current_it = new_it;
				}
				else 
					sample_surface = true;
			}
			else //TIME SAMPLING
			{
				TimeSampler<D> time_sampler (this->m_world, new_ray.get_refraction_index());
				Real ta, tb;
				const LightSource<D> *light = this->m_world->sample_light(pdf_light1);
					
				VectorN<D> l, omega;
				LightSample<D> ls;
				Real pdf_light2;

				//compute ta
				light->sample(new_ray.get_origin(), ls, pdf_light2);
				l = ls.pos - new_ray.get_origin();
				omega = new_ray.get_direction();
                ta = max(offset_time - t_path, this->m_world->time_of_flight(l.length())/new_ray.get_refraction_index());
				//ta = this->m_world->time_of_flight(ls.dist)/new_ray.get_refraction_index();
				//ta = this->m_world->time_of_flight(l.length())/new_ray.get_refraction_index();

				//compute tb
				if (new_ray.did_hit())
				{
					//UPDATE MAXIMUM TIME WITH THE MINIMUM BETWEEN THE TIME RANGE OF THE FILM AND THE TIME TO THE CLOSEST SURFACE
					light->sample(new_ray.get_position(), ls, pdf_light2);
					T_MAX = std::min(this->m_world->time_of_flight(new_ray.get_parameter()+ls.dist)/new_ray.get_refraction_index(), T_FILM - t_path);
					//T_MAX = this->m_world->time_of_flight(ls.dist + (ls.dist*ls.dist - 2*ls.dist*dot(new_ray.get_direction(), l)));
					//T_MAX = std::min(ta + this->m_world->time_of_flight(new_ray.get_parameter()), T_FILM - t_path);
				}
				else T_MAX = T_FILM - t_path;

				tb = T_MAX;
				
				Real sampled_dist, pdf_time;
				VectorN<D> p;
				if (tb < ta || time > T_FILM)
				{
					break;
				}	
#if MIS
				//SAMPLE MEAN FREE PATH EITHER WAY
				epsilon1 = RNG::StdRandom::get_real();
				while (epsilon1 < 1e-5) 
					epsilon1 = RNG::StdRandom::get_real();

				sampled_mfp = -logf(epsilon1)/u_t.avg();

				//MIS ==> DECIDE IF USING TIME SAMPLING OR MFP SAMPLING
				Real prob_mfp, prob_ts, sum_prob_pdf, w;
				prob_mfp = 1- pow(time/T_FILM,1);
				//if (prob_mfp < 0) printf("asdf\n");
				prob_ts = 1. - prob_mfp;
				
				if (RNG::StdRandom::get_real() < prob_mfp)
				{
					if (ray_dist < 0.) ray_dist = std::min(sampled_mfp, new_ray.get_parameter());	
					//CHECK IF REACHING SURFACE
					
					if (!current_it.did_hit() || sampled_mfp < new_ray.get_parameter()) //SCATTERING EVENT
					{
						pdf_time = time_sampler.pdf(sampled_mfp,ta,tb,l,omega);
						sum_prob_pdf = prob_mfp*exp(-u_t.avg()*sampled_mfp) + prob_ts*pdf_time;
						w = prob_mfp/sum_prob_pdf;
						p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
						if( light->sample(p, ls, pdf_light2) )
						{
							Spectrum radiance = w*attenuation*exp(-u_t*(ls.dist+sampled_mfp))*scat_albedo*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
							radiance /= pdf_light1*pdf_light2*prob_mfp;
							Real sample_time = time + this->m_world->time_of_flight(sampled_mfp+ls.dist);// + ls.instant;
							//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

							//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
							samples.push_back( RadianceSample(radiance, 
																sample_time,  
																new_ray.get_level(),
																this->m_world->time_of_flight(ray_dist)));
						}
					}
					else{ 
						/*pdf_time = time_sampler.pdf(sampled_mfp,ta,tb,l,omega);
						//printf("%f\n", pdf_time); fgetc(stdin);
						sum_prob_pdf = prob_mfp*exp(-u_t.avg()*sampled_mfp) + prob_ts*pdf_time;
						w = prob_mfp/sum_prob_pdf;*/
						mis_surface = 1./prob_mfp;
						sample_surface = true;
					}
				}
				else
				{
					
					time_sampler(l,omega,ta,tb,sampled_dist,pdf_time);
					sum_prob_pdf = prob_mfp*exp(-u_t.avg()*sampled_dist) + prob_ts*pdf_time;
					w = prob_ts/sum_prob_pdf;
					mis_surface = w/prob_ts;
					if (sampled_dist < new_ray.get_parameter())
					{	
						//printf("w = %f\n", w);
						p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
						if( light->sample(p, ls, pdf_light2) )
						{
							Spectrum radiance = w*attenuation*exp(-u_t*(ls.dist+sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
							radiance /= pdf_light1*pdf_light2*prob_ts;
							Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
							//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

							//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
							samples.push_back( RadianceSample(radiance, 
																sample_time,  
																new_ray.get_level(),
																this->m_world->time_of_flight((ray_dist < 0.)?sampled_dist:ray_dist)));
						}
						if (ray_dist < 0.) ray_dist = sampled_mfp;
					}
					else break;
				}
#elif STRAT
				if (ta > tb) break;

				Real t_exp = this->m_world->time_of_flight(1/u_t.avg());
				//Real t_exp = delta_time/static_cast<Real>(nb_time_samples);
				Real ta_i, tb_i;
				ta_i = ta;
				tb_i = ta + t_exp;
				//printf("\n");
				while (tb_i < tb + t_exp){
					//printf("tb = %.3f, tb_i = %.3f\n",tb,tb_i);
					time_sampler(l,omega,ta_i, tb_i, sampled_dist, pdf_time);
						
					p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
					if( light->sample(p, ls, pdf_light2))
					{
						//printf("sampled_dist = %f, ta = %f, tb = %f\n", sampled_dist + ls.dist,ta,tb);
						Spectrum radiance = attenuation*exp(-u_t*(ls.dist + sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
						radiance /= pdf_light1*pdf_light2*pdf_time;
						//if (radiance.avg() < 1e-5) printf("attenuation = %f, u_s = %f, pdfs: %f %f %f\n", attenuation.avg(), u_s.avg(), pdf_light1, pdf_light2, pdf_time);

						Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
						//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

						//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
						samples.push_back( RadianceSample(radiance, 
															sample_time,  
															new_ray.get_level(),
															this->m_world->time_of_flight(ray_dist)));
					}
					ta_i += t_exp;
					tb_i += t_exp;
				}
#elif TIME
				//if (ta > tb) break;
				//sample per frame

				time_sampler(l,omega,ta,tb,sampled_dist,pdf_time);
				p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
				if( light->sample(p, ls, pdf_light2))
				{
					Spectrum radiance = attenuation*exp(-u_t*(ls.dist + sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
					radiance /= pdf_light1*pdf_light2*pdf_time;
					//if (radiance.avg() < 1e-5) printf("attenuation = %f, u_s = %f, pdfs: %f %f %f\n", attenuation.avg(), u_s.avg(), pdf_light1, pdf_light2, pdf_time);

					Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
					//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

					//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
					samples.push_back( RadianceSample(radiance, 
														sample_time,  
														new_ray.get_level(),
														this->m_world->time_of_flight(ray_dist)));
                    
				}

                epsilon1 = RNG::StdRandom::get_real();
				while (epsilon1 < 1e-5) 
					epsilon1 = RNG::StdRandom::get_real();

				sampled_mfp = -logf(epsilon1)/u_t.avg();
                if (current_it.did_hit() && sampled_mfp > new_ray.get_parameter())
                    sample_surface = true;
                
#else				//LINE TO LINE
				if (ta > tb) break;
				//printf("%f\n",ta);
				// Line-to-line sampling, j subscript is for light path, i subscript is for camera path
				Real pdf_tj, rj, rj_max, ri, ri_max, ra, rb;
				VectorN<D> wj, xj;
				// 1. Sample light direction
				Real pdf_wj;
				LightSample<D> ls_j;
				light->sample(ls_j, pdf_wj);
				wj = ls_j.dir;
				xj = ls_j.pos;
				rj_max =(T_MAX - t_path)/this->m_world->time_of_flight(1.);

				//2. compute interval [ra,rb] of first sampled path (we choose camera path)
				ri_max = (T_MAX -t_path)/this->m_world->time_of_flight(1.);
				ra = time_sampler.rx(ta,rj_max, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
				rb = time_sampler.rx(tb,0.0f, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
					
				//Real taux = time_sampler.t(ra, rj_max, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
				//if (fabs(taux - ta) > 100) break;// printf("rj = %f ra = %f, taux, ta = %f %f\n",rj, ra, taux, ta);
				//taux = time_sampler.t(rb, 0, new_ray.get_origin(), xj, new_ray.get_direction(), wj);
					
				ra = max(ra,0.0f);
				rb = min(rb, ri_max);
				//3. Uniformly sample a position inside [ra,rb]
				ri = RNG::StdRandom::get_real()*(rb-ra);
				//printf("%f\n",ri);
				Real pdf_rba = 1./(rb-ra);
				VectorN<D> x_ri = new_ray.get_origin() + new_ray.get_direction()*ri;
					
				time_sampler(x_ri - xj, wj, ta, tb, rj, pdf_tj);
				VectorN<D> x_rj = xj + wj*rj;
				Real rij = (x_ri - x_rj).length();
				Real att_dist = 1./((D==3)?(rij*rij):(rij));
				//Assuming visibility, this should be changed
				Spectrum radiance = attenuation*att_dist*exp(-u_t*(ri+rj+rij))
									*u_s*u_s
									*ls_j.irradiance
									*this->m_world->get_medium()->f(x_rj, ls_j.dir, x_ri - x_rj)
									*this->m_world->get_medium()->f(x_ri, x_rj - x_ri, -new_ray.get_direction())
									*m_inv_incoming_samples;	
				radiance /= pdf_light1*pdf_light2*pdf_tj*pdf_rba*pdf_wj;
				//if (radiance.avg() < 1e-5) printf("attenuation = %f, u_s = %f, pdfs: %f %f %f\n", attenuation.avg(), u_s.avg(), pdf_light1, pdf_light2, pdf_time);

				Real sample_time = time + this->m_world->time_of_flight(ri+rj+rij);// + ls.instant;
				//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

				//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
				samples.push_back( RadianceSample(radiance, 
													sample_time,  
													new_ray.get_level(),
													this->m_world->time_of_flight(ray_dist)));
#endif
				if(!sample_surface){
					p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp;
					time += this->m_world->time_of_flight(sampled_mfp)/new_ray.get_refraction_index();
					//time += this->m_world->time_of_flight(sampled_dist)/new_ray.get_refraction_index();
					t_path = time;
					Real pdf;
					//Random walk's next step
					// Get sampled direction plus pdf, and update attenuation
					VectorN<D> new_omega_i;
					Medium<D> *m = new_ray.get_medium();
				
					attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
					attenuation /= pdf;

					new_ray = Ray<D>(p, new_omega_i, new_ray.get_level()+1, DEFAULT_REFRACTION_INDEX, m);
				
					Intersection<D> new_it;
					this->m_world->first_intersection(new_ray, new_it);
				
					current_it = new_it;
				}
			}
			if (sample_surface) //SURFACE HIT
			{
				time += this->m_world->time_of_flight(current_it.get_ray().get_parameter());
				
				const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);
				
				//Get direct contribution
				//if( !current_it.delta_material() )
				LightSample<D> light_sample;
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
					Spectrum radiance = mis_surface*light_sample.irradiance*current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())
					*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation*exp(-u_t*light_sample.dist+current_it.get_ray().get_parameter())/(pdf_light1*pdf_light2)*m_inv_incoming_samples;
					
					Real sample_time = time + this->m_world->time_of_flight(light_sample.dist)/current_it.get_ray().get_refraction_index() + light_sample.instant;

					/*if ( multiple_time_samples )
					{
						printf("What?\n");
						if( sample_time < max_time )
						{	
							Real time_sample = sample_time * inv_delta_time;
							int min_sample = static_cast<int>(time_sample);
							
							time_samples[min_sample] += radiance*(1.-time_sample+min_sample);
							time_samples[min_sample+1] += radiance*(time_sample-min_sample);
						}
					}
					else*/
						samples.push_back( RadianceSample(radiance, 
														  sample_time,  
														  current_it.get_ray().get_level(),
														  this->m_world->time_of_flight(ray_dist)));
				}
				
				//Test termination condition
				Real pdf;
				if ( m_termCriteria(current_it, pdf))
					break;
				attenuation /= pdf;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				
				if(it.material()->is_type(Reflectance::DELTA_TRANSMISSION))
				{
					current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
					attenuation *= dot_abs(current_it.get_normal(), new_omega_i);
				}
				else{
					attenuation *= current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
					attenuation *= dot_abs(current_it.get_normal(), new_omega_i)/pdf; 
				}
				//attenuation*=mis_surface;
				// Throw ray and update current_it
				if (dot(current_it.get_normal(), new_omega_i) >= 0.) //new ray outside the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
				else //new ray into the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
				
				new_ray.shift();
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				
				if (!new_it.did_hit() && !new_ray.get_medium()) //singularities when new out direction and normal are perpendicular
				{	
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
					new_ray.shift();
				}
				
				current_it = new_it;
			}
		}
	}	
	
	// Finally, if we're getting multiple time samples, we push their reconstructed 
	// radiance into the output smples.
	/*if( multiple_time_samples )
		for( int i=0; i < nb_time_samples; ++i )
			if( !time_samples[i].isZero() )
				samples.push_back( RadianceSample( time_samples[i], 
												  i*delta_time,  
												  0, ray_dist));*/
}

template<class TermCriteria, int D>
void VolumePathTracing<TermCriteria,D>::operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	Real max_time = (nb_time_samples-1) * delta_time;
	Real inv_delta_time = 1./delta_time, ray_dist(-1.);
	bool multiple_time_samples = (nb_time_samples > 1);
	
	std::vector<Spectrum> time_samples;
	if( multiple_time_samples )
		time_samples.resize(nb_time_samples);
	
	//srand(time(NULL));
	for(int sample = 0; sample < this->m_incoming_samples; sample++)
	{
		//NEW RAY
		Real pdf;
		VectorN<D> new_omega_i;
		//Sampling::circular(new_omega_i,pdf);
		Sampling::spherical(new_omega_i,pdf);
		Ray<D> new_ray(p, new_omega_i, 1, DEFAULT_REFRACTION_INDEX, this->m_world->get_medium()); //ASSUMES POINT P IS INSIDE THE MEDIUM
		
		Intersection<D> current_it;
		this->m_world->first_intersection(new_ray, current_it);
		if (current_it.did_hit() && dot(current_it.get_normal(), -new_ray.get_direction()) < 0.0f)
			break; //point is inside object, careful with open geometry!
		//SAMPLE MEAN FREE PATH
		Real epsilon1, epsilon2;
		Spectrum u_t, u_s, scat_albedo;
		Real avg_scat_albedo, pdf_light1, pdf_light2;
		Real sampled_mfp = std::numeric_limits<Real>::infinity(); //Sampled mean free path in free-space
		
		Spectrum attenuation(1./pdf);
		Real time = 0.;
		ray_dist = -1.;

		//TIME SAMPLING VARIABLES
		Real t_path = 0.0f, T_MAX = delta_time;
		while(1)
		{
			if(new_ray.get_medium())
			{	
				u_t = new_ray.get_medium()->get_extinction(new_ray.get_origin());
				u_s = new_ray.get_medium()->get_scattering(new_ray.get_origin());
				
				epsilon1 = RNG::StdRandom::get_real();
				while (epsilon1 < 1e-5) 
					epsilon1 = RNG::StdRandom::get_real();

				sampled_mfp = -logf(epsilon1)/u_t.avg();
				scat_albedo = u_s/u_t;
				avg_scat_albedo = scat_albedo.avg();
			}
			else
			{
				u_t = Spectrum(0.);
				u_s = Spectrum(0.);
				sampled_mfp = std::numeric_limits<Real>::infinity();
			}
			
			//only in the beginning of the path starting in the camera ray
			if (ray_dist < 0.) ray_dist = std::min(sampled_mfp, new_ray.get_parameter());
			
			if (!current_it.did_hit() || sampled_mfp < new_ray.get_parameter()) //SCATTERING EVENT
			{
				if (nb_time_samples == 1) //MFP SAMPLING
				{
					time += this->m_world->time_of_flight(sampled_mfp);

					VectorN<D> p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp; //point of scattering event
					
					//Get light source 
					const LightSource<D> *l = this->m_world->sample_light(pdf_light1);			
					LightSample<D> ls;
					Real pdf_light2;
					if(l->sample(p, ls, pdf_light2))
					{	
						Spectrum radiance = (u_s)*attenuation*exp(-u_t*ls.dist)*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())/(pdf_light1*pdf_light2)*m_inv_incoming_samples;
						
						//sample_time takes into account the length of the path followed by several bounces and the time the light takes to travel until the scattering event
						Real sample_time = time + this->m_world->time_of_flight(ls.dist)/new_ray.get_refraction_index() + ls.instant;
						
						if ( multiple_time_samples )
						{
							if( sample_time < max_time )
							{	
								Real time_sample = sample_time * inv_delta_time;
								int min_sample = static_cast<int>(time_sample);
								
								time_samples[min_sample] += radiance*(1.-time_sample+min_sample);
								time_samples[min_sample+1] += radiance*(time_sample-min_sample);
							}
						}
						else
							samples.push_back( RadianceSample(radiance, 
															  sample_time,  
															  new_ray.get_level(),
															  this->m_world->time_of_flight(ray_dist)));
					}
					
					/*if ( m_termCriteria(current_it, pdf))
					 break;
					 attenuation /= pdf;
					 */
					epsilon2 = RNG::StdRandom::get_real();
					/*epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
					while (epsilon1 < 1e-5)
						epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
					*/
					if (epsilon2 > avg_scat_albedo) 
						break;
					
					//Random walk's next step
					// Get sampled direction plus pdf, and update attenuation
					VectorN<D> new_omega_i;
					Medium<D> *m = new_ray.get_medium();
					
					attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
					attenuation /= pdf*avg_scat_albedo;
					
					new_ray = Ray<D>(p, new_omega_i,new_ray.get_level()+1, DEFAULT_REFRACTION_INDEX, m);
					
					Intersection<D> new_it;
					this->m_world->first_intersection(new_ray, new_it);
					
					current_it = new_it;
				}
				else //TIME SAMPLING
				{
					TimeSampler<D> time_sampler (this->m_world, new_ray.get_refraction_index());
					Real ta, tb;
					const LightSource<D> *light = this->m_world->sample_light(pdf_light1);
					
					VectorN<D> l, omega;
					
					LightSample<D> ls;
					Real pdf_light2;
					light->sample(new_ray.get_origin(), ls, pdf_light2);
					
					l = ls.pos - new_ray.get_origin();
					omega = new_ray.get_direction();
					ta = this->m_world->time_of_flight(l.length())/new_ray.get_refraction_index();
					tb = T_MAX - t_path;
					if (tb < ta) break;
					//printf("ta, tb = %f %f\n", ta, tb);
					Real sampled_dist, pdf_time;
					time_sampler(l,omega,ta,tb,sampled_dist,pdf_time);

					VectorN<D> p = new_ray.get_origin() + new_ray.get_direction()*sampled_dist; //point of scattering event
					if( light->sample(p, ls, pdf_light2))
					{
						Spectrum radiance = attenuation*exp(-u_t*(ls.dist + sampled_dist))*u_s*ls.irradiance*this->m_world->get_medium()->f(p,ls.dir,-new_ray.get_direction())*m_inv_incoming_samples;	
						radiance /= pdf_light1*pdf_light2*pdf_time;
						Real sample_time = time + this->m_world->time_of_flight(sampled_dist+ls.dist);// + ls.instant;
						//printf("ls.dist = %f, time = %f\n", ls.dist, sample_time);

						//if (sampled_dist <= 0.0) printf("sample_distance = %f\n", sampled_dist);
						samples.push_back( RadianceSample(radiance, 
															sample_time,  
															new_ray.get_level(),
															this->m_world->time_of_flight(ray_dist)));
					}
					//p = new_ray.get_origin() + new_ray.get_direction()*sampled_mfp;
					//time += this->m_world->time_of_flight(sampled_dist)/new_ray.get_refraction_index();
					time += this->m_world->time_of_flight(sampled_mfp)/new_ray.get_refraction_index();
					t_path = time;
					Real pdf;
					//Random walk's next step
					// Get sampled direction plus pdf, and update attenuation
					VectorN<D> new_omega_i;
					Medium<D> *m = new_ray.get_medium();
				
					attenuation *= scat_albedo*m->sample_direction(p, new_ray.get_direction(), new_omega_i, pdf);
					attenuation /= pdf;
				
					new_ray = Ray<D>(p, new_omega_i, new_ray.get_level()+1, DEFAULT_REFRACTION_INDEX, m);
				
					Intersection<D> new_it;
					this->m_world->first_intersection(new_ray, new_it);
				
					current_it = new_it;
				}
			}
			else //SURFACE HIT
			{
				time += this->m_world->time_of_flight(current_it.get_ray().get_parameter());
				
				const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);
				
				//Get direct contribution
				//if( !current_it.delta_material() )
				LightSample<D> light_sample;
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
					Spectrum radiance = light_sample.irradiance*current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())
					*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation*exp(-u_t*light_sample.dist)/(pdf_light1*pdf_light2)*m_inv_incoming_samples;
					
					Real sample_time = time + this->m_world->time_of_flight(light_sample.dist)/current_it.get_ray().get_refraction_index() + light_sample.instant;

					if ( multiple_time_samples )
					{
						if( sample_time < max_time )
						{	
							Real time_sample = sample_time * inv_delta_time;
							int min_sample = static_cast<int>(time_sample);
							
							time_samples[min_sample] += radiance*(1.-time_sample+min_sample);
							time_samples[min_sample+1] += radiance*(time_sample-min_sample);
						}
					}
					else
						samples.push_back( RadianceSample(radiance, 
														  sample_time,  
														  current_it.get_ray().get_level(),
														  this->m_world->time_of_flight(ray_dist)));
				}
				
				//Test termination condition
				Real pdf;
				if ( m_termCriteria(current_it, pdf))
					break;
				attenuation /= pdf;
				
				//Random walk's next step
				// Get sampled direction plus pdf, and update attenuation
				VectorN<D> new_omega_i;
				
				if(current_it.material()->is_type(Reflectance::DELTA_TRANSMISSION))
				{
					current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
					attenuation *= dot_abs(current_it.get_normal(), new_omega_i);
				}
				else{
					attenuation *= current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, current_it.get_normal(),current_it.get_uv(), pdf );
					attenuation *= dot_abs(current_it.get_normal(), new_omega_i)/pdf; 
				}
				
				// Throw ray and update current_it
				if (dot(current_it.get_normal(), new_omega_i) >= 0.) //new ray outside the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
				else //new ray into the object
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
				
				new_ray.shift();
				
				Intersection<D> new_it;
				this->m_world->first_intersection(new_ray, new_it);
				
				if (!new_it.did_hit() && !new_ray.get_medium()) //singularities when new out direction and normal are perpendicular
				{	
					new_ray = Ray<D>(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index(), this->m_world->get_medium());
					new_ray.shift();
				}
				
				current_it = new_it;
			}
		}
	}
	
	// Finally, if we're getting multiple time samples, we push their reconstructed 
	// radiance into the output smples.
	/*if( multiple_time_samples )
		for( int i=0; i < nb_time_samples; ++i )
			if( !time_samples[i].isZero() )
				samples.push_back( RadianceSample( time_samples[i], 
												  i*delta_time,  
												  0, ray_dist));*/
}

#endif //_PATH_TRACING_INTEGRATOR_H_
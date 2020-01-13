#ifndef _PATH_TRACING_INTEGRATOR_H_
#define _PATH_TRACING_INTEGRATOR_H_

#include "Integrator/SurfaceIntegrator.h"
#include "RayTracing/World.h"

template<class TermCriteria, int D>
class PathTracing: public SurfaceIntegrator<D>
{
	TermCriteria m_termCriteria;
public:
	PathTracing(World<D> *w, const TermCriteria &termination,  unsigned int incoming_samples = 1);
	~PathTracing();

	virtual Spectrum operator()(const Intersection<D> &it) const;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
};


template<class TermCriteria, int D>
PathTracing<TermCriteria,D>::PathTracing(World<D> *w, const TermCriteria &_termCriteria, unsigned int incoming_samples)
	:SurfaceIntegrator<D>(w,incoming_samples), m_termCriteria(_termCriteria)
{}

template<class TermCriteria, int D>
PathTracing<TermCriteria,D>::~PathTracing()
{}

template<class TermCriteria, int D>
Spectrum PathTracing<TermCriteria,D>::operator()(const Intersection<D> &it) const
{
	Spectrum reflected_radiance;	

	// Just to obtain depth maps... need to be redone!
	//return Vector3(it.material()->f(Vector3(1,0,0),Vector3(-1,0,0),Vector3(1,0,0),Vector2(0)));
	//return Vector3(it.get_ray().get_parameter());
	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		bool end = false;
		Spectrum attenuation(1.);
		Intersection<D> current_it(it);
		while(1)
		{
			if( current_it.material()->is_type(Reflectance::EMISSIVE) )
			{
				Spectrum radiance = current_it.material()->emit_light(-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())*attenuation;
				reflected_radiance += radiance;
			}

			//Get light source 
			Real pdf_light1;
			const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);

			//Get direct contribution
			//if( !current_it.delta_material() )
			//if( current_it.get_ray().get_level() > 0 )
			{
				LightSample<D> light_sample;
				Real pdf_light2;
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
					/*Spectrum radiance = light_sample.irradiance*
						current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),it.get_normal(),it.get_uv())
						*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation
						/(pdf_light1*pdf_light2);*/
					
					Spectrum radiance = light_sample.irradiance*current_it.material()->f(light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv())
						*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation
						/(pdf_light1*pdf_light2);

					reflected_radiance += radiance;
				}
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
			Ray<D> new_ray(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
			new_ray.shift();

			Intersection<D> new_it;
			this->m_world->first_intersection(new_ray, new_it);
			
			if( !new_it.did_hit()) 
				break;

			current_it = new_it;
		}
	}
	return reflected_radiance*this->m_inv_incoming_samples;
}

template<class TermCriteria, int D>
void PathTracing<TermCriteria,D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	Real max_time = (nb_time_samples-1) * delta_time;
	Real inv_delta_time = 1./delta_time;
	bool multiple_time_samples = (nb_time_samples > 1);
	
	std::vector<Spectrum> time_samples;
	if( multiple_time_samples )
		time_samples.resize(nb_time_samples);


	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		bool end = false;
		
		Spectrum attenuation(1.);
		//Real time = 0; //picoseconds(10^-12)
		Real time = this->m_world->time_of_flight(it.get_ray().get_parameter());
		Intersection<D> current_it(it);
		while(1)
		{
			// If hitting an emitting surface
			if( current_it.material()->is_type(Reflectance::EMISSIVE) )
			{
				Real delta_time_material, pdf_emission;
				Spectrum radiance = current_it.material()->emit_light(-current_it.get_ray().get_direction(),
							current_it.get_normal(),current_it.get_uv(), delta_time_material, pdf_emission)*attenuation;

				radiance /= (pdf_emission)*this->m_inv_incoming_samples;
				Real sample_time = time + delta_time_material;
					

				samples.push_back( RadianceSample(radiance, 
						sample_time,  
						current_it.get_ray().get_level(),
						it.get_ray().get_parameter()));
			}

			//Get light source 
			Real pdf_light1;
			const LightSource<D> *l = this->m_world->sample_light(current_it, pdf_light1);

			//Get direct contribution
			{
				LightSample<D> light_sample;

				Real pdf_light2; 
				if( l->sample( current_it, light_sample, pdf_light2 ) )
				{
									
					// If any light is reaching the sample point, then accomulate radiance
					Real delta_time_material, pdf_material;
					Spectrum radiance = 
						current_it.material()->f(light_sample.irradiance, light_sample.dir,-current_it.get_ray().get_direction(),current_it.get_normal(),current_it.get_uv(),
												 delta_time_material, pdf_material)
						*dot_clamped(-light_sample.dir,current_it.get_normal())*attenuation;

					//if(pdf_material < 1.)
					//	printf("PDF: %f, Radiance: [%f, %f, %f] Time: %f Delta Time: %f\n", pdf_material, radiance[0],radiance[1],radiance[2],time + this->m_world->time_of_flight(light_sample.dist)/current_it.get_ray().get_refraction_index() + light_sample.instant, delta_time_material);
	
					radiance /=(pdf_light1*pdf_light2*pdf_material)*this->m_inv_incoming_samples;
					
					Real sample_time = time + this->m_world->time_of_flight(light_sample.dist)*current_it.get_ray().get_refraction_index() + light_sample.instant + delta_time_material;
					
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
							it.get_ray().get_parameter()));
				}
			}	
			
			//Test termination condition
			Real pdf;
			if ( m_termCriteria(current_it, pdf))
				break;
			attenuation /= pdf;

			//Random walk's next step
			// Get sampled direction plus pdf, and update attenuation
			/*VectorN<D> new_omega_i; 
			Real delta_time_material;
			attenuation *= current_it.material()->sample_direction(current_it.get_ray().get_direction(), new_omega_i, 
																   current_it.get_normal(),current_it.get_uv(), delta_time_material, pdf );
			attenuation *= dot_abs(current_it.get_normal(), new_omega_i)/pdf; 


			// Throw ray and update current_it
			Ray<D> new_ray(current_it.get_position(), new_omega_i, current_it.get_ray().get_level()+1,  current_it.get_ray().get_refraction_index());
			new_ray.shift();*/
			VectorN<D> new_omega_i; 
			Ray<D> new_ray;

			Real delta_time_material;
			attenuation *= current_it.material()->sample_outgoing_ray(current_it, new_ray, delta_time_material, pdf );
			if( !current_it.material()->is_type(Reflectance::DELTA))
				attenuation *= dot_abs(current_it.get_normal(), new_ray.get_direction())/pdf; 
			else
				attenuation *= 1./pdf;

			// Throw ray and update current_it
			new_ray.shift();
			

			Intersection<D> new_it;
			this->m_world->first_intersection(new_ray, new_it);
			
			if( !new_it.did_hit() )
				break;
			
			// Update time, bounce, omega_i and current_it (take care with index of refraction!)
			time += this->m_world->time_of_flight(new_it.get_ray().get_parameter())*new_ray.get_refraction_index() + delta_time_material;
			current_it = new_it;
		}
	}
	
	// Finally, if we're getting multiple time samples, we push their reconstructed 
	// radiance into the output smples.
	if( multiple_time_samples )
		for( int i=0; i < nb_time_samples; ++i )
			if( !time_samples[i].isZero() )
				samples.push_back( RadianceSample( time_samples[i], 
							i*delta_time,  
							0, it.get_ray().get_parameter()));
}


#endif //_PATH_TRACING_INTEGRATOR_H_
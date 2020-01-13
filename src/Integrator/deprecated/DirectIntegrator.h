#ifndef _DIRECT_INTEGRATOR_H_
#define _DIRECT_INTEGRATOR_H_

#include "Integrator/SurfaceIntegrator.h"
#include "RayTracing/World.h"

/** Direct illumination integrator based on Turner Whitted's
	recursive ray-tracing lighting. It iteratively shades the
	delta recursive BRDFs (transmissive and specular), while
	computes the illumination from direct lights only.
	As opposed to Whitted's, the integrator supports area light,
	by (stochastically) sampling them. 
	Whitted, T. 1980. "An improved illumination model for shaded
	display" In Communications of the ACM.
	http://dl.acm.org/citation.cfm?id=358882	*/
template<int D>
class DirectIntegrator: public SurfaceIntegrator<D>
{
public:
	DirectIntegrator(World<D> *w):SurfaceIntegrator<D>(w){}
	~DirectIntegrator(){}

	virtual Spectrum operator()(const Intersection<D> &it) const;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
};

template<int D>
Spectrum DirectIntegrator<D>::operator()(const Intersection<D> &it) const
{
	Spectrum reflected_radiance;

	//Get direct contribution
	if( !it.material()->is_type(Reflectance::DELTA) )
	{
		//Get light source 
		for( int light = 0; light < Integrator<D>::m_world->nb_lights(); ++light)
		{
			const LightSource<D> *l = this->m_world->light(light);
			LightSample<D> light_sample;
			
			Real pdf_light; 
			if( l->sample( it, light_sample, pdf_light ) )
			{
				Spectrum radiance = light_sample.irradiance*
					it.material()->f(light_sample.dir,-it.get_ray().get_direction(),it.get_normal(),it.get_uv())
					*dot_clamped(-light_sample.dir,it.get_normal())
					/pdf_light;

				reflected_radiance += radiance;
			}
		}
	}
	return reflected_radiance;
}

template<int D>
void DirectIntegrator<D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
	//Get direct contribution
	if( !it.material()->is_type(Reflectance::DELTA) )
	{
		//Get light source 
		for( int light = 0; light < this->m_world->nb_lights(); ++light)
		{
			const LightSource<D> *l = this->m_world->light(light);
			LightSample<D> light_sample;
			
			Real pdf_light; 
			if( l->sample( it, light_sample, pdf_light ) )
			{
				Spectrum radiance = light_sample.irradiance*
					it.material()->f(light_sample.dir,-it.get_ray().get_direction(),it.get_normal(),it.get_uv())
						*dot_clamped(-light_sample.dir,it.get_normal())
						/pdf_light;

				samples.push_back( RadianceSample(radiance, 
					this->m_world->time_of_flight(light_sample.dist)/it.get_ray().get_refraction_index(),  
					it.get_ray().get_level()));
			}
		}
	}
}


#endif //_DIRECT_INTEGRATOR_H_
#ifndef _SINGLE_SCATTERING_INTEGRATOR_H_
#define _SINGLE_SCATTERING_INTEGRATOR_H_

#include "Integrator/VolumeIntegrator.h"
#include "RayTracing/World.h"
#include "Sampling/SamplingTechniques.h"
#include "Path.h"

//=============================================================
template<class TermCriteria, int D>
class SingleScatteringIntegrator: public VolumeIntegrator<D>
{
public:
	enum VertexSampling{ MeanFreePath, Equiangular, TimeSampling};

private:
	VertexSampling m_vertex_sampling;
	Real m_beta_MIS;

    Film *film;
	TermCriteria m_termCriteria;
	unsigned int m_incoming_samples;
	Real m_inv_incoming_samples;
	FILE *f_log;


	Real sample_distance(const Ray<D> &r, const LightSample<D> light, Real &p, bool &media_interaction )const;
	bool connect_vertex_with_light( const typename Path<D>::Vertex* v, const LightSource<D> *light, Spectrum &f, Real &p, Real *time=0)const;

public:
	SingleScatteringIntegrator(World<D> *w, const TermCriteria &termination,  unsigned int incoming_samples = 1, FILE *_f_log = NULL, Film *_film = NULL)
	:VolumeIntegrator<D>(w), m_termCriteria(termination), m_incoming_samples(incoming_samples),
	m_inv_incoming_samples(1./static_cast<Real> (m_incoming_samples)), f_log(_f_log), film(_film),
	m_vertex_sampling(VertexSampling::MeanFreePath)
	{}

	~SingleScatteringIntegrator()
	{}
	void preprocess() 
	{
		if (f_log)
		{
			fprintf(f_log, "Ray Trace Single Scattering\n");
			fprintf(f_log, "=======================\n");
			fprintf(f_log, "# Samples: %d\n", m_incoming_samples);
			fprintf(f_log, "\n");
		}
	}


	
	
	void set_sampling_technique(unsigned int tech)
	{
		switch(tech)
		{
		case 0:
			m_vertex_sampling = VertexSampling::MeanFreePath;
			break;
		case 1:
			m_vertex_sampling = VertexSampling::Equiangular;
			break;
		case 2:
			m_vertex_sampling = VertexSampling::TimeSampling;
		}
	}
	
	void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const;

	void operator()(const Ray<D> &r0, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Spectrum operator()(const Ray<D> &r0) const;// {return Spectrum(0.);};
	
	void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Spectrum operator()(const Intersection<D> &it) const;
}; //SingleScatteringIntegrator

//=============================================================================
template<class TermCriteria, int D>
Real SingleScatteringIntegrator<TermCriteria,D>::sample_distance(const Ray<D> &r, const LightSample<D> light, Real &p, bool &media_interaction )const
{
	Real t=-1.;
	media_interaction = true;

	//printf("\n- Sampling\n");
	switch(m_vertex_sampling)
	{
	case VertexSampling::TimeSampling:
		{
			Real ta = m_world->time_of_flight((light.pos-r.get_origin()).length())*r.get_refraction_index();
			Real tb = film->get_time_resolution()*film->get_exposure_time()+film->get_time_offset();
			Probability::TimeLineToPoint<D> st(ta, tb, r.get_direction(), light.pos-r.get_origin(), r.get_refraction_index(), 1./m_world->time_of_flight(1));
			while(t<0) t = st.cdf_inv(RNG::StdRandom::get_real());
			p = st.pdf(t);

			if( r.did_hit() && (r.get_parameter() < t ) )
			{
				t = r.get_parameter();
				p = st.cdf(t);
				media_interaction = false;
			}
		}
		break;
	case VertexSampling::Equiangular:
		{
			Probability::Equiangular<D> st( r.get_direction(),light.pos-r.get_origin());
			while(t<0) t = st.cdf_inv(RNG::StdRandom::get_real());
			t = st.cdf_inv(RNG::StdRandom::get_real());
			p = st.pdf(t);

			if( r.did_hit() && (r.get_parameter() < t ) )
			{
				t = r.get_parameter();
				p = st.cdf(t);
				media_interaction = false;
			}
		}
		break;
	case VertexSampling::MeanFreePath:
		{
			Real u_t = r.get_medium()->get_extinction(r.get_origin()).avg();
			t = -logf(RNG::StdRandom::get_real())/u_t;
			p = exp(-t*u_t)*u_t;
			// Check if before ray intersection point
			if( r.did_hit() && (r.get_parameter() < t ) )
			{
				t = r.get_parameter();
				p = exp(-r.get_parameter()*u_t);
				media_interaction = false; // Check
			}	

			
		}
		break;
	}

	//printf("--- %d -> T: %f; P: %f; Media: %d\n", m_vertex_sampling,t,p,media_interaction);
	return t;
}



template<class TermCriteria, int D>
bool SingleScatteringIntegrator<TermCriteria,D>::connect_vertex_with_light( const typename Path<D>::Vertex* v, const LightSource<D> *light, Spectrum &f, Real &p, Real *time)const
{
	LightSample<D> light_sample;
	Spectrum fe;
	Real pl = 1.;
	Real pe = 1., te = 0.;
	
	f = Spectrum(0.);
	p = 1.;

	if (light->sample( v->get_vertex_position(), light_sample, pl ))
	{	
		if( time )
		{
			v->compute_scattering(-light_sample.dir, fe, te, pe);
			(*time)= m_world->time_of_flight(light_sample.dist)*m_world->get_ior() + te;

		}
		else
			v->compute_scattering(-light_sample.dir, fe );

		Spectrum att(1.);
		if( m_world->get_medium()) 
			att = exp(-m_world->get_medium()->get_extinction(v->get_vertex_position())*light_sample.dist);
		
		p = pe*pl;
		
		
		f = fe*light_sample.irradiance*att;
		//f = att;
		return true;
	}	
	else
	{
		f = Spectrum();
		p = 1.;
		if( time )
			(*time)= m_world->time_of_flight(light_sample.dist)*m_world->get_ior();
		
		return false;
	}
}

//=============================================================================
// Computes the light transport from an initial camera-sampled ray r.
template<class TermCriteria, int D>
Spectrum SingleScatteringIntegrator<TermCriteria,D>::operator()(const Ray<D> &r0) const
{
	Spectrum reflected_radiance(0.);
	
	Intersection<D> curr_it; 
	Ray<D> curr_ray(r0.get_origin(),r0.get_direction(), false, r0.get_level(), r0.get_refraction_index(), r0.get_medium()); 
	this->m_world->first_intersection(curr_ray, curr_it);

	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		// Sample light
		LightSample<D> light_sample;
		const LightSource<D> *light;
		Real pl1, pl2, pl3;
		light = m_world->sample_light(pl1);
		light->sample( light_sample, pl2 );
		
		// Sample eye vertex
		typename Path<D>::Vertex* v;
		Real pv, r;
		bool media_interaction = true;
		r = sample_distance(curr_ray, light_sample, pv, media_interaction);

		VectorN<D> position;
		if( media_interaction )
		{	
			position = curr_ray.get_origin() + curr_ray.get_direction()*r;
			v = new typename Path<D>::VolumeVertex(position, curr_ray.get_direction(), curr_ray.get_medium(), 
							Spectrum(1), pv, pv);
		}
		else
		{
			position = curr_it.get_position();
			v = new typename Path<D>::SurfaceVertex(position, curr_ray.get_direction(), curr_it.get_normal(), curr_it.get_uv(), curr_it.material(),
							Spectrum(1), pv, pv);
		}

		// Connect with light and attenuate
		Spectrum f;
		if( connect_vertex_with_light( v, light, f, pl3 ) )
		{	
			Spectrum T = exp(-r*curr_ray.get_medium()->get_extinction(position));
			f *= T;
			reflected_radiance += f/(pl1*pv*pl3);
		}	


		// Remove vertex
		delete(v);
	}
	return reflected_radiance*m_inv_incoming_samples;
}


template<class TermCriteria, int D>
Spectrum SingleScatteringIntegrator<TermCriteria,D>::operator()(const Intersection<D> &it) const
{
	return (*this)(it.get_ray());
	
}

template<class TermCriteria, int D>
void SingleScatteringIntegrator<TermCriteria,D>::operator()(const Ray<D> &r0, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	Spectrum reflected_radiance(0.);
	
	Intersection<D> curr_it; 
	Ray<D> curr_ray(r0.get_origin(),r0.get_direction(), false, r0.get_level(), r0.get_refraction_index(), r0.get_medium()); 
	this->m_world->first_intersection(curr_ray, curr_it);

	for(int sample = 0; sample < this->m_incoming_samples; ++sample )
	{
		Real t;
		// Sample light
		LightSample<D> light_sample;
		const LightSource<D> *light;
		Real pl1, pl2, pl3;
		light = m_world->sample_light(pl1);
		light->sample( light_sample, pl2 );
		
		// Sample eye vertex
		typename Path<D>::Vertex* v;
		Real pv, r;
		bool media_interaction = true;
		r = sample_distance(curr_ray, light_sample, pv, media_interaction);

		VectorN<D> position;
		if( media_interaction )
		{	

			position = curr_ray.get_origin() + curr_ray.get_direction()*r;
			v = new typename Path<D>::VolumeVertex(position, curr_ray.get_direction(), curr_ray.get_medium(), 
							Spectrum(1), pv, pv);
		}
		else
		{
			position = curr_it.get_position();
			v = new typename Path<D>::SurfaceVertex(position, curr_ray.get_direction(), curr_it.get_normal(), curr_it.get_uv(), curr_it.material(),
							Spectrum(1), pv, pv);
		}

		// Connect with light and attenuate
		Spectrum f;
		if( connect_vertex_with_light( v, light, f, pl3, &t ) )
		{	
			Spectrum T = exp(-r*curr_ray.get_medium()->get_extinction(position));
			f *= T;
			Real p = (pl1*pl2*pv*pl3*m_incoming_samples);
			t += m_world->time_of_flight(r)*m_world->get_ior();

			samples.push_back( RadianceSample(f/p, t)); 
		}
		else
		{
			t += m_world->time_of_flight(r)*m_world->get_ior();
			samples.push_back( RadianceSample(Spectrum(0.), t)); 
		}
		// Remove vertex
		delete(v);		
	}
}




template<class TermCriteria, int D>
void SingleScatteringIntegrator<TermCriteria,D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	(*this)(it.get_ray(), samples, delta_time, nb_time_samples);
}

template<class TermCriteria, int D>
void SingleScatteringIntegrator<TermCriteria,D>::operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{}

#endif //_BIDIRECTIONAL_PATH_TRACING_INTEGRATOR_H_
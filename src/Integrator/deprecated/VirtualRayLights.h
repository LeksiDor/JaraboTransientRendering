/*
 *  VirtualRayLights.h
 *  transient
 */

#ifndef _VRL_H_
#define _VRL_H_

#include "Photon.h"
#include "BeamTracing.h"
#include "RayTracing/Intersection.h"
#include "RayTracing/Ray.h"
#include "Color/Spectrum.h"
#include "DataStructures/KDTree.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
//#include <tgmath.h>
#include <list>
#include <stdio.h>

#ifndef _SAMPLES_VRL_
#define _SAMPLES_VRL_ 2
#endif

template<class TermCriteria, int D, BeamBlurRegion BR>
class VirtualRayLights: public BeamTracing<TermCriteria,D,BR>
{
	std::vector< CosineLightSource<D> > m_vpls;

	KDTree<Photon<D>, D> m_caustics_tree;
	BVH<Beam<D,BR>, D> *m_caustics_bvh;

	std::vector<Beam<D,BR> > m_global_vector, m_caustics_vector;
	Real m_max_dist, m_inv_max_area, m_inv_global_region_measure, m_inv_caustics_region_measure; 
	Real m_global_beam_radius, m_caustic_beam_radius; 
	int m_nb_photons;
	bool m_final_gathering, frozen, rm_single_scattering;
	int m_nb_samples_vrl;

	typedef typename KDTree<Photon<D>, D>::Node KDTreeNode;
	FILE *f_log;
public:
	VirtualRayLights(World<D> *w, Real max_dist, Real global_beam_radius, Real caustic_beam_radius, int nb_photons, bool final_gathering,
				int nb_global_photons = std::numeric_limits<int>::max(),
				int nb_caustics_photons = std::numeric_limits<int>::max(),
				int nb_global_beams = std::numeric_limits<int>::max(),
				int nb_caustics_beams = std::numeric_limits<int>::max(),
				int max_nb_global_shots = 5000000, int _b_subdivide = 0, 
				bool _rm_single_scattering = false, Scattering::Level _scat_level = Scattering::ALL, FILE *_f_log = NULL,
				unsigned int incoming_samples = 1)
	:BeamTracing<TermCriteria,D,BR>(w, global_beam_radius, caustic_beam_radius, incoming_samples,
									nb_global_photons, nb_caustics_photons,
									nb_global_beams, nb_caustics_beams,
									max_nb_global_shots, _b_subdivide, false, _scat_level),
	m_max_dist(max_dist), m_inv_max_area(1./(M_PI*m_max_dist*m_max_dist)), m_nb_photons(nb_photons), m_final_gathering(final_gathering), 
	m_global_beam_radius(global_beam_radius), m_caustic_beam_radius(caustic_beam_radius),
	m_inv_global_region_measure((D==2)?(0.5/global_beam_radius):((BR==BR1D)?(0.5/global_beam_radius):(1./(M_PI*global_beam_radius*global_beam_radius)))),  
	m_inv_caustics_region_measure((D==2)?(0.5/caustic_beam_radius):((BR==BR1D)?(0.5/caustic_beam_radius):(1./(M_PI*caustic_beam_radius*caustic_beam_radius)))), 
	rm_single_scattering(_rm_single_scattering), frozen(false), m_caustics_bvh(NULL), f_log(_f_log), m_nb_samples_vrl(_SAMPLES_VRL_)
	{	}
	
	VirtualRayLights(World<D> *w, const std::string &namefile, 
				int _b_subdivide = 0, 
				unsigned int incoming_samples = 1)
	:BeamTracing<TermCriteria,D,BR>(w,1.,1.,incoming_samples,std::numeric_limits<int>::max(),
									std::numeric_limits<int>::max(),5000000, _b_subdivide, false, Scattering::ALL),
	m_max_dist(10.),m_inv_max_area(1./(M_PI*m_max_dist*m_max_dist)), m_nb_photons(100), m_final_gathering(false),
	m_inv_global_region_measure((D==2)?(0.5):((BR==BR1D)?(0.5):(1./M_PI))),  
	m_inv_caustics_region_measure((D==2)?(0.5):((BR==BR1D)?(0.5):(1./M_PI))),
	rm_single_scattering(false),  frozen(false), m_nb_samples_vrl(_SAMPLES_VRL_)
	{
		load(namefile);
		frozen = true;
	}
	~VirtualRayLights(){}
	
	// Functions to set irradiance estimation parameters
	void set_max_dist(const Real max_dist){ m_max_dist = max_dist; m_inv_max_area = 1./(M_PI*m_max_dist*m_max_dist);}
	void set_nb_photons(const int nb_photons){ m_nb_photons = nb_photons;}
	
	
	void preprocess();
	
	void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const;

	void operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
	Spectrum operator()(const Ray<D> &r) const 
	{
		std::list<RadianceSample> samples;
		(*this)(r,samples);

		std::list<RadianceSample>::const_iterator s;
		Spectrum R(0.);

		for(s=samples.begin(); s!=samples.end(); s++)
			R += s->radiance;

		return R;
	}
	
	void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
	Spectrum operator()(const Intersection<D> &it) const 
	{
		std::list<RadianceSample> samples;
		(*this)(it,samples);

		std::list<RadianceSample>::const_iterator s;
		Spectrum R(0.);

		for(s=samples.begin(); s!=samples.end(); s++)
			R += s->radiance;

		return R;			
	};

private:
	// Auxiliary functions needed to sample the VRL
	Real distance( const Ray<D> &r1, const Ray<D> &r2, Real &t1, Real &t2 )const
	{
		Real a = dot(r1.get_direction(),r1.get_direction());
		Real b = dot(r1.get_direction(),r2.get_direction());
		Real c = dot(r2.get_direction(),r2.get_direction());
		Real d = dot(r1.get_direction(),r1.get_origin());
		Real e = dot(r2.get_direction(),r2.get_origin());
		
		t1 = (b*e-c*d)/(a*c-b*b);
		t2 = (a*e-b*d)/(a*c-b*b);
		
		VectorN<D> h =	r1.get_origin() + t1*r1.get_direction() -
						r2.get_origin() + t2*r2.get_direction();

		return h.length();
	}

	// Function that computes distance from the INFINITE line where
	// the ray lays on to the point.
	Real distance( const Ray<D> &r, const VectorN<D> &p, Real &t)const
	{
		VectorN<D> ap = p-r.get_origin();
		
		t = dot(ap,r.get_direction());
		
		return (r.get_origin()+r.get_direction()*t-p).length();
	}
	
	// Sampling functions
	Real sample_vrl( const Real A0, const Real A1, const Real h, 
					 const Real sinthita, Real &pdf)const;
	Real sample_camera_ray( const Real B0, const Real B1, const Real h, Real &pdf)const;
	void compute_single_scattering(	const Ray<D> &r, std::list<RadianceSample> &samples, 
									const Real delta_time, const unsigned int nb_time_samples)const;
	void print_info()const;

}; //VirtualRayLights

template<class TermCriteria, int D, BeamBlurRegion BR>
void VirtualRayLights<TermCriteria,D,BR>::print_info()const
{
	fprintf(f_log, "\n");
	fprintf(f_log, "PHOTON BEAMS\n");
	fprintf(f_log, "============\n");
	fprintf(f_log, "Global beam radius: %f\n", m_global_beam_radius);
	fprintf(f_log, "Caustic beam radius: %f\n", m_caustic_beam_radius);
	fprintf(f_log, "Global photons: %d\n", global_photons.size());
	fprintf(f_log, "Caustic photons: %d\n", caustic_photons.size());
	fprintf(f_log, "Global beams: %d\n", global_beams.size());
	fprintf(f_log, "Caustic beams: %d\n", caustic_beams.size());
	
	if (D == 2 || BR == BR1D)
		fprintf(f_log, "Blur: 1D\n");	
	else
		fprintf(f_log, "Blur: 2D\n");
	
	if (rm_single_scattering)
		fprintf(f_log, "Ray marching single scattering: TRUE\n");	
	else
		fprintf(f_log, "Ray marching single scattering: FALSE\n");	
	fprintf(f_log, "\n");
	printf("Total Photons: %7d global photons, %7d caustic photons\n", global_photons.size(), caustic_photons.size());
	printf("Total Beams: %7d global beams, %7d caustic beams\n", global_beams.size(), caustic_beams.size());
}

// Preprocess functions
template<class TermCriteria, int D, BeamBlurRegion BR>
void VirtualRayLights<TermCriteria,D,BR>::preprocess()
{
	if (frozen) return;
	//Clear trees
	m_caustics_tree.clear();
	
	//Trace photons
	std::list< Photon<D> > global_photons, caustic_photons;
	std::vector< Beam<D,BR> > global_beams, caustic_beams;

	Scattering::Level aux_scat_level = scat_level; scat_level = Scattering::ALL;

	global_beams.clear();
	caustic_beams.clear();
	global_photons.clear();
	caustic_photons.clear();
	trace_beams(global_beams, caustic_beams, global_photons, caustic_photons);
	
	scat_level = aux_scat_level;


	m_global_vector = global_beams;
	m_caustics_vector = caustic_beams;
	
	// Store photons into KDTrees and balance KDTrees
	typename std::list<Photon<D> >::iterator it;

	// Store caustic beams, so we can apply photon beams
	if (caustic_beams.size()) 	
	{
		m_caustics_bvh = new BVH<Beam<D,BR>, D>(caustic_beams);
	}

	// Need to store VPLs/VSLs	
	if( global_photons.size() )
	{
	}

	if( caustic_photons.size() )
	{
		FILE *fphotons = fopen("../../ws/caustic_photons.txt", "w+");
		for( it = caustic_photons.begin(); it != caustic_photons.end(); it ++)
		{		
			m_caustics_tree.store(std::vector<Real>(it->m_position.m_data, 
													it->m_position.m_data+D), *it);
			//fprintf(fphotons, "[%f %f %f] power: %f %f %f\n", it->m_position[0], it->m_position[1], it->m_position[2], it->m_power[0], it->m_power[1], it->m_power[2]);
		}
		fclose(fphotons);
		caustic_photons.clear();
		m_caustics_tree.balance();
	}
	// Log:

	fprintf(f_log, "\n");
	fprintf(f_log, "VIRTUAL RAYLIGHTS\n");
	fprintf(f_log, "============\n");
	fprintf(f_log, "Global beam radius: %f\n", m_global_beam_radius);
	fprintf(f_log, "Caustic beam radius: %f\n", m_caustic_beam_radius);
	fprintf(f_log, "Global photons: %d\n", global_photons.size());
	fprintf(f_log, "Caustic photons: %d\n", caustic_photons.size());
	fprintf(f_log, "Global beams: %d\n", global_beams.size());
	fprintf(f_log, "Caustic beams: %d\n", caustic_beams.size());
	
	if (D == 2 || BR == BR1D)
		fprintf(f_log, "Blur: 1D\n");	
	else
		fprintf(f_log, "Blur: 2D\n");
	
	if (rm_single_scattering)
		fprintf(f_log, "Ray marching single scattering: TRUE\n");	
	else
		fprintf(f_log, "Ray marching single scattering: FALSE\n");	
	fprintf(f_log, "\n");
	printf("Total Photons: %7d global photons, %7d caustic photons\n", global_photons.size(), caustic_photons.size());
	printf("Total Beams: %7d global beams, %7d caustic beams\n", global_beams.size(), caustic_beams.size());
}


// Integration functions
namespace VRLaux
{
	Real rng()
	{
		Real epsilon = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		while( epsilon < 0.0 || epsilon >= 1.)
			epsilon = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);

		return epsilon;
	}

	Real asinh(Real value)
	{   
		if(value>0)
			return log(value + sqrt(value * value + 1));
		else
			return -log(-value + sqrt(value * value + 1));
	}
};

using namespace VRLaux;

// Sample the Virtual Ray Light assuming an isotropic medium. Note that 
// this function samples according to the attenuation IN 3D.
template<class TermCriteria, int D, BeamBlurRegion BR>
Real VirtualRayLights<TermCriteria,D,BR>::sample_vrl( const Real A0, const Real A1, 
													  const Real h, const Real sinthita, 
													  Real &pdf)const
{
	// Compute vt according to Eq. 13 in [1].
	Real epsilon = rng();
	Real vt = h*sinh(A1*epsilon+A0*(1.f-epsilon))/sinthita;

	// Compute the pdf according to Eq.9 in [1]. 
	// This is equivalent to (Eq.10)/(Eq.11) in [1]. Note that
	// the terms have been already simplified.
	pdf = sinthita / (sqrtf(h*h+vt*vt*sinthita*sinthita)*(A1-A0));

	return vt;
}

// Sample the camera ray assuming an isotropic medium, using Kulla and Fajardo's
// equiangular sampling [2]. Note that this function samples according to the 
// attenuation IN 3D.
template<class TermCriteria, int D, BeamBlurRegion BR>
Real VirtualRayLights<TermCriteria,D,BR>::sample_camera_ray( const Real B0, const Real B1, 
													  const Real h, Real &pdf)const
{
	// Compute vt according to Eq. 14 in [1], which bases
	// on Eq. 6 in [2].
	Real epsilon = rng();
	Real ut = h*tan(B1*epsilon+B0*(1.f-epsilon));

	// Compute the pdf according to a modified version of 
	// Eq.5 in [2].
	pdf = h / ((B1-B0)*(h*h + ut*ut));

	return ut;
}

//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
// Should be changed for a more general implementation or move it to another volume integrator	
template<class TermCriteria, int D, BeamBlurRegion BR>
void VirtualRayLights<TermCriteria,D,BR>::compute_single_scattering(const Ray<D> &r, std::list<RadianceSample> &samples, 
													 const Real delta_time, const unsigned int nb_time_samples ) const
{
	Spectrum u_t = (m_world->get_medium()?m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
	u_s=(m_world->get_medium()?m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
	Real ior = r.get_refraction_index();
	
	if (m_world->get_medium())
	{
		Real stepSize = 0.001, step, dist = min(r.get_parameter(),100.0f), t0 = 0., pdf_light1;
		unsigned int nSamples = ceil(dist/stepSize), s; //max boundary of scene to 100.
		step = dist/nSamples;
		VectorN<D> p = r.get_origin(), pPrev;
		Spectrum Tr(1.);
		t0 += 0.01*stepSize;
		for (s = 0; s < nSamples; s++) 
		{
			pPrev = p;
			p = r.get_origin() + r.get_direction()*t0;
			const LightSource<D> *l = m_world->sample_light(pdf_light1);
			LightSample<D> lsample;
			Real pdf_light2;
			Tr *= exp(-m_world->get_medium()->get_extinction(VectorN<D>(0.))*step);
			if(l->sample(p, lsample, pdf_light2)) 
			{
				samples.push_back(RadianceSample(step*u_s*Tr*m_world->get_medium()->f(p,lsample.dir,-r.get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2),
												 m_world->time_of_flight(lsample.dist + t0)*ior, r.get_level(), this->m_world->time_of_flight(t0)*ior));
			}
			t0 += step; 
		}
	}
}

template<class TermCriteria, int D, BeamBlurRegion BR>
void VirtualRayLights<TermCriteria,D,BR>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
#ifdef _USED_
	Spectrum reflected_radiance(0.), attenuation(1.), global_beams_contribution(0.), caustic_beams_contribution(0.);
	Spectrum total_single_scattering(0.), single_scattering(0.);
	std::vector<const Beam<D,BR>*> intersected_global_beams, intersected_caustics_beams;
	//Get the contribution of all the beams until the first intersection
	typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;
	
	if (rm_single_scattering && this->is_level(Scattering::SINGLE) ){
		//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
		//Should be changed for a more general implementation or move it to another volume integrator
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
				 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
		Real ior = it.get_ray().get_refraction_index();
		if(this->m_world->get_medium())
		{	
			Real step, stepSize = 0.001, dist = min(it.get_ray().get_parameter(),100.0f), t0 = 0., pdf_light1, pdf_light2;
			unsigned int nSamples = ceil(dist/stepSize), s; //max boundary of scene to 100.
			step = dist/nSamples;
			VectorN<D> p = it.get_ray().get_origin(), pPrev;
			Spectrum Tr0(1.); //extinction along the ray marching ray, while Tr is the accumulated extinction of the previous paths
			t0 += 0.01*stepSize;
			single_scattering = Spectrum(0.);
			for (s = 0; s < nSamples; s++) {
				pPrev = p;
				p = it.get_ray().get_origin() + it.get_ray().get_direction()*t0;
				const LightSource<D> *l = this->m_world->sample_light(pdf_light1);
				LightSample<D> lsample;
				Tr0 = exp(-u_t*t0);
				if(l->sample(p, lsample, pdf_light2)) 
				{
					samples.push_back(RadianceSample(step*u_s*Tr0*this->m_world->get_medium()->f(p,lsample.dir,-it.get_ray().get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2),
													 this->m_world->time_of_flight(lsample.dist + t0)*ior, it.get_ray().get_level(), this->m_world->time_of_flight(t0)*ior));
				}
				t0 += step; 
			}
		}
	}
	
	//GLOBAL
	if (m_global_bvh)
	{
		m_global_bvh->intersect_all(it.get_ray(), it.get_ray().get_parameter(), intersected_global_beams);
		for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
		{
			RadianceSample sample;
			(*ptr_beam_it)->contribution(it.get_ray(), sample);
			sample.radiance *= (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_global_region_measure;
			if (BR != BR2D) 
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, delta_time);
		}
	}
	//CAUSTICS
	if (m_caustics_bvh)
	{
		m_caustics_bvh->intersect_all(it.get_ray(), it.get_ray().get_parameter(), intersected_caustics_beams);
		for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
		{	
			RadianceSample sample;
			(*ptr_beam_it)->contribution(it.get_ray(), sample);
			sample.radiance *= (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_caustics_region_measure;
			if (BR != BR2D) 
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, delta_time);
		}
		
	}
	intersected_global_beams.clear();
	intersected_caustics_beams.clear();
	
	if (!it.did_hit())
		return;
	
	Spectrum Tr(1.);
	Medium<D> *m = this->m_world->get_medium(); //assuming infinite homogeneous medium // it.get_ray().get_medium();
	//extinction along the camera ray
	if(m) 
	{
		Tr = exp(-m->get_extinction(it.get_ray().get_origin())*it.get_ray().get_parameter());
		//printf("Tr = %f %f %f\n", Tr[0], Tr[1], Tr[2]);
	}
	
	//Handle delta materials (e.g. perfect reflection/refraction)
	Intersection<D> shaded_it(it);
	//camera time
	Real time =this->m_world->time_of_flight(it.get_ray().get_parameter())*it.get_ray().get_refraction_index();
	
	while( shaded_it.material()->is_type(Reflectance::DELTA) )
	{
		Real pdf;
		//Check termination condition...
		if ( m_termCriteria(shaded_it, pdf))
			return;	
		
		//Get output direction, plus the attenuation
		VectorN<D> new_omega_i;
		Ray<D> new_ray;
		
		attenuation *= 1. - shaded_it.material()->get_absorption(it.get_uv());
		shaded_it.material()->sample_outgoing_ray(shaded_it, new_ray, pdf);
		
		bool outside = dot(new_ray.get_direction(),shaded_it.get_normal()) > 0.;
		
		Intersection<D> new_it;
		//Get new intersection
		this->m_world->first_intersection(new_ray, new_it);
		
		if (rm_single_scattering && outside && this->is_level(Scattering::SINGLE) ){
			//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
			//Should be changed for a more general implementation or move it to another volume integrator
			Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
					 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
			Real ior = new_ray.get_refraction_index();
			if(this->m_world->get_medium())
			{	
				Real step, stepSize = 0.001, dist = min(new_ray.get_parameter(),100.0f), t0 = 0., pdf_light1, pdf_light2;
				unsigned int nSamples = ceil(dist/stepSize), s; //max boundary of scene to 100.
				step = dist/nSamples;
				VectorN<D> p = new_ray.get_origin(), pPrev;
				Spectrum Tr0(1.);
				t0 += 0.01*stepSize;
				single_scattering = Spectrum(0.);
				for (s = 0; s < nSamples; s++) {
					pPrev = p;
					p = new_ray.get_origin() + new_ray.get_direction()*t0;
					const LightSource<D> *l = this->m_world->sample_light(pdf_light1);
					LightSample<D> lsample;
					Tr0 = exp(-this->m_world->get_medium()->get_extinction(VectorN<D>(0.))*t0);
					if(l->sample(p, lsample, pdf_light2)) 
					{
						samples.push_back(RadianceSample(Tr*attenuation*step*u_s*Tr0*this->m_world->get_medium()->f(p,lsample.dir,-new_ray.get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2),
														 time + this->m_world->time_of_flight(lsample.dist + t0)*ior, new_ray.get_level(), this->m_world->time_of_flight(t0)*ior + time));
					}
					t0 += step; 
				}
			}
		}
		
		if (m_global_bvh && outside)
		{
			m_global_bvh->intersect_all(new_ray, new_ray.get_parameter(), intersected_global_beams);
			for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
			{
				RadianceSample sample;
				(*ptr_beam_it)->contribution(new_ray, sample);
				sample.distance += time; //add the distance to camera for unwarping
				sample.time += time;
				sample.radiance *= Tr*attenuation*(*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_global_region_measure;
				if (BR != BR2D) 
					samples.push_back(sample);
				else if (sample.radiance.avg() != 0.)
					time_blur(sample, samples, delta_time);
			}
			
		}
		//CAUSTICS
		if (m_caustics_bvh && outside)
		{
			m_caustics_bvh->intersect_all(new_ray, new_ray.get_parameter(), intersected_caustics_beams);
			for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
			{	
				RadianceSample sample;
				(*ptr_beam_it)->contribution(new_ray, sample);
				sample.distance += time; //add the distance to camera for unwarping
				sample.time += time;
				sample.radiance *= Tr*attenuation*(*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_caustics_region_measure;
				if (BR != BR2D) 
					samples.push_back(sample);
				else if (sample.radiance.avg() != 0.)
					time_blur(sample, samples, delta_time);
			}
			
		}
		intersected_global_beams.clear();
		intersected_caustics_beams.clear();		
		
		if( !new_it.did_hit() )
			return;
		
		if (outside) //ray outside of object
			Tr *= exp(-m->get_extinction(new_ray.get_origin())*new_ray.get_parameter());
		
		shaded_it = new_it;
		time += this->m_world->time_of_flight(new_ray.get_parameter())*new_ray.get_refraction_index();
	}
	
	if (rm_single_scattering && shaded_it.did_hit() && this->is_level(Scattering::SINGLE) )
	{
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
				 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
		Real ior = shaded_it.get_ray().get_refraction_index();
		
		//DIRECT LIGHT ON SURFACE
		Real pdf_light1, pdf_light2;
		const LightSource<D> *l = this->m_world->sample_light(shaded_it, pdf_light1);
		LightSample<D> light_sample;
		if( l->sample(shaded_it, light_sample, pdf_light2 ) )
		{
			Spectrum radiance = light_sample.irradiance*shaded_it.material()->f(light_sample.dir,-shaded_it.get_ray().get_direction(),shaded_it.get_normal(),shaded_it.get_uv())
			*dot_clamped(-light_sample.dir,shaded_it.get_normal())*exp(-u_t*light_sample.dist)/(pdf_light1*pdf_light2);
			
			samples.push_back(RadianceSample(Tr*attenuation*radiance, time + this->m_world->time_of_flight(light_sample.dist)*ior, shaded_it.get_ray().get_level(), time));
		}
	}
	
	//	Get all samples where the radiance estimation is needed.
	//	If we are performing final gathering, compute the reflected
	//	samples. Else, just add the current intersection.
	list<RadianceEstimationSample<D> > estimation_samples;
	
	if( m_final_gathering )
	{
		for(unsigned int sample = 0; sample < this->m_incoming_samples; ++sample )
		{
			//Get output direction, plus the attenuation
			Real pdf;
			VectorN<D> new_omega_i;
			Spectrum form_factor(1.);
			
			form_factor = shaded_it.material()->sample_direction(shaded_it.get_ray().get_direction(),
																 new_omega_i, shaded_it.get_normal(), shaded_it.get_uv(), pdf );
			
			form_factor *= dot_abs(shaded_it.get_normal(), new_omega_i)/pdf;
			
			//Create outgoing ray
			Ray<D> new_ray(shaded_it.get_position(), new_omega_i, shaded_it.get_ray().get_level()+1,  shaded_it.get_ray().get_refraction_index());
			new_ray.shift();
			
			
			//Get new intersection
			Intersection<D> new_it;
			this->m_world->first_intersection(new_ray, new_it);
			
			//Again, take care with the index of refraction!!
			if( new_it.did_hit() )
				estimation_samples.push_back(RadianceEstimationSample<D>(form_factor,new_it,
																		 time + this->m_world->time_of_flight(new_it.get_ray().get_parameter()), new_ray.get_level()));
		}
	}
	else
		estimation_samples.push_back(RadianceEstimationSample<D>(Spectrum(1.),shaded_it, time, shaded_it.get_ray().get_level()));
	
	
	//Compute radiance estimation! 
	typename std::list<RadianceEstimationSample<D> >::iterator cit;
	for( cit = estimation_samples.begin(); cit != estimation_samples.end(); cit++ )
	{
		Spectrum sample_radiance(0.);
		VectorN<D> p = cit->it.get_position();
		VectorN<D> d = -cit->it.get_ray().get_direction();
		
		std::vector<const KDTreeNode*> photons;
		typename std::vector<const KDTreeNode*>::iterator ph;
		Real max_distance;
		
		// Compute the contribution of the sample
		Spectrum factor_contribution = attenuation * cit->form_factor * (m_final_gathering? this->m_inv_incoming_samples:1.);
		// Get global contribution
		//	Find closest global photons...
		m_global_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_photons,  photons, max_distance);
		Spectrum fc_global;
		if (D == 3)
			fc_global = factor_contribution / (M_PI*max_distance*max_distance);
		else
			fc_global = factor_contribution / (2*max_distance);
		
		//	... and compute their contribution
		/*if ( multiple_time_samples )
		{
			// First we order all close photons according to their time
			sort(photons.begin(), photons.end());
			
			// And compute their power.
			std::vector<Spectrum> photons_power; photons_power.reserve(m_nb_photons); //Candidate to go outside the loop...
			for( ph = photons.begin(); ph != photons.end(); ph++ )
				photons_power.push_back(cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv())
										* ((*ph)->data().m_power)
										* dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal())
										* fc_global);
			
			
			// For each time sample...
			for( unsigned int ts = 0; ts < nb_time_samples; ++ts )
			{	
				// Get the contribution for the sample using the first technique presented
				// by Cammarano and Jensen [2002].
				
				
			}
		}
		else*/
			for( ph = photons.begin(); ph != photons.end(); ph++ )
				samples.push_back(RadianceSample( Tr*cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv())
												 * ((*ph)->data().m_power)
												 * dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal())
												 * fc_global,
												 cit->time + (*ph)->data().m_time, cit->bounce + (*ph)->data().m_level, cit->time)); 
		
		// Get caustics contribution
		//	Find closest caustic photons...
		photons.clear();
		m_caustics_tree.find(std::vector<Real>(p.m_data,p.m_data+D), m_nb_photons,  photons, max_distance);
		Spectrum fc_caustics;
		if (D == 3)
			fc_caustics = factor_contribution / (M_PI*max_distance*max_distance);
		else
			fc_caustics = factor_contribution / (2*max_distance);

		
		/*if ( multiple_time_samples )
		{
			
		}
		else*/
			//	... and compute their contribution
			for( ph = photons.begin(); ph != photons.end(); ph++ )
				samples.push_back(RadianceSample( Tr*cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv())
												 * ((*ph)->data().m_power)
												 * dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal())
												 * fc_caustics,
												 cit->time + (*ph)->data().m_time, cit->bounce + (*ph)->data().m_level, cit->time)); 
	}
	return;
#endif
}










// Integration functions
template<class TermCriteria, int D, BeamBlurRegion BR>
void VirtualRayLights<TermCriteria,D,BR>::operator()(const Ray<D> &r, std::list<RadianceSample> &samples, 
													 const Real delta_time, const unsigned int nb_time_samples ) const
{
	// Check if there's any participating media at all...
	if(!m_world->get_medium())
		return;

	// We'll be assuming a homogeneous medium
	Spectrum u_t= m_world->get_medium()->get_extinction(VectorN<D>(0.)), 
			 u_s= m_world->get_medium()->get_scattering(VectorN<D>(0.));
	Real ior = r.get_refraction_index();

	// Compute single scattering along the ray
	if (rm_single_scattering && this->is_level(Scattering::SINGLE) ) 
		compute_single_scattering(r, samples, delta_time, nb_time_samples );
	
	
	// Compute radiance due to VRLs
	if(~m_global_vector.empty())
	{	
		for( unsigned int i=0; i<m_global_vector.size(); ++i )
		{
			// Get u_h, v_h, and h (Fig. 3 in [1])
			Real uh, vh;
			Real h = distance( r, m_global_vector[i], uh, vh);
			

			// Compute u0t, u1t, v0t and v1t as the change of
			// variables explained in Sec. 4.1 in [1]).
			Real u0t = -uh, u1t = (r.get_parameter()>m_max_dist?m_max_dist:r.get_parameter())-uh,
				 v0t = -vh,
				 v1t = (m_global_vector[i].get_parameter()>m_max_dist?m_max_dist:m_global_vector[i].get_parameter())-vh;

			// Compute values needed to apply VRL and Camera Ray sampling
			// according to [1] and [2].
			Real costhita = dot(r.get_direction(), m_global_vector[i].get_direction());
			Real sinthita = sqrt(1-costhita*costhita);
			Real Av0 = asinh(v0t/h*sinthita),
				 Av1 = asinh(v1t/h*sinthita);
			Real Bu0 = atan2(u0t,h),
				 Bu1 = atan2(u1t,h);
			
			Spectrum flux = m_global_vector[i].get_photon()->m_power;
			Real t_vrl = m_global_vector[i].get_time_of_flight();
			
			// And then compute the contribution of the VRL by solving
			// Eq. 7 [1] using Eq. 8 [1].
			for( unsigned int s=0; s< m_nb_samples_vrl; ++s )
			{
				Real pdf_v, pdf_u;
				// Sample v according to Eq.13 in [1]. Note that we are 
				// assuming an isotropic medium.
				Real vt = sample_vrl( Av0, Av1, h, sinthita, pdf_v );
				
				// Sample v according to Eq.14 in [1], formulated in [2]. 
				Real ut = sample_camera_ray( Bu0, Bu1, h, pdf_u );

				// Evaluate function applying the integrated term in Eq. 7 in [1].
				// Note that we need to unmake the variable changes performed before
				// (Sec. 4.1 in [1]).
				VectorN<D>	u = r.get_origin()+r.get_direction()*(ut+uh),
							v = m_global_vector[i].get_origin() + 
								m_global_vector[i].get_direction()*(vt+vh);
				
				
				// Compute distance
				Real W = (u-v).length();
				
				VectorN<D> v2u = (u-v).normalize();
				// Evaluate the visibility term
				Ray<D> vray( v, v2u );
				if (m_world->intersects(vray, W+1.e-8))
					continue;

				Spectrum Luv(1.f/(D==3?(W*W):W));
				
				// Compute phase functions
				Luv *= m_world->get_medium()->f( u, v2u,-r.get_direction())*
					   m_world->get_medium()->f( v, m_global_vector[i].get_direction(),-v2u);

				// Compute scattering term (assuming homogeneous medium)
				Luv *= u_s*u_s;

				// Finally, compute extinction (assuming homogeneous medium)
				Real c2u2v2l = (ut+uh)+(vt+vh)+W;
				Luv *= exp(-u_t*c2u2v2l);

				// And store the sample
				samples.push_back(RadianceSample( 
							Luv*flux/(pdf_u*pdf_v*(Real)m_nb_samples_vrl), 
							(m_world->time_of_flight(c2u2v2l)+t_vrl)*ior, 
							r.get_level()+m_global_vector[i].get_level(), 
							m_world->time_of_flight((ut+uh)*ior)));
			}
		}
	}


	// Compute radiance due to caustics. Note that this contribution is computed using photon beams,
	// which are more addecuate to this type of phenomena...
	if (m_caustics_bvh)
	{
		std::vector<const Beam<D,BR>*> intersected_caustics_beams;
		typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;

		m_caustics_bvh->intersect_all(r, r.get_parameter(), intersected_caustics_beams);
		for( ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
		{
			/*RadianceSample sample;
			(*ptr_beam_it)->contribution(r, sample);
			sample.radiance *= (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_caustics_region_measure;
			if (BR != BR2D) 
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, delta_time);*/
		}
	}
}


// Integration functions
template<class TermCriteria, int D, BeamBlurRegion BR>
void VirtualRayLights<TermCriteria,D,BR>::operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
	// Check if there's any participating media at all...
	if(!m_world->get_medium())
		return;

	// We'll be assuming a homogeneous medium
	Spectrum u_t= m_world->get_medium()->get_extinction(VectorN<D>(0.)), 
			 u_s= m_world->get_medium()->get_scattering(VectorN<D>(0.));
	Real ior = m_world->get_ior();

	
	//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
	//Should be changed for a more general implementation or move it to another volume integrator
	
	if (rm_single_scattering && this->is_level(Scattering::SINGLE) ) 
	{
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.));
		if (this->m_world->get_medium())
		{
			Real pdf_light1;
			const LightSource<D> *l = this->m_world->sample_light(pdf_light1);
			LightSample<D> lsample;
			Real pdf_light2;
			if(l->sample(p, lsample, pdf_light2)) 
			{
				samples.push_back(RadianceSample(lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2),
												this->m_world->time_of_flight(lsample.dist), 0, 0.));
			}
		}
	}
	
	// Compute radiance due to VRLs. Note that when computing a flux in a point, we 
	// are in a case similar to the described in Sec.4.3 in [1]. Therefore in this 
	// case we are only sampling the VRL, so we are only solving Eq. 6 [1]. 
	if(~m_global_vector.empty())
	{	
		for( unsigned int i=0; i<m_global_vector.size(); ++i )
		{
			// Get u_h, v_h, and h (Fig. 2 in [2], with the names
			// proposed in [1])
			Real vh;
			Real h = distance( m_global_vector[i], p, vh);
			
			// Compute u0t, u1t, v0t and v1t as the change of
			// variables explained in Sec. 4.1 in [1]) 
			Real v0t = -vh,
				 v1t = (m_global_vector[i].get_parameter()>m_max_dist?m_max_dist:m_global_vector[i].get_parameter())-vh;

			// Compute values needed to apply VRL and Camera Ray sampling
			// according to [1] and [2].
			Real Bu0 = atan2(v0t,h),
				 Bu1 = atan2(v1t,h);
			
			Spectrum flux = m_global_vector[i].get_photon()->m_power;
			Real t_vrl = m_global_vector[i].get_time_of_flight();
			
			// And then compute the contribution of the VRL by solving
			// Eq. 6 [1] using Eq. 8 [1].
			for( unsigned int s=0; s< m_nb_samples_vrl; ++s )
			{
				Real pdf_v;
				// Sample v according to Eq.14 in [1], formulated in [2]. 
				Real vt = sample_camera_ray( Bu0, Bu1, h, pdf_v );

				// Evaluate function applying the integrated term in Eq. 6 in [1].
				// Note that we need to unmake the variable changes performed before
				// (Sec. 4.1 in [1]).
				VectorN<D>	v = m_global_vector[i].get_origin() + 
								m_global_vector[i].get_direction()*(vt+vh);
				
				
				// Compute distance
				Real W = (p-v).length();
				
				VectorN<D> v2p = (p-v).normalize();
				// Evaluate the visibility term
				Ray<D> vray( v, v2p );
				if (m_world->intersects(vray, W+1.e-8))
					continue;

				Spectrum Lpv(1.f/(D==3?(W*W):W));
				
				// Compute phase functions
				Lpv *= m_world->get_medium()->f( v, m_global_vector[i].get_direction(),-v2p);

				// Compute scattering term (assuming homogeneous medium)
				Lpv *= u_s;

				// Finally, compute extinction (assuming homogeneous medium)
				Real p2v2l = (vt+vh)+W;
				Lpv *= exp(-u_t*p2v2l);

				// And store the sample
				samples.push_back(RadianceSample( 
							Lpv*flux/(pdf_v*(Real)m_nb_samples_vrl), 
							(m_world->time_of_flight(p2v2l)+t_vrl)*ior, 
							m_global_vector[i].get_level(), 
							0.));
			}
		}
	}
	
	// Compute radiance due to caustics. Note that this contribution is computed using photon beams,
	// which are more addecuate to this type of phenomena...
	if (m_caustics_bvh)
	{
		std::vector<const Beam<D,BR>*> intersected_caustics_beams;
		typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;

		m_caustics_bvh->intersect_all(p, intersected_caustics_beams);
		for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
		{
			/*RadianceSample sample;
			(*ptr_beam_it)->contribution(p, sample);
			sample.radiance *= (*ptr_beam_it)->power_at(0.)*m_inv_caustics_region_measure;
			if (BR == BR2D) 
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, delta_time);*/
		}
	}
}

#endif //_VRL_H_
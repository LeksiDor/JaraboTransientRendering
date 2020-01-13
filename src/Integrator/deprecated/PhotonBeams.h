/*
 *  PhotonBeams.h
 *  transient
 *
 *  Created by Julio on 04/07/2013.
 *
 */

#ifndef _PHOTON_BEAMS_H_
#define _PHOTON_BEAMS_H_
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
#include <list>
#include <stdio.h>

#ifndef _TEMPORAL_KERNEL_SIZE_ 
// In frames
#define _TEMPORAL_KERNEL_SIZE_ 7 
#endif 

#define FIXTIMEBLUR 10

template<class TermCriteria, int D, BeamBlurRegion BR>
class PhotonBeams: public BeamTracing<TermCriteria,D,BR>
{
	KDTree<Photon<D>, D> m_global_tree, m_caustics_tree;
	BVH<Beam<D,BR>, D> *m_global_bvh, *m_caustics_bvh;
	std::vector<Beam<D,BR> > m_global_vector, m_caustics_vector;
	Real m_max_dist, m_inv_max_area, m_inv_global_region_measure, m_inv_caustics_region_measure; 
	Real m_global_beam_radius, m_caustic_beam_radius; 
	int m_nb_photons;
	bool m_final_gathering, frozen, rm_single_scattering;
	typedef typename KDTree<Photon<D>, D>::Node KDTreeNode;
    const unsigned int time_blur_width;
	FILE *f_log;
    Film *film;
public:
	PhotonBeams(World<D> *w, Real max_dist, Real global_beam_radius, Real caustic_beam_radius, int nb_photons, bool final_gathering,
				int nb_global_photons = std::numeric_limits<int>::max(),
				int nb_caustics_photons = std::numeric_limits<int>::max(),
				int nb_global_beams = std::numeric_limits<int>::max(),
				int nb_caustics_beams = std::numeric_limits<int>::max(),
				int max_nb_global_shots = 5000000, int _b_subdivide = 0, 
				bool _rm_single_scattering = false, Scattering::Level _scat_level = Scattering::ALL, FILE *_f_log = NULL,
				unsigned int incoming_samples = 1, Real TMAX = std::numeric_limits<Real>::infinity(), Film *_film = NULL, unsigned int _time_blur_width = 0)
	:BeamTracing<TermCriteria,D,BR>(w, global_beam_radius, caustic_beam_radius, incoming_samples,
									nb_global_photons, nb_caustics_photons,
									nb_global_beams, nb_caustics_beams,
									max_nb_global_shots, _b_subdivide, _rm_single_scattering, _scat_level, TMAX, _film),
	m_max_dist(max_dist), m_inv_max_area(1./(M_PI*m_max_dist*m_max_dist)), m_nb_photons(nb_photons), m_final_gathering(final_gathering), 
	m_global_beam_radius(global_beam_radius), m_caustic_beam_radius(caustic_beam_radius),
	m_inv_global_region_measure((D==2)?(0.5/global_beam_radius):((BR==BR1D)?(0.5/global_beam_radius):(1./(M_PI*global_beam_radius*global_beam_radius)))),  
	m_inv_caustics_region_measure((D==2)?(0.5/caustic_beam_radius):((BR==BR1D)?(0.5/caustic_beam_radius):(1./(M_PI*caustic_beam_radius*caustic_beam_radius)))), 
	rm_single_scattering(_rm_single_scattering), frozen(false), m_global_bvh(NULL), m_caustics_bvh(NULL), f_log(_f_log), film(_film), time_blur_width(_time_blur_width)
	{}
	
	PhotonBeams(World<D> *w, const std::string &namefile, 
				int _b_subdivide = 0, 
				unsigned int incoming_samples = 1)
	:BeamTracing<TermCriteria,D,BR>(w,1.,1.,incoming_samples,std::numeric_limits<int>::max(),
									std::numeric_limits<int>::max(),5000000, _b_subdivide, false, Scattering::ALL),
	m_max_dist(10.),m_inv_max_area(1./(M_PI*m_max_dist*m_max_dist)), m_nb_photons(100), m_final_gathering(false),
	m_inv_global_region_measure((D==2)?(0.5):((BR==BR1D)?(0.5):(1./M_PI))),  
	m_inv_caustics_region_measure((D==2)?(0.5):((BR==BR1D)?(0.5):(1./M_PI))),
	rm_single_scattering(false),  frozen(false)
	{
		load(namefile);
		frozen = true;
	}
	~PhotonBeams(){}
	
	// Functions to set irradiance estimation parameters
	void set_max_dist(const Real max_dist){ m_max_dist = max_dist; m_inv_max_area = 1./(M_PI*m_max_dist*m_max_dist);}
	void set_nb_photons(const int nb_photons){ m_nb_photons = nb_photons;}
	
	
	void preprocess();
	virtual void store(const std::string &namefile)const;
	void load(const std::string &namefile);
	
	void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const;

	void operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Spectrum operator()(const Ray<D> &r) const;// {return Spectrum(0.);};
	
	void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
	Spectrum operator()(const Intersection<D> &it) const;
	
	void time_blur (RadianceSample &sample, std::list<RadianceSample> &samples, Real delta_time) const
	{
		Real N;
        
        Real cam_tof, beam_tof, total_tof;
        
		if(time_blur_width == 0)
        {    
            //BASED ON THE SPATIAL BLUR
            cam_tof = this->m_world->time_of_flight(sample.camera_range);
            //beam_tof = this->m_world->time_of_flight(sample.beam_range);
            total_tof = this->m_world->time_of_flight(sample.beam_range);//cam_tof	+ beam_tof;
            
            N = std::max(1.0f,total_tof/delta_time); //# of pixels covered by the range of the sample
             
            sample.radiance /= N;
            Real f, step = 1./N;
            for (f = -0.5; f <= 0.5; f+=step)
            {
                RadianceSample subsample = RadianceSample(sample);
                subsample.time += f*total_tof;
                //subsample.distance += f*cam_tof;
                samples.push_back(subsample);
            }
        }
        else 
        {
            //FIXED WIDTH: 1, 2 .. N PIXELS
            N = time_blur_width;
            cam_tof = this->m_world->time_of_flight(sample.camera_range);
            total_tof = film->get_exposure_time()*N;
            sample.radiance /= N;
            Real f, step = 1./N;
            for (f = -0.5; f <= 0.5; f+=step)
            {
                RadianceSample subsample = RadianceSample(sample);
                subsample.time += f*total_tof;
                //subsample.distance += f*cam_tof;
                samples.push_back(subsample);
            }
        }
	}
	
}; //PhotonBeams


// Preprocess functions
template<class TermCriteria, int D, BeamBlurRegion BR>
void PhotonBeams<TermCriteria,D,BR>::preprocess()
{
	if (frozen) return;
	//Clear trees
	m_global_tree.clear();
	m_caustics_tree.clear();
	
	//Trace photons
	std::list<Photon<D> > global_photons, caustic_photons;
	std::vector<Beam<D,BR> > global_beams, caustic_beams;
	global_beams.clear();
	caustic_beams.clear();
	global_photons.clear();
	caustic_photons.clear();
	this->trace_beams(global_beams, caustic_beams, global_photons, caustic_photons);
	
	m_global_vector = global_beams;
	m_caustics_vector = caustic_beams;
	
	//Store photons into KDTrees and balance KDTrees
	typename std::list<Photon<D> >::iterator it;
	//global_photons.clear();
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
	//BeamTracing<TermCriteria, D, BR>::store(global_beams, "../../../ws/global_beams.txt");
	//BeamTracing<TermCriteria, D, BR>::store(caustic_beams, "../../../ws/caustic_beams.txt");
	if (global_beams.size()) {
		
		m_global_bvh = new BVH<Beam<D,BR>, D>(global_beams);
		//global_beams.clear();
	}
	
	if (caustic_beams.size()) {
		m_caustics_bvh = new BVH<Beam<D,BR>, D>(caustic_beams);
		//caustic_beams.clear();
	}
	
	if( global_photons.size() )
	{
		for( it = global_photons.begin(); it != global_photons.end(); it ++)
		{	
		 m_global_tree.store(std::vector<Real>(it->m_position.m_data, 
												  it->m_position.m_data+D), *it);
		}
		global_photons.clear();
		m_global_tree.balance();
	}
	if( caustic_photons.size() )
	{
		for( it = caustic_photons.begin(); it != caustic_photons.end(); it ++)
		{		
			m_caustics_tree.store(std::vector<Real>(it->m_position.m_data, 
													it->m_position.m_data+D), *it);
		}
		caustic_photons.clear();
		m_caustics_tree.balance();
	}
}

template<class TermCriteria, int D, BeamBlurRegion BR>
void PhotonBeams<TermCriteria,D,BR>::load(const std::string &namefile)
{
	FILE *f = fopen(namefile.c_str(), "r");
	if (!f)
	{
		fprintf(stderr, "Failed loading beams from file: File %s not found.\n", namefile.c_str());
		exit(0);
	}
	VectorN<D> o, w;
	Real d, r, power, time;
	int level;
	
	std::vector<Beam<D,BR> > global_beams;
	std::list<Photon<D> > global_photons;
	
	if (D == 3) 
	{	
		while (!feof(f))
		{
			fscanf(f, "[%f %f %f] [%f %f %f] [%f %f] power: %f level: %d time: %f\n", 
				   &o[0], &o[1], &o[2], &w[0], &w[1], &w[2], &r, &d, &power, &level, &time);
			Photon<D> p (o, w, power, time, level);
			Beam<D,BR> b (p, r, d, true, level, DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
			global_beams.push_back(b);
			global_photons.push_back(p);
			
		}
		m_inv_global_region_measure = 1./(M_PI*r*r);
	}
	else if (D == 2)
	{	
		while (!feof(f))
		{
			fscanf(f, "[%f %f] [%f %f] [%f %f] power: %f level: %d time: %f\n", 
				   &o[0], &o[1], &w[0], &w[1], &r, &d, &power, &level, &time);
			Photon<D> p (o, w, power, time, level);
			Beam<D,BR> b (p, r, d, true, level, DEFAULT_REFRACTION_INDEX, this->m_world->get_medium());
			global_beams.push_back(b);
			global_photons.push_back(p);
		}
		m_inv_global_region_measure = 0.5/r;
	}
	fclose(f);
	
	m_global_bvh = new BVH<Beam<D,BR>, D>(global_beams);
	global_beams.clear();
	
	typename std::list<Photon<D> >::iterator it;
	if( global_photons.size() )
	{
		for( it = global_photons.begin(); it != global_photons.end(); it ++)
		{	m_global_tree.store(std::vector<Real>(it->m_position.m_data, 
												  it->m_position.m_data+D), *it);
		}
		global_photons.clear();
		m_global_tree.balance();
	}
}


template<class TermCriteria, int D, BeamBlurRegion BR>
void PhotonBeams<TermCriteria,D,BR>::store(const std::string &namefile)const
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


// Integration functions
template<class TermCriteria, int D, BeamBlurRegion BR>
Spectrum PhotonBeams<TermCriteria,D,BR>::operator()(const Ray<D> &r) const
{
	Spectrum global_beams_contribution(0.), caustic_beams_contribution(0.), single_scattering(0.);
	std::vector<const Beam<D,BR>*> intersected_global_beams, intersected_caustics_beams;
	std::vector<const Beam<D,BR>*> v_intersected_global_beams, v_intersected_caustics_beams;
	std::vector<Beam<D,BR> > beams_to_write;
	
	static Real bvh_total_secs = 0., linear_total_secs = 0.;
	//static Real total_secs = 0.;
	static int count = 0;
	Timer timer;
	
	//Get the contribution of all the beams until the first intersection
	typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;
	typename std::vector<Beam<D,BR> >::const_iterator beam_it;
	
	typename std::vector<Intersection<D>*>::iterator ptr_it_it;	
	
	if ( rm_single_scattering && this->is_level(Scattering::SINGLE) ) {
		//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
		//Should be changed for a more general implementation or move it to another volume integrator
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
				 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
		
		if (this->m_world->get_medium())
		{
			Real stepSize = 0.001, step, dist = min(r.get_parameter(),100.0f), t0 = 0., pdf_light1;
			unsigned int nSamples = ceil(dist/stepSize), s; //max boundary of scene to 100.
			step = dist/nSamples;
			VectorN<D> p = r.get_origin(), pPrev;
			Spectrum Tr(1.);
			t0 += 0.01*stepSize;
			for (s = 0; s < nSamples; s++) {
				pPrev = p;
				p = r.get_origin() + r.get_direction()*t0;
				const LightSource<D> *l = this->m_world->sample_light(pdf_light1);
				LightSample<D> lsample;
				Real pdf_light2;
				Tr = exp(-u_t*t0);
				if(l->sample(p, lsample, pdf_light2)) 
				{
					single_scattering += Tr*this->m_world->get_medium()->f(p,lsample.dir,-r.get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2);
				}
				t0 += step; 
			}
			single_scattering *= step*u_s;
		}
	}
	std::vector<Intersection<D>*> intersections;
	//GLOBAL BEAMS
	if (m_global_bvh)
	{				
		//Intersect the beams directly
		m_global_bvh->intersect_all(r, r.get_parameter(), intersected_global_beams);
		
		for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
		{
			//if ((*ptr_beam_it)->get_level() != 0) 
			global_beams_contribution += (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())
			*(*ptr_beam_it)->contribution(r)*(*ptr_beam_it)->power_at(0.);
			//beams_to_write.push_back(**ptr_beam_it);
		}
		/*
		m_global_bvh->intersect_all(r, r.get_parameter(), intersections);
		for(ptr_it_it = intersections.begin(); ptr_it_it != intersections.end(); ptr_it_it++)
		{
			//if ((*ptr_beam_it)->get_level() != 0) 
			const Beam<D,BR> *b = (*ptr_it_it)->b;
			global_beams_contribution += b->get_medium()->get_scattering(b->get_origin())*b->contribution(*ptr_it_it)*b->power_at(0.);
		}*/
	}
	intersections.clear();
	//CAUSTICS
	if (m_caustics_bvh)
	{
		m_caustics_bvh->intersect_all(r, r.get_parameter(), intersected_caustics_beams);
		for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
			caustic_beams_contribution += (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())
											*(*ptr_beam_it)->contribution(r)*(*ptr_beam_it)->power_at(0.);		 
		
		/*
		 m_caustics_bvh->intersect_all(r, r.get_parameter(), intersections);
		for(ptr_it_it = intersections.begin(); ptr_it_it != intersections.end(); ptr_it_it++)
		{
			//if ((*ptr_beam_it)->get_level() != 0) 
			const Beam<D,BR> *b = (*ptr_it_it)->b;
			caustic_beams_contribution += b->get_medium()->get_scattering(b->get_origin())*b->contribution(*ptr_it_it)*b->power_at(0.);
		}*/
	}
	intersections.clear();
	intersected_global_beams.clear();
	intersected_caustics_beams.clear();
	return global_beams_contribution*m_inv_global_region_measure + caustic_beams_contribution*m_inv_caustics_region_measure + single_scattering;
}

template<class TermCriteria, int D, BeamBlurRegion BR>
Spectrum PhotonBeams<TermCriteria,D,BR>::operator()(const Intersection<D> &it) const
{
	Spectrum reflected_radiance(0.), attenuation(1.), global_beams_contribution(0.), caustic_beams_contribution(0.);
	Spectrum total_single_scattering(0.), single_scattering(0.);
	std::vector<const Beam<D,BR>*> intersected_global_beams, intersected_caustics_beams;
	
	//Get the contribution of all the beams until the first intersection
	typename std::vector<Beam<D,BR> >::const_iterator beam_it;
	typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;

	Spectrum Tr(1.); //accumulated extinction along the entire ray path (only for photons, not beams)
	
	if (rm_single_scattering && this->is_level(Scattering::SINGLE) ){
		//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
		//Should be changed for a more general implementation or move it to another volume integrator
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
				 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
		
		if(this->m_world->get_medium())
		{	
			Real step, stepSize = 0.001, dist = min(it.get_ray().get_parameter(),100.0f), t0 = 0., pdf_light1, pdf_light2;
			unsigned int nSamples = ceil(dist/stepSize), s; //max boundary of scene to 100.
			step = dist/nSamples;
			VectorN<D> p = it.get_ray().get_origin(), pPrev;
			Spectrum Tr0(1.);
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
					single_scattering += Tr0*this->m_world->get_medium()->f(p,lsample.dir,-it.get_ray().get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2);
				}
				t0 += step; 
			}
			single_scattering *= step*u_s;
		}
		total_single_scattering += single_scattering;
	}

	//GLOBAL
	if (m_global_bvh)
	{
		m_global_bvh->intersect_all(it.get_ray(), it.get_ray().get_parameter(), intersected_global_beams);
		for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
			global_beams_contribution += (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())
											*(*ptr_beam_it)->contribution(it.get_ray())*(*ptr_beam_it)->power_at(0.);
		
	}
	//CAUSTICS
	if (m_caustics_bvh)
	{
		m_caustics_bvh->intersect_all(it.get_ray(), it.get_ray().get_parameter(), intersected_caustics_beams);
		for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
			caustic_beams_contribution += (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())
											*(*ptr_beam_it)->contribution(it.get_ray())*(*ptr_beam_it)->power_at(0.);
		
	}
	intersected_global_beams.clear();
	intersected_caustics_beams.clear();
	
	if (!it.did_hit())
		return global_beams_contribution*m_inv_global_region_measure + caustic_beams_contribution*m_inv_caustics_region_measure + total_single_scattering;
	
	Medium<D> *m = this->m_world->get_medium(); //assuming infinite homogeneous medium
	if(m) 
	{
		Tr = exp(-m->get_extinction(it.get_ray().get_origin())*it.get_ray().get_parameter());
		//printf("Tr = %f %f %f\n", Tr[0], Tr[1], Tr[2]);
	}
	
	//Handle delta materials (e.g. perfect reflection/refraction)
	Intersection<D> shaded_it(it);
	while(shaded_it.material()->is_type(Reflectance::DELTA) )
	{
		Real pdf;
		//Check termination condition...
		if ( m_termCriteria(shaded_it, pdf))
			return global_beams_contribution*m_inv_global_region_measure + caustic_beams_contribution*m_inv_caustics_region_measure + total_single_scattering;
		
		//Get output direction, plus the attenuation
		VectorN<D> new_omega_i;
		Ray<D> new_ray;
		attenuation *= 1. - shaded_it.material()->get_absorption(it.get_uv());
		shaded_it.material()->sample_outgoing_ray(shaded_it, new_ray, pdf);
		
		//attenuation *= dot_abs(shaded_it.get_normal(), new_ray.get_direction())/pdf;
		
		Intersection<D> new_it;
		//Get new intersection
		this->m_world->first_intersection(new_ray, new_it);
		bool outside = dot(new_ray.get_direction(), shaded_it.get_normal()) > 0.;
		
		if (rm_single_scattering && outside && this->is_level(Scattering::SINGLE) ){
			
			Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
					 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
			
			if(this->m_world->get_medium())
			{
				//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
				//Should be changed for a more general implementation or move it to another volume integrator
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
					Tr0 = exp(-u_t*t0);
					if(l->sample(p, lsample, pdf_light2)) 
					{
						single_scattering += Tr0*this->m_world->get_medium()->f(p,lsample.dir,-new_ray.get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2);
					}
					t0 += step; 
				}
				single_scattering *= step*u_s;
				total_single_scattering += Tr*attenuation*single_scattering;
			}
		}
		
		//Now we know where the new ray hits => Get the contribution of all the beams until the intersection 
		if (m_global_bvh && outside)
		{
			m_global_bvh->intersect_all(new_ray, new_ray.get_parameter(), intersected_global_beams);
			for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
				global_beams_contribution += Tr*attenuation *(*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())
				*(*ptr_beam_it)->contribution(new_ray)*(*ptr_beam_it)->power_at(0.);
		}
		if (m_caustics_bvh && outside)
		{	
			m_caustics_bvh->intersect_all(new_ray, new_ray.get_parameter(), intersected_caustics_beams);
			for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
				caustic_beams_contribution += Tr*attenuation*(*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())
				*(*ptr_beam_it)->contribution(new_ray)*(*ptr_beam_it)->power_at(0.);
		}
		intersected_global_beams.clear();
		intersected_caustics_beams.clear();
		
		
		if( !new_it.did_hit() )
			return global_beams_contribution*m_inv_global_region_measure + caustic_beams_contribution*m_inv_caustics_region_measure + total_single_scattering;
		
		if (outside) //ray outside the object
			Tr *= exp(-m->get_extinction(new_ray.get_origin())*new_ray.get_parameter());
		
		shaded_it = new_it;
	}
	
	//REACHED A NON-DELTA SURFACE -> Compute radiance on surface
	if (rm_single_scattering && shaded_it.did_hit() && this->is_level(Scattering::SINGLE) )
	{
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
				 u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
		
		//DIRECT LIGHT ON SURFACE
		if (shaded_it.did_hit()){
			Real pdf_light1, pdf_light2;
			const LightSource<D> *l = this->m_world->sample_light(shaded_it, pdf_light1);
			LightSample<D> light_sample;
			if( l->sample(shaded_it, light_sample, pdf_light2 ) )
			{
				Spectrum radiance = light_sample.irradiance*shaded_it.material()->f(light_sample.dir,-shaded_it.get_ray().get_direction(),shaded_it.get_normal(),shaded_it.get_uv())
				*dot_clamped(-light_sample.dir,shaded_it.get_normal())*exp(-u_t*light_sample.dist)/(pdf_light1*pdf_light2);
				
				total_single_scattering += Tr*attenuation*radiance;
			}
		}
	}
	
	//	Get all samples where the radiance estimation is needed.
	//	If we are performing final gathering, compute the reflected
	//	samples. Else, just add the current intersection.
	list<pair<Spectrum, Intersection<D> > > samples;
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
			
			if( new_it.did_hit() )
				samples.push_back(pair<Spectrum,Intersection<D> >(form_factor,new_it));
		}
	}
	else
		samples.push_back(pair<Spectrum,Intersection<D> >(Spectrum(1),shaded_it));
	
	
	//Compute radiance estimation! 
	typename std::list<pair< Spectrum,Intersection<D> > >::iterator cit;
	for( cit = samples.begin(); cit != samples.end(); cit++ )
	{
		Spectrum sample_radiance(0.);
		VectorN<D> p = cit->second.get_position();
		VectorN<D> d = -cit->second.get_ray().get_direction();
		
		typename KDTree<Photon<D>,D>::Node blah;
		//std::vector<const typename KDTree<Photon<D>, D>::Node*> photons;KDTreeNode
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
		if (D == 2)
			reflected_radiance += sample_radiance*cit->first*(1./(2*max_distance));
		else // D == 3
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
		if (D == 2)
			reflected_radiance += sample_radiance*cit->first*(1./(2*max_distance));
		else // D == 3
			reflected_radiance += sample_radiance*cit->first*(1./(M_PI*max_distance*max_distance));
		
	}
	
	return Tr*reflected_radiance * attenuation * (m_final_gathering? this->m_inv_incoming_samples:1.) 
			+ global_beams_contribution*m_inv_global_region_measure + caustic_beams_contribution*m_inv_caustics_region_measure + total_single_scattering;
}


// Auxiliar structure to store samples where we
// will estimate the radiance
template<int D>
struct RadianceEstimationSample{
	Spectrum form_factor;
	Intersection<D> it;
	Real time;
	int bounce;
	RadianceEstimationSample(const Spectrum &ff, const Intersection<D> &i, Real t, int b=0)
	:form_factor(ff), it(i), time(t), bounce(b){}
};

template<class TermCriteria, int D, BeamBlurRegion BR>
void PhotonBeams<TermCriteria,D,BR>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	
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
				time_blur(sample, samples, delta_time/nb_time_samples);
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
				time_blur(sample, samples, delta_time/nb_time_samples);
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
					time_blur(sample, samples, delta_time/nb_time_samples);
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
					time_blur(sample, samples, delta_time/nb_time_samples);
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
		// Just use max distance...
		int tm_kernel = _TEMPORAL_KERNEL_SIZE_>1?_TEMPORAL_KERNEL_SIZE_:1;
		float inv_tm_kernel = 1./static_cast<Real>(tm_kernel);
		float init_tm_kernel = tm_kernel>1?(-tm_kernel/2.):0.,
			  c_tm_kernel;

		//	... and compute their contribution
		for( ph = photons.begin(); ph != photons.end(); ph++ )
			for( int ts=0, c_tm_kernel = init_tm_kernel; ts < tm_kernel; ++ts, c_tm_kernel += delta_time )
			{
				samples.push_back(RadianceSample( Tr*cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv())
											* ((*ph)->data().m_power)
											* dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal())
											* fc_global * inv_tm_kernel,
											cit->time + (*ph)->data().m_time + c_tm_kernel, 
											cit->bounce + (*ph)->data().m_level,
											cit->time));

			}
		/*
			for( ph = photons.begin(); ph != photons.end(); ph++ )
				samples.push_back(RadianceSample( Tr*cit->it.material()->f((*ph)->data().m_direction, d, cit->it.get_normal(), cit->it.get_uv())
												 * ((*ph)->data().m_power)
												 * dot_clamped(-(*ph)->data().m_direction, cit->it.get_normal())
												 * fc_global,
												 cit->time + (*ph)->data().m_time, cit->bounce + (*ph)->data().m_level, cit->time)); 
		*/
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
}

// Integration functions
template<class TermCriteria, int D, BeamBlurRegion BR>
void PhotonBeams<TermCriteria,D,BR>::operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	Spectrum global_beams_contribution(0.), caustic_beams_contribution(0.);
	std::vector<const Beam<D,BR>*> intersected_global_beams, intersected_caustics_beams;
	std::vector<const Beam<D,BR>*> v_intersected_global_beams, v_intersected_caustics_beams;
	std::vector<Beam<D,BR> > beams_to_write;
	
	static Real bvh_total_secs = 0., linear_total_secs = 0.;
	//static Real total_secs = 0.;
	static int count = 0;
	Timer timer;
	
	//Get the contribution of all the beams until the first intersection
	typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;
	typename std::vector<Beam<D,BR> >::const_iterator beam_it;
	
	if (rm_single_scattering && this->is_level(Scattering::SINGLE) ) {
		//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
		//Should be changed for a more general implementation or move it to another volume integrator
		Spectrum u_t=(this->m_world->get_medium()?this->m_world->get_medium()->get_extinction(VectorN<D>(0.)):Spectrum(0.)), 
		u_s=(this->m_world->get_medium()?this->m_world->get_medium()->get_scattering(VectorN<D>(0.)):Spectrum(0.));
		Real ior = r.get_refraction_index();
		if (this->m_world->get_medium())
		{
			Real stepSize = 0.001, step, dist = min(r.get_parameter(),100.0f), t0 = 0., pdf_light1;
			unsigned int nSamples = ceil(dist/stepSize), s; //max boundary of scene to 100.
			step = dist/nSamples;
			VectorN<D> p = r.get_origin(), pPrev;
			Spectrum Tr(1.);
			t0 += 0.01*stepSize;
			for (s = 0; s < nSamples; s++) {
				pPrev = p;
				p = r.get_origin() + r.get_direction()*t0;
				const LightSource<D> *l = this->m_world->sample_light(pdf_light1);
				LightSample<D> lsample;
				Real pdf_light2;
				Tr *= exp(-this->m_world->get_medium()->get_extinction(VectorN<D>(0.))*step);
				if(l->sample(p, lsample, pdf_light2)) 
				{
					samples.push_back(RadianceSample(step*u_s*Tr*this->m_world->get_medium()->f(p,lsample.dir,-r.get_direction())*lsample.irradiance*exp(-u_t*lsample.dist)/(pdf_light1*pdf_light2),
													 this->m_world->time_of_flight(lsample.dist + t0)*ior, r.get_level(), this->m_world->time_of_flight(t0)*ior));
				}
				t0 += step; 
			}
		}
	}
	
	//GLOBAL BEAMS
	if (m_global_bvh)
	{	
		m_global_bvh->intersect_all(r, r.get_parameter(), intersected_global_beams);
		
		for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
		{
			RadianceSample sample;
			(*ptr_beam_it)->contribution(r, sample);
			sample.radiance *= (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_global_region_measure;
			if (BR != BR2D)
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, film->get_exposure_time());
		}
	}
	//CAUSTICS
	if (m_caustics_bvh)
	{
		m_caustics_bvh->intersect_all(r, r.get_parameter(), intersected_caustics_beams);
		for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
		{
			RadianceSample sample;
			(*ptr_beam_it)->contribution(r, sample);
			sample.radiance *= (*ptr_beam_it)->get_medium()->get_scattering((*ptr_beam_it)->get_origin())*(*ptr_beam_it)->power_at(0.)*m_inv_caustics_region_measure;
            if (BR != BR2D) 
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, film->get_exposure_time());
		}
	}
	intersected_global_beams.clear();
	intersected_caustics_beams.clear();
}

// Integration functions
template<class TermCriteria, int D, BeamBlurRegion BR>
void PhotonBeams<TermCriteria,D,BR>::operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples ) const
{
	//FUNCTION FOR INTEGRATING THE INCOMING RADIANCE IN A POINT OF THE SCENE
	Spectrum global_beams_contribution(0.), caustic_beams_contribution(0.);
	std::vector<const Beam<D,BR>*> intersected_global_beams, intersected_caustics_beams;
	std::vector<const Beam<D,BR>*> v_intersected_global_beams, v_intersected_caustics_beams;
	std::vector<Beam<D,BR> > beams_to_write;
	
	static Real bvh_total_secs = 0., linear_total_secs = 0.;
	//static Real total_secs = 0.;
	static int count = 0;
	Timer timer;
	
	//Get the contribution of all the beams until the first intersection
	typename std::vector<const Beam<D,BR>*>::iterator ptr_beam_it;
	typename std::vector<Beam<D,BR> >::const_iterator beam_it;
	
	if (rm_single_scattering && this->is_level(Scattering::SINGLE) ) {
		//COMPUTE SINGLE SCATTERING, ASSUMING HOMOGENEOUS ISOTROPIC MEDIUM IN THE WHOLE SCENE
		//Should be changed for a more general implementation or move it to another volume integrator
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
	
	//GLOBAL BEAMS
	if (m_global_bvh)
	{	
		m_global_bvh->intersect_all(p, intersected_global_beams);
		
		for(ptr_beam_it = intersected_global_beams.begin(); ptr_beam_it != intersected_global_beams.end(); ptr_beam_it++)
		{
			RadianceSample sample;
			
			//POINT QUERY x BEAM DATA
			(*ptr_beam_it)->contribution(p, sample);
			sample.radiance *= (*ptr_beam_it)->power_at(0.)*m_inv_global_region_measure;
			if (BR == BR2D)
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, delta_time/nb_time_samples);
		}
	}
	//CAUSTICS
	if (m_caustics_bvh)
	{
		m_caustics_bvh->intersect_all(p, intersected_caustics_beams);
		for(ptr_beam_it = intersected_caustics_beams.begin(); ptr_beam_it != intersected_caustics_beams.end(); ptr_beam_it++)
		{
			RadianceSample sample;
			(*ptr_beam_it)->contribution(p, sample);
			sample.radiance *= (*ptr_beam_it)->power_at(0.)*m_inv_caustics_region_measure;
			if (BR == BR2D) 
				samples.push_back(sample);
			else if (sample.radiance.avg() != 0.)
				time_blur(sample, samples, delta_time/nb_time_samples);
		}
	}
	intersected_global_beams.clear();
	intersected_caustics_beams.clear();
}

#endif //_PHOTON_BEAMS_H_